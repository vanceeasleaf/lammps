/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "bond_fpu.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondFPU::BondFPU(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondFPU::~BondFPU()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(alpha);
    memory->destroy(beta);
    memory->destroy(g);
  }
}

/* ---------------------------------------------------------------------- */

void BondFPU::compute(int eflag, int vflag)
{
  int i1,i2,n,type;
  double delx,dely,delz,ebond,fbond;
  double rsq,r,dr,rk;

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < nbondlist; n++) {
    i1 = bondlist[n][0];
    i2 = bondlist[n][1];
    type = bondlist[n][2];

    delx = x[i1][0] - x[i2][0];
    dely = x[i1][1] - x[i2][1];
    delz = x[i1][2] - x[i2][2];

    rsq = delx*delx + dely*dely + delz*delz;
    r = sqrt(rsq);
    dr = r - r0[type];
    rk = k[type] * dr + alpha[type]* dr*dr+ beta[type]* dr*dr*dr;

    double rp=sqrt(delx*delx+dely*+dely);
    double rv=fabs(delz);
    double c=rp-r0[type];
    double d=fabs(c);
    double e=1.0;
    if(d<1e-8)e=0.0;
    if(c<0.0)e=-1.0;
    double fx=-.25*g[type]*(delz*delz)*c*delx/rp;
    double fy=-.25*g[type]*(delz*delz)*c*dely/rp;
    double fz=-.25*g[type]*delz*c*c;
    // force & energy

    if (r > 0.0) fbond = -rk/r;
    else fbond = 0.0;

    if (eflag) ebond = (.5*k[type] * dr * dr+.33333333333*alpha[type]* dr*dr*dr+ .25*beta[type]* dr*dr*dr* dr+.25*g[type]*delz*delz*c*c);

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx*fbond+fx;
      f[i1][1] += dely*fbond+fy;
      f[i1][2] += delz*fbond+fz;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx*fbond+fx;
      f[i2][1] -= dely*fbond+fy;
      f[i2][2] -= delz*fbond+fz;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx,dely,delz);
  }
}

/* ---------------------------------------------------------------------- */

void BondFPU::allocate()
{
  allocated = 1;
  int n = atom->nbondtypes;

  memory->create(k,n+1,"bond:k");
  memory->create(r0,n+1,"bond:r0");
  memory->create(beta,n+1,"bond:alpha");
  memory->create(beta,n+1,"bond:beta");
  memory->create(g,n+1,"bond:g");
  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondFPU::coeff(int narg, char **arg)
{
  if (narg != 5) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double r0_one = force->numeric(FLERR,arg[2]);
  double alpha_one = force->numeric(FLERR,arg[3]);
  double beta_one = force->numeric(FLERR,arg[4]);
  double g_one = force->numeric(FLERR,arg[5]);
  rabs = force->numeric(FLERR,arg[6]);
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    alpha[i] = alpha_one;
    beta[i] = beta_one;
    g[i] = g_one;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondFPU::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondFPU::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&alpha[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&beta[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&g[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondFPU::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
    fread(&alpha[1],sizeof(double),atom->nbondtypes,fp);
    fread(&beta[1],sizeof(double),atom->nbondtypes,fp);
    fread(&g[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&alpha[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&beta[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&g[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondFPU::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,k[i],r0[i],alpha[i],beta[i],g[i]);
}

/* ---------------------------------------------------------------------- */

double BondFPU::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr +alpha[type]* dr*dr+ beta[type]* dr*dr*dr;
  fforce = 0;
  if (r > 0.0) fforce = -rk/r;
  return .5*k[type] * dr *dr+.33333333333*alpha[type]* dr*dr*dr+ .25*beta[type]* dr*dr*dr*dr;
}
