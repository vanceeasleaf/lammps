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

#include "comm.h"
#include "math.h"
#include "string.h"
#include "stdio.h"
#include "stdlib.h"
#include "compute_tc.h"
#include "atom.h"
#include "input.h"
#include "variable.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INVOKED_SCALAR 1
#define INVOKED_VECTOR 2

/* ---------------------------------------------------------------------- */

ComputeTC::ComputeTC(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 11 && narg != 12 && narg != 13) error->all(FLERR,"Illegal compute tc command");

  scalar_flag = 1;
  extscalar = 0; 

  int n;
  if (strncmp(arg[3],"c_",2) == 0) {
    n = strlen(arg[3]);
    id_temp = new char[n];
    strcpy(id_temp,&arg[3][2]);
  } else {
    error->all(FLERR,"Could not find c_ for compute temp ID in compute tc command");
  }  

  if (strncmp(arg[4],"c_",2) == 0) {
    n = strlen(arg[4]);
    id_flux = new char[n];
    strcpy(id_flux,&arg[4][2]);
  } else {
    error->all(FLERR,"Could not find c_ for compute heat/flux ID in compute tc command");
  }  

  if (strncmp(arg[5],"v_",2) == 0) {
    n = strlen(arg[5]);
    id_variable1 = new char[n];
    strcpy(id_variable1,&arg[5][2]);
  } else {
    error->all(FLERR,"Could not find v_ for variable1 ID in compute tc command");
  }

  if (strncmp(arg[6],"v_",2) == 0) {
    n = strlen(arg[6]);
    id_variable2 = new char[n];
    strcpy(id_variable2,&arg[6][2]);
  } else {
    error->all(FLERR,"Could not find v_ for variable2 ID in compute tc command");
  }

  n = strlen(arg[7]) + 1;
  direction = new char[n];
  strcpy(direction,arg[7]);
  if (strcmp(direction,"x") != 0 && strcmp(direction,"y") != 0 && strcmp(direction,"z") != 0 && strcmp(direction,"iso") != 0) error->all(FLERR,"Wrong direction value for compute tc command");

  n = strlen(arg[8]) + 1;
  portion = new char[n];
  strcpy(portion,arg[8]);
  if (strcmp(portion,"first") != 0 && strcmp(portion,"second") != 0) error->all(FLERR,"Wrong heat flux portion value for compute tc command");

  mbig = atoi(arg[9]);
  if (mbig <= 0) error->all(FLERR,"Wrong mbig value for compute tc command"); 

  tc_ntimestep_start = atoi(arg[10]);
  if (tc_ntimestep_start < 0) error->all(FLERR,"Wrong tc_ntimestep_start value for compute tc command");

  if(narg >= 12) tc_ntimestep_frequency = atoi(arg[11]);
  else tc_ntimestep_frequency = 100000;
  if(narg >= 13) tc_interval = atoi(arg[12]);
  else tc_interval = 1;
  int itemp = modify->find_compute(id_temp);
  if (itemp < 0)
    error->all(FLERR,"Could not find compute temp ID for compute tc command");
  if (strcmp(modify->compute[itemp]->style,"temp") != 0)
    error->all(FLERR,"Compute tc compute ID does not compute temp");

  int iflux = modify->find_compute(id_flux);
  if (iflux < 0)
    error->all(FLERR,"Could not find compute heat/flux ID for compute tc command");
  if (strcmp(modify->compute[iflux]->style,"heat/flux") != 0)
    error->all(FLERR,"Compute tc compute ID does not compute heat/flux");

  int ivariable1 = input->variable->find(id_variable1);
  if (ivariable1 < 0) error->all(FLERR,"Could not find variable1 ID for compute tc command");
  if (input->variable->equalstyle(ivariable1) == 0)
    error->all(FLERR,"Compute tc variable1 is not equal-style variable");

  int ivariable2 = input->variable->find(id_variable2);
  if (ivariable2 < 0) error->all(FLERR,"Could not find variable2 ID for compute tc command");
  if (input->variable->equalstyle(ivariable2) == 0)
    error->all(FLERR,"Compute tc variable2 is not equal-style variable");

  ac = NULL;
  acc = NULL;
  tcc = NULL;
  jflux = NULL;

  thermal_conductivity = 0.0;
}

/* ---------------------------------------------------------------------- */

ComputeTC::~ComputeTC()
{ 
  delete [] id_temp;
  delete [] id_flux;
  delete [] id_variable1;
  delete [] id_variable2;
  delete [] direction;
  delete [] portion;
  memory->sfree(ac);
  memory->sfree(acc);
  memory->sfree(tcc);
  //memory->destroy_2d_double_array(jflux);
      for(int i=0;i<mbig+1;i++)
    delete jflux[i];
	delete [] jflux;
}

/* ---------------------------------------------------------------------- */

void ComputeTC::init()
{
  // error checks

  int itemp = modify->find_compute(id_temp);
  if (itemp < 0)
    error->all(FLERR,"Could not find compute temp ID for compute tc command");
  c_temp = modify->compute[itemp];

  int iflux = modify->find_compute(id_flux);
  if (iflux < 0)
    error->all(FLERR,"Could not find compute heat/flux ID for compute tc command");
  c_flux = modify->compute[iflux];

  int ivariable1 = input->variable->find(id_variable1);
  if (ivariable1 < 0) error->all(FLERR,"Could not find variable1 ID for compute tc command");

  int ivariable2 = input->variable->find(id_variable2);
  if (ivariable2 < 0) error->all(FLERR,"Could not find variable2 ID for compute tc command");

}

/* ---------------------------------------------------------------------- */

double ComputeTC::compute_scalar()
{
  invoked_scalar = update->ntimestep;

  int ntimestep = update->ntimestep;
  int mtimestep = ntimestep/tc_interval;
  int nsteps = update->nsteps;
  if (mtimestep == 0) {
    memory->sfree(ac);
    memory->sfree(acc);
    memory->sfree(tcc);
    //memory->destroy_2d_double_array(jflux);
    
    //for(int i=0;i<mbig+1;i++)
    //delete jflux[i];
	//delete [] jflux;
    ac  = (double *) memory->smalloc((mbig+1)*sizeof(double),"compute tc:ac");
    acc = (double *) memory->smalloc((mbig+1)*sizeof(double),"compute tc:acc");
    tcc = (double *) memory->smalloc((mbig+1)*sizeof(double),"compute tc:tcc");
    //jflux = memory->create_2d_double_array(mbig+1,6,"compute tc:jflux");
    jflux=new double*[mbig+1];
    for(int i=0;i<mbig+1;i++)
    jflux[i]=new double[6];
    for (int m = 0; m <= mbig; m++) {
      ac[m]  = 0.0;
      acc[m] = 0.0;
      tcc[m] = 0.0;
    }
    ac0 = 0.0;
  }

  // invoke 2 computes if they haven't been already

  if (!(c_temp->invoked_flag & INVOKED_SCALAR)) {
    c_temp->compute_scalar();
    c_temp->invoked_flag |= INVOKED_SCALAR;
  }
  if (!(c_flux->invoked_flag & INVOKED_VECTOR)) {
    c_flux->compute_vector();
    c_flux->invoked_flag |= INVOKED_VECTOR;
  }

  double temp = c_temp->scalar;
  double *jt  = c_flux->vector;

  factor1 = input->variable->compute_equal(0);
  factor2 = input->variable->compute_equal(1);

  int s, t;
  double xtmp,ytmp,ztmp;

  if (strcmp(portion,"first") == 0) {
    if (strcmp(direction,"x") == 0) {
      xtmp = jt[0];
      ac0 += xtmp*xtmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][0] = xtmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][0]*xtmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][0] = xtmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][0]*xtmp;
        }
      }
    } else if (strcmp(direction,"y") == 0) {
      ytmp = jt[1];
      ac0 += ytmp*ytmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][1] = ytmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][1]*ytmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][1] = ytmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][1]*ytmp;
        }
      }
    } else if (strcmp(direction,"z") == 0) {
      ztmp = jt[2];
      ac0 += ztmp*ztmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][2] = ztmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][2]*ztmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][2] = ztmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][2]*ztmp;
        }
      }
    } else if (strcmp(direction,"iso") == 0) {
      xtmp = jt[0];
      ytmp = jt[1];
      ztmp = jt[2];
      ac0 += xtmp*xtmp + ytmp*ytmp + ztmp*ztmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][0] = xtmp;
        jflux[mtimestep][1] = ytmp;
        jflux[mtimestep][2] = ztmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][0]*xtmp + jflux[mtimestep-m][1]*ytmp + jflux[mtimestep-m][2]*ztmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][0] = xtmp;
        jflux[s][1] = ytmp;
        jflux[s][2] = ztmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][0]*xtmp + jflux[t][1]*ytmp + jflux[t][2]*ztmp;
        }
      }
    }
  } else if (strcmp(portion,"second") == 0) {
    if (strcmp(direction,"x") == 0) {
      xtmp = jt[3];
      ac0 += xtmp*xtmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][3] = xtmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][3]*xtmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][3] = xtmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][3]*xtmp;
        }
      }
    } else if (strcmp(direction,"y") == 0) {
      ytmp = jt[4];
      ac0 += ytmp*ytmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][4] = ytmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[4][mtimestep-m]*ytmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][4] = ytmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][4]*ytmp;
        }
      }
    } else if (strcmp(direction,"z") == 0) {
      ztmp = jt[5];
      ac0 += ztmp*ztmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][5] = ztmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][5]*ztmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][5] = ztmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][5]*ztmp;
        }
      }
    } else if (strcmp(direction,"iso") == 0) {
      xtmp = jt[3];
      ytmp = jt[4];
      ztmp = jt[5];
      ac0 += xtmp*xtmp + ytmp*ytmp + ztmp*ztmp;
      if (mtimestep <= mbig) {
        jflux[mtimestep][3] = xtmp;
        jflux[mtimestep][4] = ytmp;
        jflux[mtimestep][5] = ztmp;
        for (int m = 1; m <= mtimestep; m++) {    
          ac[m] += jflux[mtimestep-m][3]*xtmp + jflux[mtimestep-m][4]*ytmp + jflux[mtimestep-m][5]*ztmp;
        }
      } else {
        s = mtimestep%(mbig+1);
        jflux[s][3] = xtmp;
        jflux[s][4] = ytmp;
        jflux[s][5] = ztmp;
        for (int m = 1; m <= mbig; m++) {
          t = (mbig+1+s-m)%(mbig+1);
          ac[m] += jflux[t][3]*xtmp + jflux[t][4]*ytmp + jflux[t][5]*ztmp;
        }
      }
    }
  }

  for (int m = 1; m <= mbig; m++) {
    acc[m] = ac[m];
  }

  double volume = domain->xprd * domain->yprd * domain->zprd;
  if(strcmp(direction,"iso") == 0) {
    factor = update->dt*tc_interval/3.0/force->boltz/volume/temp/temp;
  } else {
    factor = update->dt*tc_interval/force->boltz/volume/temp/temp;
  }

  if (mtimestep <= mbig) {
    for (int m = 1; m <= mtimestep; m++) {
      acc[m] /= (mtimestep-m+1);
      if(m==1) tcc[m] = acc[m] * factor;
      if(m>1) tcc[m] = tcc[m-1] + acc[m] * factor;
    }
    thermal_conductivity = tcc[mtimestep];
  } else {
    for (int m = 1; m <= mbig; m++) {
      acc[m] /= (mtimestep-m+1);
      if(m==1) tcc[m] = acc[m] * factor;
      if(m>1) tcc[m] = tcc[m-1] + acc[m] * factor;
    }
    thermal_conductivity = tcc[mbig];
  }

  if (comm->me == 0) {
    if(ntimestep>=tc_ntimestep_start && ntimestep%tc_ntimestep_frequency==0 || ntimestep==nsteps) {
      FILE *fp;
      char *str1 = new char[10];
      char *space = new char[2];
      sprintf(str1, "%s%d%s","ac",ntimestep,".dat");    
      sprintf(space, "\t"); 
      int n = strlen(str1) + 1;
      char *str = new char[n];    
      strcpy(str,str1);
      fp = fopen(str,"w");
      if (fp == NULL) error->one(FLERR,"Cannot open ac.dat");
      for (int m = 1; m <= mbig; m++) { 
        fprintf(fp,"%f%s%20.10f\n",m*update->dt*tc_interval,space,acc[m]/(ac0/(mtimestep+1))*factor1);
      }
      fclose(fp);

      sprintf(str1, "%s%d%s","tc",ntimestep,".dat"); 
      strcpy(str,str1);
      fp = fopen(str,"w");
      if (fp == NULL) error->one(FLERR,"Cannot open tc.dat");
      for (int m = 1; m <= mbig; m++) { 
        fprintf(fp,"%f%s%20.10f\n",m*update->dt*tc_interval,space,tcc[m]*factor2);
      }
      fclose(fp);
    }
  }
  scalar = thermal_conductivity*factor2;
  return scalar;
}
