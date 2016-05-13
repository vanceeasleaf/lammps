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

#ifdef COMPUTE_CLASS

ComputeStyle(tc,ComputeTC)

#else

#ifndef LMP_COMPUTE_TC_H
#define LMP_COMPUTE_TC_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeTC : public Compute {
 public:
  ComputeTC(class LAMMPS *, int, char **);
  ~ComputeTC();
  void init();
  double compute_scalar();

 private:
  char   *id_temp,*id_flux,*id_variable1,*id_variable2,*direction,*portion;
  class  Compute *c_temp,*c_flux;

  int    mbig,tc_ntimestep_start,tc_ntimestep_frequency,tc_interval;
  double factor,factor1,factor2,thermal_conductivity;
  double *ac,*acc,*tcc,**jflux,ac0;
};

}

#endif
#endif
