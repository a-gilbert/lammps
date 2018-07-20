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

#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include "fix_ltemp.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "pair.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLTemp::FixLTemp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->temp_flag != 1))
    error->all(FLERR,
        "fix ltemp command requires atom_style with temperature.");

  if (narg != 4)
    error->all(FLERR,"Illegal number of arguments for fix meso command");

  cut = atof(arg[3]);

}

/* ---------------------------------------------------------------------- */

int FixLTemp::setmask() {
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}


void FixLTemp::setup_pre_force(int vflag)
{
  // set local temp
  int i, j, k, itype, jtype;
  double r, vrms, v2, ri, imass;
  int n_neigh;
  double kb = force->boltz;
  double **v = atom->v;
  double **x = atom->x;
  double *temp = atom->temp;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;



  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vrms = 0;
      n_neigh = 0;
      itype = type[i];
      imass = mass[itype];
      for(j = 0; j < nlocal; j++) {
          if(mask[j] & groupbit) {
            r = 0;
            for(k = 0; k < 3; k++) {
              ri = (x[i][k] - x[j][k]);
              ri = ri*ri;
              r += ri;
            }
            r = sqrt(r);
            if(r < cut) {
              n_neigh++;
              v2 = 0;
              for(k = 0; k < 3; k++) {
                v2 += v[j][k]*v[j][k];
              }
              vrms += v2;
            }
          }
        }
      vrms = imass*vrms/(3*kb*n_neigh);
      temp[i] = vrms;
    }
  }
}
