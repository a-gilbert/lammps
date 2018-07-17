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
#include "fix_meso.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "pair.h"
#include "fdint.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::FDINT;

/* ---------------------------------------------------------------------- */

FixMesoe::FixMesoe(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1))
    error->all(FLERR,
        "fix meso command requires atom_style with both energy and density");

  if (narg != 3)
    error->all(FLERR,"Illegal number of arguments for fix meso command");

  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

int FixMesoe::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMesoe::init() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}

void FixMesoe::setup_pre_force(int vflag)
{
  // set vest equal to v. Solve for heat capacity, chemical potential
  double **v = atom->v;
  double **vest = atom->vest;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];
    }
  }
}

FixMesoe::g(double cvn, double fxn, double rho, double e, double alpha,
            double nkb)
{
  double out = f(cvn + fxn, rho, e, alpha, nkb)/fxn;
  out = out - 1.0;
  return out;
}

/*--------------------------------------------
  The function that returns zero when we've reached
  the correct heat capacity
  -----------------------------------------*/
FixMesoe::f(double cv, double rho, double e, double alpha, double nkb)
{
  double x = alpha*rho*pow(cv/e, 1.5);
  double y = arcih(x);
  double out = 5.0*i3h(y)/(2.0*x);
  out = out - 4.5*x/inh(y);
  out = nkb*out;
  out = cv - out;
  return out;
}


/*-------------------------------------------------------
Utilizing a relaxed steffen step solver, calculate the electron
heat capacity in a given unit system.
 --------------------------------------------------------
 */
double FixMesoe::solve_cv(int units, double rho, double e, double mass)
{
  if (units == 0) {
    //SI
    double alpha = 1.835080110417837e-82;
    double nkb = mass * 1.515633314709168e7; //particle mass * kb/m_{electron} = N*k_{b}
    double w = 1.0;
    double cv = 1.0;
    for(int i = 1, i < 32, i++) {
      double fxn = f(cv, rho, e, alpha, nkb);
      double dx = fxn/g(cv, fxn, rho, e, alpha, nkb);
      if(dx >= )
    }
  } else if (units == 1) {
    //CGS
  } else {
    //???
  }

}

/* ----------------------------------------------------------------------
 allow for both per-type and per-atom mass
 ------------------------------------------------------------------------- */

void FixMesoe::initial_integrate(int vflag) {
  // update v and x and rho and e of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **vest = atom->vest;
  double *rho = atom->rho;
  double *drho = atom->drho;
  double *e = atom->e;
  double *de = atom->de;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int i;
  double dtfm;

  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }

      e[i] += dtf * de[i]; // half-step update of particle internal energy
      rho[i] += dtf * drho[i]; // ... and density

      // extrapolate velocity for use with velocity-dependent potentials, e.g. SPH
      vest[i][0] = v[i][0] + 2.0 * dtfm * f[i][0];
      vest[i][1] = v[i][1] + 2.0 * dtfm * f[i][1];
      vest[i][2] = v[i][2] + 2.0 * dtfm * f[i][2];

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoe::final_integrate() {

  // update v, rho, and e of atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double *e = atom->e;
  double *de = atom->de;
  double *rho = atom->rho;
  double *drho = atom->drho;
  int *type = atom->type;
  int *mask = atom->mask;
  double *mass = atom->mass;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;
  double dtfm;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      if (rmass_flag) {
        dtfm = dtf / rmass[i];
      } else {
        dtfm = dtf / mass[type[i]];
      }
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      e[i] += dtf * de[i];
      rho[i] += dtf * drho[i];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixMesoe::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
