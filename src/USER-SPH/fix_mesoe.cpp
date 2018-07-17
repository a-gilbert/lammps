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
#include "math_const.h"
#include "fdint.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace LAMMPS_NS::FDINT;

/* ---------------------------------------------------------------------- */

FixMesoe::FixMesoe(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg) {

  if ((atom->e_flag != 1) || (atom->rho_flag != 1) || (atom->smu_flag != 1)
      || (atom->cv_flag != 1))
    error->all(FLERR,
        "fix mesoe command requires atom_style with energy, density, chemical potential, and heat capacity.");

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


/*--------------------------------------------
  The function that returns zero when we've reached
  the correct heat capacity
  -----------------------------------------*/
FixMesoe::f(double cv, double rho, double e, double nkb, double alpha)
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
Utilizing a false position method, calculate the electron
heat capacity in a given unit system.
 --------------------------------------------------------*/
double FixMesoe::solve_cv(double rho, double e, double mass)
{
  double out, alpha, nkb, cv_u, cv_l, fu, fl, cv_r, fr, cond;

  alpha = Force->hplanck/MY_2PI;
  alpha = alpha*alpha;
  alpha = 0.5*alpha/(Force->e_mass * Force->boltz);
  alpha = pow(alpha, 1.5);
  alpha = 2*MY_PI*MY_PI*alpha;

  cv_u = 10;
  cv_l = 1e-10;
  fu = f(cv_u, rho, e, nkb, alpha);
  fl = f(cv_l, rho, e, nkb, alpha);
  if (fu < 0) {
    while (fu < 0) {
      cv_u = cv_u*cv_u;
      fu = f(cv_u, rho, e, nkb, alpha);
    }
  }
  if (fl > 0) {
    while (fl > 0) {
      cv_l = cv_l*cv_l;
      fl = f(cv_l, rho, e, nkb, alpha);
    }
  }

  cond = 1
  for(int i = 1, i < 20 && cond != 0, i++) {
    cv_r = cv_u*fl - cv_l*fu;
    cv_r = xr/(fl - fu);
    fr = f(cv_r, rho, e, nkb, alpha);
    cond = fl*fr;
    if (cond < 0) {
      cv_u = cv_r;
      fu = fr;
    } else if (cond > 0) {
      cv_l = cv_r;
      fl = fr;
    } else if (cond == 0) {
      cv_l = cv_r;
      fl = fr;
      cv_u = cv_r;
      fu = fr;
    }
  }
  out = 0.5*(cv_u + cv_l);
  return out;
}

//Get the reduced chemical potential mu/kb*T
double FixMesoe::get_mu(double rho, double e, double cv) {
  double mu, alpha, x;

  alpha = Force->hplanck/MY_2PI;
  alpha = alpha*alpha;
  alpha = 0.5*alpha/(Force->e_mass * Force->boltz);
  alpha = pow(alpha, 1.5);
  alpha = 2*MY_PI*MY_PI*alpha;

  x = alpha * rho * pow(cv/e, 1.5);
  mu = arcih(x);
  return mu;
}


void FixMesoe::setup_pre_force(int vflag)
{
  // set vest equal to v. Solve for heat capacity, set chemical potential.
  double **v = atom->v;
  double **vest = atom->vest;
  double *rho = atom->rho;
  double *e = atom->e;
  double *cv = atom->cv;
  double *smu = atom->smu;
  int *mask = atom->mask;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int rmass_flag = atom->rmass_flag;
  int nlocal = atom->nlocal;
  double pmass;
  if (igroup == atom->firstgroup)
    nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (rmass_flag) {
        pmass = rmass[i];
      } else {
        pmass = mass[type[i]];
      }
      vest[i][0] = v[i][0];
      vest[i][1] = v[i][1];
      vest[i][2] = v[i][2];
      cv[i] = solve_cv(rho[i], e[i], pmass);
      smu[i] = get_mu(0, rho[i], e[i], cv[i]);
    }
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
