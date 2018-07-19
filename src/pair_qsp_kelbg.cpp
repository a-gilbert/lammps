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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_qsp_kelbg.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairQspKelbg::PairQspKelbg(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairQspKelbg::~PairQspKelbg()
{
  if (!copymode) {
    if (allocated) {
      memory->destroy(setflag);

      memory->destroy(cut_kelbgsq);
      memory->destroy(on);
      memory->destroy(style);
    }
  }
}

double PairQspKelbg::j1(double x, double xi) {
  double out = exp(-1*x);
  out = pow(out, -1*x -1);
  double t1 = 1-exp(-1*MY_2PI*MY_PIS*xi/x);
  out = out/t1;
  return out;
}

double PairQspKelbg::get_j1(double xi) {
  static const double nodes[] = {0.6170308532782703957143,
                                 2.112965958578524151141,
                                 4.610833151017532413683,
                                 8.399066971204842190477,
                                 14.26010306592083084898};

  static const double weights[] = {0.348014540023348861432,
                                   0.50228067413249296034,
                                   0.1409159194944725554027,
                                   0.00871989302609998253045,
                                   6.89733235856402952577e-5};
  double out = 0;
  for (int i = 0; i < 5; i++) {
    out += weights[i]*j1(nodes[i], xi);
  }
  return out;
}

double PairQspKelbg::get_lsee(double xi) {
  double out = 2*MY_2PI*xi*get_j1(xi);
  return log(out);
}

/*----------------------------------------------------------------------- */
double PairQspKelbg::get_lsei(double) {
  return 1;
}

/* ---------------------------------------------------------------------- */

double PairQspKelbg::get_nu(int style, double xi) {
  double out = 1.0;
  if (style == 0) {
    //standard kelbg
  } else if (style == 1) {
    //electron-electron, with laguerre quadrature
    out = -MY_PI*MY_SQRT2*xi/get_lsee(xi);
  } else if (style == 2) {
    //not implemented yet. Anyone know a numerical way around
    //riemann zeta normalization??
  }
  return out;
}

void PairQspKelbg::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum, itype,jtype;
  double qtmp,xtmp,ytmp,ztmp, delx,dely,delz, fpair;
  double rsq, r2inv, rinv, itemp, teff, s; //extra double s just in case.
  double xi, lambdasq, nusq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *temp = atom->temp;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;
  double hbar = force->hplanck/MY_2PI;
  double kb = force->boltz;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    itemp = temp[i];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      if(on[itype][jtype]) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_kelbgsq[itype][jtype]) {
          teff = 0.5*(temp[j] + itemp);
          r2inv = 1.0/rsq;
          rinv = sqrt(r2inv);
          lambdasq = 1.0/mass[itype];
          lambdasq += 1.0/mass[jtype];
          lambdasq = lambdasq/(kb*teff);
          lambdasq = lambdasq*MY_2PI*hbar*hbar;
          xi = q[i]*q[j]*qqrd2e/(kb*teff);
          xi = xi/sqrt(lambdasq);
          nusq = get_nu(style[itype][jtype], xi);
          nusq = nusq*nusq;
          s = rsq/lambdasq;
          fpair = -2*MY_2PI - s;
          s = -1*s*MY_2PI;
          fpair = fpair*exp(s);
          fpair = fpair + 2*MY_2PI*exp(s*nusq);
          fpair = q[i]*q[j]*qqrd2e*fpair*rinv/lambdasq;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;

          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }

        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQspKelbg::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;


  memory->create(cut_kelbgsq, n+1,n+1,"pair:cut_kelbgsq");
  memory->create(on,n+1,n+1,"pair:on");
  memory->create(style,n+1,n+1,"pair:style");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQspKelbg::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
  //cutoffs will be set per a pair type.
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQspKelbg::coeff(int narg, char **arg)
{
  if (narg != 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_kelbgsq1 = force->numeric(FLERR,arg[2]);
  int ion = force->inumeric(FLERR, arg[3]);
  int istyle = force->inumeric(FLERR, arg[4]);

  cut_kelbgsq1 = cut_kelbgsq1*cut_kelbgsq1;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_kelbgsq[i][j] = cut_kelbgsq1;
      on[i][j] = ion;
      style[i][j] = istyle;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQspKelbg::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cut_kelbgsq[j][i] = cut_kelbgsq[i][j];
  on[j][i] = on[i][j];
  style[j][i] = style[i][j];
  setflag[j][i] = 1;




  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  return cut_kelbgsq[i][j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairQspKelbg::init_style()
{
  if (!atom->q_flag && !atom->temp_flag)
    error->all(FLERR,"Pair style qsp/kelbg requires atom attribute q and temp.");

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQspKelbg::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut_kelbgsq[i][j],sizeof(double),1,fp);
        fwrite(&on[i][j],sizeof(int),1,fp);
        fwrite(&style[i][j],sizeof(int),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQspKelbg::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&cut_kelbgsq[i][j],sizeof(double),1,fp);
          fread(&on[i][j],sizeof(int),1,fp);
          fread(&style[i][j],sizeof(int),1,fp);
        }
        MPI_Bcast(&cut_kelbgsq[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&on[i][j],1,MPI_INT,0,world);
        MPI_Bcast(&style[i][j],1,MPI_INT,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQspKelbg::write_restart_settings(FILE *fp)
{
  fwrite(&cut_kelbgsq,sizeof(double),1,fp);
  fwrite(&on,sizeof(int),1,fp);
  fwrite(&style,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQspKelbg::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_kelbgsq,sizeof(double),1,fp);
    fread(&on,sizeof(int),1,fp);
    fread(&style,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_kelbgsq,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&on,1,MPI_INT,0,world);
  MPI_Bcast(&style,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairQspKelbg::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %d %d\n", i, cut_kelbgsq[i][i], on[i][i], style[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairQspKelbg::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %d %d\n",i,j,
              cut_kelbgsq[i][j], on[i][j], style[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairQspKelbg::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_kelbgsq") == 0) return (void *) &cut_kelbgsq;
  if (strcmp(str,"on") == 0) return (void *) &on;
  if (strcmp(str, "style") == 0) return (void *) &style;
  return NULL;
}
