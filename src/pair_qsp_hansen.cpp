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
#include "pair_qsp_hansen.h"
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

PairQspHansen::PairQspHansen(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = msmflag = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairQspHansen::~PairQspHansen()
{
  if (!copymode) {
    if (allocated) {
      memory->destroy(setflag);

      memory->destroy(cut_hansensq);
      memory->destroy(on);
    }
  }
}


void PairQspHansen::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum, itype,jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, fpair, imass, jmass;
  double rsq, r2inv, rinv, itemp, teff, lambda, en, s;
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
    imass = mass[itype];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      jtype = type[j];
      jmass = mass[jtype];
      if(on[itype][jtype]) {
        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];

        rsq = delx*delx + dely*dely + delz*delz;
        if (rsq < cut_hansensq[itype][jtype]) {
          teff = 0.5*(temp[j] + itemp);
          r2inv = 1.0/rsq;
          rinv = sqrt(r2inv);
          lambda = 1.0/mass[itype];
          lambda += 1.0/mass[jtype];
          lambda = lambda/(kb*teff);
          lambda = lambda*MY_2PI*hbar*hbar;
          lambda = sqrt(lambda);
          s = MY_2PI*sqrt(rsq)/lambda;
          fpair = 1 + s;
          fpair = exp(-1*s)*fpair;
          fpair = -1*qtmp*q[j]*qqrd2e*r2inv*fpair;
          fpair = rinv*fpair;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;

          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }

          if (eflag) {
            en = exp(-1*s);
            en = -1*qtmp*q[j]*qqrd2e*rinv*en;
          }

          if (evflag)
            ev_tally(i, j, nlocal, newton_pair, 0.0, en, fpair,
                     delx, dely, delz);


        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQspHansen::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");
  memory->create(cut_hansensq, n+1,n+1,"pair:cut_kelbgsq");
  memory->create(on,n+1,n+1,"pair:on");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQspHansen::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
  //cutoffs will be set per a pair type.
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQspHansen::coeff(int narg, char **arg)
{
  if (narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_hansensq1 = force->numeric(FLERR,arg[2]);
  int ion = force->inumeric(FLERR, arg[3]);

  cut_hansensq1 = cut_hansensq1*cut_hansensq1;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_hansensq[i][j] = cut_hansensq1;
      on[i][j] = ion;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairQspHansen::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cut_hansensq[j][i] = cut_hansensq[i][j];
  on[j][i] = on[i][j];
  setflag[j][i] = 1;




  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  return cut_hansensq[i][j];
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairQspHansen::init_style()
{
  if (!atom->q_flag && !atom->temp_flag)
    error->all(FLERR,"Pair style qsp/kelbg requires atom attribute q and temp.");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQspHansen::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut_hansensq[i][j],sizeof(double),1,fp);
        fwrite(&on[i][j],sizeof(int),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQspHansen::read_restart(FILE *fp)
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
          fread(&cut_hansensq[i][j],sizeof(double),1,fp);
          fread(&on[i][j],sizeof(int),1,fp);
        }
        MPI_Bcast(&cut_hansensq[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&on[i][j],1,MPI_INT,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQspHansen::write_restart_settings(FILE *fp)
{
  fwrite(&cut_hansensq,sizeof(double),1,fp);
  fwrite(&on,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQspHansen::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_hansensq,sizeof(double),1,fp);
    fread(&on,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_hansensq,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&on,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairQspHansen::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %d\n", i, cut_hansensq[i][i], on[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairQspHansen::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %d\n",i,j,
              cut_hansensq[i][j], on[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairQspHansen::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_hansensq") == 0) return (void *) &cut_hansensq;
  if (strcmp(str,"on") == 0) return (void *) &on;
  return NULL;
}
