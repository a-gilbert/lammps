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
#include "pair_qsp_hmsym.h"
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

PairQspHMSym::PairQspHMSym(LAMMPS *lmp) : Pair(lmp)
{
  //ewaldflag = pppmflag = msmflag = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairQspHMSym::~PairQspHMSym()
{
  if (!copymode) {
    if (allocated) {
      memory->destroy(setflag);

      memory->destroy(cut_hmsq);
      memory->destroy(temp);
      memory->destroy(on);
    }
  }
}


void PairQspHMSym::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum, itype,jtype;
  double xtmp,ytmp,ztmp, delx,dely,delz,fpair;
  double rsq, lambdasq, teff, s, en;
  int *ilist, *jlist, *numneigh, **firstneigh;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
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
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
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
        if (rsq < cut_hmsq[itype][jtype]) {
          teff = temp[i][j];
          lambdasq = 1.0/mass[itype];
          lambdasq += 1.0/mass[jtype];
          lambdasq = lambdasq/(MY_2PI*kb*teff);
          lambdasq = lambdasq*hbar*hbar;
          s = -1*rsq/(MY_PI*lambdasq*log(2));
          fpair = 2*kb*teff*exp(s);
          fpair = fpair/(MY_PI*lambdasq);

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;

          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }

          if (eflag) {
            en = kb*teff*log(2)*exp(s);
          }

          if (evflag)
            ev_tally(i, j, nlocal, newton_pair, 0.0, en,
                     fpair, delx, dely, delz);

        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairQspHMSym::allocate() {
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq, n+1, n+1, "pair:cutsq");
  memory->create(cut_hmsq, n+1,n+1,"pair:cut_kelbgsq");
  memory->create(temp,n+1,n+1,"pair:temp");
  memory->create(on,n+1,n+1,"pair:on");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairQspHMSym::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
  //cutoffs will be set per a pair type.
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairQspHMSym::coeff(int narg, char **arg)
{
  if (narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double cut_hmsq1 = force->numeric(FLERR,arg[2]);
  int ion = force->inumeric(FLERR, arg[3]);

  cut_hmsq1 = cut_hmsq1*cut_hmsq1;

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_hmsq[i][j] = cut_hmsq1;
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

double PairQspHMSym::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  cut_hmsq[j][i] = cut_hmsq[i][j];
  on[j][i] = on[i][j];
  setflag[j][i] = 1;




  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  return sqrt(cut_hmsq[i][j]);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairQspHMSym::init_style()
{
  if (!atom->q_flag && !atom->temp_flag)
    error->all(FLERR,"Pair style qsp/kelbg requires atom attribute q and temp.");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQspHMSym::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut_hmsq[i][j],sizeof(double),1,fp);
        fwrite(&on[i][j],sizeof(int),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQspHMSym::read_restart(FILE *fp)
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
          fread(&cut_hmsq[i][j],sizeof(double),1,fp);
          fread(&on[i][j],sizeof(int),1,fp);
        }
        MPI_Bcast(&cut_hmsq[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&on[i][j],1,MPI_INT,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairQspHMSym::write_restart_settings(FILE *fp)
{
  fwrite(&cut_hmsq,sizeof(double),1,fp);
  fwrite(&on,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairQspHMSym::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_hmsq,sizeof(double),1,fp);
    fread(&on,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_hmsq,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&on,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairQspHMSym::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %d\n", i, cut_hmsq[i][i], on[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairQspHMSym::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %d\n",i,j,
              cut_hmsq[i][j], on[i][j]);
}

/* ---------------------------------------------------------------------- */

void *PairQspHMSym::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"cut_fsq") == 0) return (void *) &cut_hmsq;
  if (strcmp(str,"on") == 0) return (void *) &on;
  return NULL;
}
