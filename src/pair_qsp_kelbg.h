/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/*
  A short range correction to coulomb interactions that accounts for the distinguishable
  nature of particles in quantum statistics at a given temperature. Designed only for
  ion-electron and electron-electron interactions.

  See 
  Jones, Christopher S., and Michael S. Murillo. “Analysis of Semi-Classical 
  Potentials for Molecular Dynamics and Monte Carlo Simulations of Warm 
  Dense Matter.” High Energy Density Physics 3, no. 3–4 (October 2007): 379–94.
   https://doi.org/10.1016/j.hedp.2007.02.038.

 */
#ifdef PAIR_CLASS

PairStyle(qsp/kelbg,PairQspKelbg)

#else

#ifndef LMP_PAIR_QSP_KELBG_LONG_H
#define LMP_PAIR_QSP_KELBG_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairQspKelbg : public Pair {
 public:
  PairQspKelbg(class LAMMPS *);
  virtual ~PairQspKelbg();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *);
  void write_data_all(FILE *);
  virtual void *extract(const char *, int &);

 protected:
  double **cut_kelbgsq;
  int **on, **style;

  double *cut_respa;

  virtual void allocate();

 private:
  double j1(double, double);
  double get_j1(double);
  double get_lsei(double);
  double get_lsee(double);
  double get_nu(int, double);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Pair style buck/coul/long requires atom attribute q

The atom style defined does not have these attributes.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
