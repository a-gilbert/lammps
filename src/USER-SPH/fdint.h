/* Basic pade approximants for fermi dirac integrals and inverse fermi dirac
   integrals from Karasiev et Al. 2014, and Antia 1993.
 */


#ifndef LMP_FDINT_H
#define LMP_FDINT_H

namespace LAMMPS_NS {
  namespace FDINT {
    //Approximate value for ArcI_{1/2}(y)
    double arcih(double y);
    //Approximate value for I_{1/2}(y)
    double ih(double y);
    //Approximate value for I_{-1/2}(y)
    double inh(double y);
    //Approximate value for I_{3/2}(y)
    double i3h(double y);

  }
}

#endif
