#include <cmath>
#include "fdint.h"

using namespace LAMMPS_NS;
using namespace FDINT;

static const double arcih_a[] = {-1.2582793945794000, //a2.5
                                 0.1207822376352453,
                                 0.0233056178489510,
                                 1.0911595094936000,
                                 -0.2993063964300200,
                                 -0.00286186596155192,
                                 0.5051953653801600,
                                 0.0419579806591870,
                                 1.3695261714367000,
                                 0.2685157355131100};

static const double arcih_b[] = {0.0813113962506270,
                                 1.1903358203098999,
                                 1.1445576113258000,
                                 0.2049158578610270};

static const double ih_a[] = {7.940307136e2,
                              1.072518408e3,
                              4.106017002e2,
                              4.607473842e1,
                              1.0};

static const double ih_b[] = {8.959677183e2,
                              1.526979592e3,
                              8.307577602e2,
                              1.638158630e2,
                              9.960923624,
                              1.047712331e-1};

static const double ih_c[] = {7.265461948e-8,
                              1.678032858e-5,
                              1.365376899e-3,
                              4.647886226e-2,
                              5.231390123e-1,
                              1.567714263,
                              1.0};

static const double ih_d[] = {1.089819298e-7,
                              2.503603684e-5,
                              2.017068914e-3,
                              6.719888328e-2,
                              7.001197731e-1,
                              1.309399040,
                              1.727377764e-1};



double arcih(double y) {
  double num = log(y) + arcih_a[0]*pow(y, 5/2.0) + arcih_a[1];
  double den = 1;
  double temp_y = y;
  int i = 1;
  while(i < 5) {
    num += arcih_a[i+1]*temp_y;
    den += arcih_b[i-1]*temp_y;
    temp_y = temp_y*y;
    i++;
  }
  while(i < 8) {
    num += arcih_a[i+1]*temp_y;
    temp_y = temp_y*y;
  }
  num += arcih_a[9]*temp_y*y*y;
  temp_y = num/den;
  return temp_y;
}

double ih(double y) {
  double temp_y;
  if(y < 2)
    {
      y = exp(y);
      temp_y = y;
      double num = ih_a[0];
      double den = ih_b[0];
      for(int i = 1; i < 4; i++) {
        num += ih_a[i]*temp_y;
        den += ih_b[i]*temp_y;
        temp_y = temp_y*y;
      }
      den += ih_b[5]*temp_y;
      temp_y = num/den;
    }
  else {

  }
  return temp_y;
}
