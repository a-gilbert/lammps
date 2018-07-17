#include <cmath>
#include "fdint.h"

namespace LAMMPS_NS {
  namespace FDINT {

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

    static const double i3h_a[] = {9.895512903e2,
                                   1.237156375e3,
                                   4.413986183e2,
                                   4.693212727e1,
                                   1.0};

    static const double i3h_b[] = {7.443927085e2,
                                   1.062245497e3,
                                   4.720721124e2,
                                   7.386867306e1,
                                   3.424526047,
                                   2.473929073e-2};

    static const double i3h_c[] = {6.7384341042e-8,
                                   7.4281282702e-6,
                                   4.6220789293e-4,
                                   1.1905625478e-2,
                                   1.3661062300e-1,
                                   6.5500705397e-1,
                                   1.0};

    static const double i3h_d[] = {1.6846085253e-7,
                                   1.7531170088e-5,
                                   1.0476768850e-3,
                                   2.3334235654e-2,
                                   1.9947560547e-1,
                                   4.7103657850e-1,
                                   -1.7443752246e-2};



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
        i++;
      }
      num += arcih_a[9]*temp_y*y*y;
      temp_y = num/den;
      return temp_y;
    }

    double ih(double y) {
      double temp_y;
      double num;
      double den;
      if(y < 2)
        {
          y = exp(y);
          temp_y = y;
          num = ih_a[0];
          den = ih_b[0];
          for(int i = 1; i < 5; i++) {
            num += ih_a[i]*temp_y;
            den += ih_b[i]*temp_y;
            temp_y = temp_y*y;
          }
          den += ih_b[5]*temp_y;
          temp_y = y*num/den;
        }
      else {
        temp_y = 1/(y*y);
        num = ih_c[0];
        den = ih_d[0];
        for(int i = 1; i < 7; i++) {
          num += ih_c[i]*temp_y;
          den += ih_d[i]*temp_y;
          temp_y = temp_y*temp_y;
        }
        temp_y = pow(y, 1.5)*num/den;
      }
      return temp_y;
    }


    double i3h(double y) {
      double temp_y;
      double num;
      double den;
      if(y < 2)
        {
          y = exp(y);
          temp_y = y;
          num = i3h_a[0];
          den = i3h_b[0];
          for(int i = 1; i < 5; i++) {
            num += ih_a[i]*temp_y;
            den += ih_b[i]*temp_y;
            temp_y = temp_y*y;
          }
          den += ih_b[5]*temp_y;
          temp_y = y*num/den;
        }
      else {
        temp_y = 1/(y*y);
        num = i3h_c[0];
        den = i3h_d[0];
        for(int i = 1; i < 7; i++) {
          num += ih_c[i]*temp_y;
          den += ih_d[i]*temp_y;
          temp_y = temp_y*temp_y;
        }
        temp_y = pow(y, 2.5)*num/den;
      }
      return temp_y;
    }
  }
}
