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

    static const double inh_a[] = {8.830316038e2,
                                   1.183989392e3,
                                   4.473770672e2,
                                   4.892542028e1,
                                   1.00};

    static const double inh_b[] = {4.981972343e2,
                                   1.020272984e3,
                                   6.862151992e2,
                                   1.728621255e2,
                                   1.39857990e1,
                                   2.138408204e-1};

    static const double inh_c[] = {-4.9141019880e-8,
                                   -4.2786358805e-6,
                                   -7.4382915429e-4,
                                   -3.2856045308e-2,
                                   -5.6853219702e-1,
                                   -1.9284139162,
                                   1.000};

    static const double inh_d[] = {-2.4570509894e-8,
                                   -3.6344227710e-6,
                                   -3.7345152736e-4,
                                   -1.6589736860e-2,
                                   -2.9154391835e-1,
                                   -1.1843742874,
                                   7.0985168479e-1,
                                   -6.0197789199e-2};

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


    //For public record: I loathe both this function
    //and the heat capacity of
    //fermions at finite temperature.
    double arcih(double y) {
      double u_init = pow(y, 2.0/3);
      double v_init = u_init*u_init;
      double u = u_init;
      double v = v_init;
      double num = log(y) + arcih_a[0]*pow(u, 5.0/2) + arcih_a[1];
      double den = 1;
      int i = 1;
      while(i < 5) {
        num += arcih_a[i+1]*u;
        den += arcih_b[i-1]*v;
        u = u*u_init; //u = u^{i+1}
        v = v*v_init; //v = v^{i+1}
        i++;
      }
      while(i < 8) {
        num += arcih_a[i+1]*u;
        u = u * u_init;
        i++;
      }
      num += arcih_a[9]*u*u_init*u_init;
      num = num/den;
      return num;
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


    double inh(double y) {
      double temp_y;
      double num;
      double den;
      if (y < 2) {
        y = exp(y);
        temp_y = y;
        num = inh_a[0];
        den = inh_b[0];
        for (int i = 1; i < 5; i++) {
          num += inh_a[i]*temp_y;
          den += inh_b[i]*temp_y;
          temp_y = temp_y*y;
        }
        den += inh_b[5]*temp_y;
        temp_y = y*num/den;

      } else {
        temp_y = 1/(y*y);
        num = inh_c[0];
        den = inh_d[0];
        for (int i = 1; i < 7; i++) {
          num += inh_c[i]*temp_y;
          den += inh_d[i]*temp_y;
          temp_y = temp_y*y;
        }
        den += inh_d[7]*temp_y;
        temp_y = pow(y, 0.5)*num/den;
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
