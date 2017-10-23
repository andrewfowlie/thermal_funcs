#include <stdbool.h>
#include <wstp.h>
#include <thermal_funcs.h>
#include <derivatives.h>


double J_F(double y_squared, int n) {
   if (n == 0) {
      return J_F_bessel(y_squared, 1E-7, 1E-7, 10000, false);
   } else if (n == 1) {
      return D1_J_F_bessel(y_squared, 1E-7, 1E-7, 10000);
   } else if (n == 2) {
      return D2_J_F_bessel(y_squared, 1E-7, 1E-7, 10000);
   } else {
      return -1;
   }
}

double J_B(double y_squared, int n) {
   if (n == 0) {
      return J_B_bessel(y_squared, 1E-7, 1E-7, 10000, false);
   } else if (n == 1) {
      return D1_J_B_bessel(y_squared, 1E-7, 1E-7, 10000);
   } else if (n == 2) {
      return D2_J_B_bessel(y_squared, 1E-7, 1E-7, 10000);
   } else {
      return -1;
   }
}

int main(int argc, char *argv[]) {
   return WSMain(argc, argv);
}
