#include <stdbool.h>
#include <wstp.h>
#include <thermal_funcs.h>
#include <derivatives.h>


double J_F(double y_squared, int order) {
  /**
      @returns Fermionic thermal function from Bessel functions
      @param y_squared Argument of thermal function
      @param order Order of derivative
  */
  if (order == 0) {
    return J_F_bessel(y_squared, 1E-7, 1E-7, 10000, false);
  } else if (order == 1) {
    return D1_J_F_bessel(y_squared, 1E-7, 1E-7, 10000);
  } else if (order == 2) {
    return D2_J_F_bessel(y_squared, 1E-7, 1E-7, 10000);
  } else {
    WSEvaluateString(stdlink, "JF::derivative=\"derivative must be 0, 1 or 2.\"; Message[JF::derivative]");
    return -1.;
  }
}

double J_B(double y_squared, int order) {
  /**
      @returns Bosonic thermal function from Bessel functions
      @param y_squared Argument of thermal function
      @param order Order of derivative
  */
  if (order == 0) {
    return J_B_bessel(y_squared, 1E-7, 1E-7, 10000, false);
  } else if (order == 1) {
    return D1_J_B_bessel(y_squared, 1E-7, 1E-7, 10000);
  } else if (order == 2) {
    return D2_J_B_bessel(y_squared, 1E-7, 1E-7, 10000);
  } else {
    WSEvaluateString(stdlink, "JB::derivative=\"derivative must be 0, 1 or 2.\"; Message[JB::derivative]");
    return -1.;
  }
}

int main(int argc, char *argv[]) {
  return WSMain(argc, argv);
}
