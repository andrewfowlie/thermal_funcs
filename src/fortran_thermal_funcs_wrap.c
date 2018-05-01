/**
    @file file
    @example fortran_example.f
    @brief Fortran wrapper for `thermal_funcs` library.
*/

#include <thermal_funcs.h>

double j_b(double *y_squared) {
  /**
      @returns Bosonic thermal function
      @param y_squared Argument of thermal function
  */
  return J_B_bessel(*y_squared, 1E-7, 1E-7, 10000, false);
}
