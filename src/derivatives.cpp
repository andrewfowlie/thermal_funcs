/*

  First and second derivatives of thermal functions with respect to y^2.

  The Bessel representation is differentiated analytically and numerically.

*/


#include <derivatives.h>
#include <zeta.h>
#include <thermal_funcs.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>
#include <algorithm>


// First derivatives of thermal functions with respect to y^2 at y^2 -> 0.
// Found in Mathematica:

// -Sum[1/n^2*Limit[D[xsq*BesselK[2,Sqrt[xsq]*n],xsq],xsq->0],{n,1,Infinity}]
// -Sum[(-1)^n/n^2*Limit[D[xsq*BesselK[2, Sqrt[xsq]*n], xsq], xsq -> 0], {n, 1, Infinity}]


const double D1_J_B_0 = pow(M_PI, 2) / 12.;
const double D1_J_F_0 = -pow(M_PI, 2) / 24.;
const double INF = 1E10;


// Bessel function representation of derivatives with respect to y^2.


double K1(cdouble x) {
  /**
      @returns K1 Bessel function.
  */
  gsl_set_error_handler_off();  // Default handler aborts

  if (real(x) != 0. && imag(x) != 0.) {
    #ifdef THROW
      throw std::invalid_argument("K1 only implemented for x real or imaginary");
    #endif
  } else  if (std::abs(x) == 0.) {
    #ifdef THROW
      throw std::invalid_argument("K1 diverges for |x| = 0");
    #endif
  } else if (real(x) != 0.) {
    return gsl_sf_bessel_Kn(1, real(x));
  } else if (imag(x) != 0.) {
      return -0.5 * M_PI * gsl_sf_bessel_Yn(1, imag(x));
  }
}

double K0(cdouble x) {
  /**
      @returns K0 Bessel function.

  */
  gsl_set_error_handler_off();  // Default handler aborts

  if (real(x) != 0. && imag(x) != 0.) {
    #ifdef THROW
      throw std::invalid_argument("K0 only implemented for x real or imaginary");
    #endif
  } else  if (std::abs(x) == 0.) {
    #ifdef THROW
      throw std::invalid_argument("K0 diverges for |x| = 0");
    #endif
  } else if (real(x) != 0.) {
    return gsl_sf_bessel_Kn(0, real(x));
  } else if (imag(x) != 0.) {
      return -0.5 * M_PI * gsl_sf_bessel_Yn(0, imag(x));
  }
}

double D1_bessel_sum(double y_squared, double abs_error, double rel_error,
                     int max_n, bool bosonic) {
  /**
      @returns First derivative thermal function with respect to y^2.
      Found by summing Bessel functions.

      @param y_squared Argment of thermal function
  */
  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  const double sign = 2. * static_cast<double>(bosonic) - 1.;
  const cdouble y = sqrt(cdouble(y_squared));
  const double abs_y = std::abs(y);
  double factor = 0.5 * abs_y * sign;
  double sum = factor * K1(y);

  for (int n = 2; n <= max_n; n += 1) {
    const double n_double = static_cast<double>(n);
    factor *= sign * (n_double - 1.) / n_double;
    const double term = factor * K1(n_double * y);
    sum += term;

    #ifdef DEBUG
      printf("term_%d = %e and sum = %e\n", n, term, sum);
    #endif

    if (std::abs(term) < std::max(abs_error, rel_error * std::abs(sum))) {
      #ifdef DEBUG
        printf("number of terms in sum = %d\n", n);
      #endif
      break;
    } else if (std::isinf(sum)) {
      #ifdef THROW
        throw std::runtime_error("sum diverging");
      #endif
      break;
    }

    #ifdef DEBUG
    if (n == max_n) {
      printf("reached max number of terms in sum = %d\n", n);
    }
    #endif
  }

  return sum;
}

double D2_bessel_sum(double y_squared, double abs_error, double rel_error,
                     int max_n, bool bosonic) {
  //
  /**
      @returns Second derivative thermal function with respect to y^2.
      Found by summing Bessel functions.

      @param y_squared Argment of thermal function
  */

  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  const double sign = 2. * static_cast<double>(bosonic) - 1.;
  const cdouble y = sqrt(cdouble(y_squared));
  double factor = -0.25 * sign;
  double sum = factor * K0(y);

  for (int n = 2; n <= max_n; n += 1) {
    const double n_double = static_cast<double>(n);
    factor *= sign;
    const double term = factor * K0(n_double * y);
    sum += term;

    #ifdef DEBUG
      printf("term_%d = %e and sum = %e\n", n, term, sum);
    #endif

    if (std::abs(term) < std::max(abs_error, rel_error * std::abs(sum))) {
      #ifdef DEBUG
        printf("number of terms in sum = %d\n", n);
      #endif
      break;
    } else if (std::isinf(sum)) {
      #ifdef THROW
        throw std::runtime_error("sum diverging");
      #endif
      break;
    }

    #ifdef DEBUG
    if (n == max_n) {
      printf("reached max number of terms in sum = %d\n", n);
    }
    #endif
  }

  return sum;
}

double D1_J_F_bessel(double y_squared, double abs_error, double rel_error,
                     int max_n) {
  /**
      @returns First derivative fermionic thermal function with respect to y^2.
      Found by summing Bessel functions.

      @param y_squared Argment of thermal function
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return D1_J_F_0;
  }
  return D1_bessel_sum(y_squared, abs_error, rel_error, max_n, false);
}

double D1_J_B_bessel(double y_squared, double abs_error, double rel_error,
                     int max_n) {
  /**
      @returns First derivative bosonic thermal function with respect to y^2.
      Found by summing Bessel functions.

      @param y_squared Argment of thermal function
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return D1_J_B_0;
  }
  return D1_bessel_sum(y_squared, abs_error, rel_error, max_n, true);
}

double D2_J_F_bessel(double y_squared, double abs_error, double rel_error,
                     int max_n) {
  /**
      @returns Second derivative fermionic thermal function with respect to y^2.
      Found by summing Bessel functions.

      @param y_squared Argment of thermal function
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    #ifdef THROW
      throw std::runtime_error("second derivative diverges at 0.");
    #endif
    return INF;
  }
  return D2_bessel_sum(y_squared, abs_error, rel_error, max_n, false);
}

double D2_J_B_bessel(double y_squared, double abs_error, double rel_error,
                     int max_n) {
  /**
      @returns Second derivative bosonic thermal function with respect to y^2.
      Found by summing Bessel functions.

      @param y_squared Argment of thermal function
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    #ifdef THROW
      throw std::runtime_error("second derivative diverges at 0.");
    #endif
    return INF;
  }
  return D2_bessel_sum(y_squared, abs_error, rel_error, max_n, true);
}


// Numerical derivatives of thermal functions with respect to y^2.


struct J_params {
  bool bosonic;
  double abs_error;
  double rel_error;
  int max_n;
};

double J_wrapper(double y_squared, void *p) {
  /**
     Wrapper for thermal functions with signature required by numerical
     differentiation.
  */
  const struct J_params *params = (struct J_params *)p;
  const bool bosonic = params->bosonic;
  const double abs_error = params->abs_error;
  const double rel_error = params->rel_error;
  const int max_n = params->max_n;

  if (bosonic) {
    return J_B_bessel(y_squared, abs_error, rel_error, max_n);
  } else {
    return J_F_bessel(y_squared, abs_error, rel_error, max_n);
  }
}

double D1_J_approx(double y_squared, double step, bool bosonic,
                   double abs_error, double rel_error, int max_n) {
  /**
     @returns Numerical derivative of thermal function.
  */
  gsl_function J;
  double derivative;
  double abs_err;

  J.function = &J_wrapper;
  struct J_params params = {bosonic, abs_error, rel_error, max_n};
  J.params = &params;

  gsl_deriv_central(&J, y_squared, step, &derivative, &abs_err);

  return derivative;
}

struct D_params {
  double step;
  bool bosonic;
  double abs_error;
  double rel_error;
  int max_n;
};

double D1_J_wrapper(double y_squared, void *p) {
  /**
     Wrapper for first derivative of thermal functions with signature required
     by numerical differentiation.
  */
  const struct D_params *params = (struct D_params *)p;
  const double step = params->step;
  const bool bosonic = params->bosonic;
  const double abs_error = params->abs_error;
  const double rel_error = params->rel_error;
  const int max_n = params->max_n;
  return D1_J_approx(y_squared, step, bosonic, abs_error, rel_error, max_n);
}

double D2_J_approx(double y_squared, double step, bool bosonic,
                   double abs_error, double rel_error, int max_n) {
  /**
     @returns Numerical second derivative of thermal function.
  */
  gsl_function D1_J;
  double derivative;
  double abs_err;

  D1_J.function = &D1_J_wrapper;
  struct D_params params = {step, bosonic, abs_error, rel_error, max_n};
  D1_J.params = &params;

  gsl_deriv_central(&D1_J, y_squared, step, &derivative, &abs_err);

  return derivative;
}

double D1_J_B_approx(double y_squared, double step, double abs_error,
                     double rel_error, int max_n) {
  /**
     @returns Numerical derivative of bosonic thermal function.
  */
  return D1_J_approx(y_squared, step, true, abs_error, rel_error, max_n);
}

double D1_J_F_approx(double y_squared, double step, double abs_error,
                     double rel_error, int max_n) {
  /**
     @returns Numerical derivative of fermionic thermal function.
  */
  return D1_J_approx(y_squared, step, false, abs_error, rel_error, max_n);
}

double D2_J_B_approx(double y_squared, double step, double abs_error,
                     double rel_error, int max_n) {
  /**
     @returns Numerical second derivative of bosonic thermal function.
  */
  return D2_J_approx(y_squared, step, true, abs_error, rel_error, max_n);
}

double D2_J_F_approx(double y_squared, double step, double abs_error,
                     double rel_error, int max_n) {
  /**
     @returns Numerical second derivative of fermionic thermal function.
  */
  return D2_J_approx(y_squared, step, false, abs_error, rel_error, max_n);
}
