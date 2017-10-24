/*

  Derivatives of thermal functions.

*/


#include <derivatives.h>
#include <zeta.h>
#include <thermal_funcs.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <stdexcept>
#include <algorithm>


// Derivatives of thermal functions at y_squared -> 0. Found in Mathematica:
//
// -Sum[1/n^2*Limit[D[xsq*BesselK[2,Sqrt[xsq]*n],xsq],xsq->0],{n,1,Infinity}]
// -Sum[(-1)^n/n^2*Limit[D[xsq*BesselK[2, Sqrt[xsq]*n], xsq], xsq -> 0], {n, 1, Infinity}]
//
// or from Taylor expansions.

const double D1_J_B_0 = pow(M_PI, 2) / 12.;
const double D1_J_F_0 = -pow(M_PI, 2) / 24.;
const double INF = 1E10;

// Bessel function representation

double D_K2(double x, int order) {
  // Derivative of K2 Bessel function.
  // D[BesselK[2, x], {x, n}]

  #ifdef THROW
  if (order < 0) {
    throw std::runtime_error("Derivative order must be >= 0");
  }
  #endif

  gsl_set_error_handler_off();  // Default handler aborts

  double d = 0.;

  for (int i = 0; i <= order; i += 1) {
    d += gsl_sf_bessel_Kn(2 + 2 * i - order, x) * gsl_sf_choose(order, i);
  }

  d *= pow(-0.5, order);

  return d;
}

double D_Y2(double x, int order) {
  // Derivative of Y2 Bessel function.
  // D[BesselY[2, x], {x, n}]

  #ifdef THROW
  if (order < 0) {
    throw std::runtime_error("Derivative order must be >= 0");
  }
  #endif

  gsl_set_error_handler_off();  // Default handler aborts

  double d = 0.;

  for (int i = 0; i <= order; i += 1) {
    d += pow(-1., i) * gsl_sf_bessel_Yn(2 + 2 * i - order, x) * gsl_sf_choose(order, i);
  }

  d *= pow(-0.5, order);

  return d;
}

double D_real_K2(cdouble x, int order) {
  // Derivative of real part of K2 Bessel function.
  // Utilize fact that
  // Re[BesselK[2, x * I]] = 0.5 * Pi BesselY[2, x]
  // to define K2 for imaginary arguments.

  #ifdef THROW
  if (order < 0) {
    throw std::runtime_error("Derivative order must be >= 0");
  }
  #endif

  if (real(x) != 0. && imag(x) != 0.) {
    #ifdef THROW
      throw std::invalid_argument("K2 only implemented for x real or imaginary");
    #endif
  } else if (real(x) != 0.) {
    return D_K2(real(x), order);
  } else if (imag(x) != 0.) {
    return 0.5 * M_PI *D_Y2(imag(x), order);
  }
}

double D1_bessel_sum(double y_squared, double abs_error, double rel_error, int max_n, bool bosonic) {
  // Bessel sum of derivative thermal function with respect to y^2.

  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  cdouble y = sqrt(cdouble(y_squared));
  double sign = 2. * static_cast<double>(bosonic) - 1.;
  double sign_y_squared = (y_squared >= 0.) ? 1. : -1.;
  double factor = sign;
  double sum = - factor * (D_real_K2(y, 0) + 0.5 * sign_y_squared * std::abs(y) * D_real_K2(y, 1));

  for (int n = 2; n <= max_n; n += 1) {
    factor *= sign * pow((static_cast<double>(n) - 1.) / static_cast<double>(n), 2);
    const cdouble arg = static_cast<double>(n) * y;
    const double term = - factor * (D_real_K2(arg, 0) + 0.5 * n * sign_y_squared * std::abs(y) * D_real_K2(arg, 1));
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

double D2_bessel_sum(double y_squared, double abs_error, double rel_error, int max_n, bool bosonic) {
  // Bessel sum of derivative thermal function with respect to y^2.

  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  cdouble y = sqrt(cdouble(y_squared));
  double sign_y_squared = (y_squared >= 0.) ? 1. : -1.;
  double sign = 2. * static_cast<double>(bosonic) - 1.;
  double factor = sign;
  double sum = - factor * (D_real_K2(y, 1) * 0.75 / std::abs(y)  + D_real_K2(y, 2) * sign_y_squared * 0.25);

  for (int n = 2; n <= max_n; n += 1) {
    factor *= sign * pow((static_cast<double>(n) - 1.) / static_cast<double>(n), 2);
    const cdouble arg = static_cast<double>(n) * y;
    const double term = - factor * (D_real_K2(arg, 1) / std::abs(y) * (0.5 * n + 0.25 * pow(n, 2))
                                     + D_real_K2(arg, 2) * 0.25 * sign_y_squared * pow(n, 2));
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

double D1_J_F_bessel(double y_squared, double abs_error, double rel_error, int max_n) {
  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return D1_J_F_0;
  }
  return D1_bessel_sum(y_squared, abs_error, rel_error, max_n, false);
}

double D1_J_B_bessel(double y_squared, double abs_error, double rel_error, int max_n) {
  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return D1_J_B_0;
  }
  return D1_bessel_sum(y_squared, abs_error, rel_error, max_n, true);
}

double D2_J_F_bessel(double y_squared, double abs_error, double rel_error, int max_n) {
  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    #ifdef THROW
      throw std::runtime_error("second derivative diverges at 0.");
      return INF;
    #endif
  }
  return D2_bessel_sum(y_squared, abs_error, rel_error, max_n, false);
}

double D2_J_B_bessel(double y_squared, double abs_error, double rel_error, int max_n) {
  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    #ifdef THROW
      throw std::runtime_error("second derivative diverges at 0.");
      return INF;
    #endif
  }
  return D2_bessel_sum(y_squared, abs_error, rel_error, max_n, true);
}


// Numerical derivative


struct J_params {bool bosonic; double abs_error; double rel_error; int max_n;};

double J_wrapper(double y_squared, void *p) {
  // Wrapper for J with signature required by numerical differentiation
  struct J_params *params = (struct J_params *)p;
  bool bosonic = params->bosonic;
  double abs_error = params->abs_error; 
  double rel_error = params->rel_error; 
  int max_n = params->max_n;
  
  if (bosonic) {
    return J_B_bessel(y_squared, abs_error, rel_error, max_n);
  } else {
    return J_F_bessel(y_squared, abs_error, rel_error, max_n);
  }
}

double D1_J_approx(double y_squared, double step, bool bosonic, double abs_error, double rel_error, int max_n) {
    gsl_function J;
    double derivative;
    double abs_err;

    J.function = &J_wrapper;
    struct J_params params = {bosonic, abs_error, rel_error, max_n};
    J.params = &params;

    gsl_deriv_central(&J, y_squared, step, &derivative, &abs_err);

    return derivative;
}

struct D_params {double step; bool bosonic; double abs_error; double rel_error; int max_n;};

double D1_J_wrapper(double y_squared, void *p) {
  // Wrapper for first derivative with signature required by numerical differentiation
  struct D_params *params = (struct D_params *)p;
  double step = params->step;
  bool bosonic = params->bosonic;
  double abs_error = params->abs_error; 
  double rel_error = params->rel_error; 
  int max_n = params->max_n;
  return D1_J_approx(y_squared, step, bosonic, abs_error, rel_error, max_n);
}

double D2_J_approx(double y_squared, double step, bool bosonic, double abs_error, double rel_error, int max_n) {
    gsl_function D1_J;
    double derivative;
    double abs_err;

    D1_J.function = &D1_J_wrapper;
    struct D_params params = {step, bosonic, abs_error, rel_error, max_n};
    D1_J.params = &params;

    gsl_deriv_central(&D1_J, y_squared, step, &derivative, &abs_err);

    return derivative;
}

double D1_J_B_approx(double y_squared, double step, double abs_error, double rel_error, int max_n) {
  return D1_J_approx(y_squared, step, true, abs_error, rel_error, max_n);
}

double D1_J_F_approx(double y_squared, double step, double abs_error, double rel_error, int max_n) {
  return D1_J_approx(y_squared, step, false, abs_error, rel_error, max_n);
}

double D2_J_B_approx(double y_squared, double step, double abs_error, double rel_error, int max_n) {
  return D2_J_approx(y_squared, step, true, abs_error, rel_error, max_n);
}

double D2_J_F_approx(double y_squared, double step, double abs_error, double rel_error, int max_n) {
  return D2_J_approx(y_squared, step, false, abs_error, rel_error, max_n);
}




