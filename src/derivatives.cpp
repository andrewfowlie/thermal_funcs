/*

  Derivatives of Bessel function representation of thermal functions.

*/


#include <derivatives.h>
#include <zeta.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <complex>


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
  // Derivative of sum of Bessel functions with respect to y^2.

  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  cdouble y = sqrt(cdouble(y_squared));
  double sign = 2. * static_cast<double>(bosonic) - 1.;
  double factor = sign;
  double sum = factor * D_real_K2(y, 1);

  for (int n = 2; n <= max_n; n += 1) {
    factor *= sign * (static_cast<double>(n) - 1.) / static_cast<double>(n);
    const double term = factor * D_real_K2(static_cast<double>(n) * y, 1);
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
  
  sum *= 0.5 / std::abs(y);

  return sum;
}

double D2_bessel_sum(double y_squared, double abs_error, double rel_error, int max_n, bool bosonic) {
  // Second derivative of sum of Bessel functions with respect to y^2.

  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  cdouble y = sqrt(cdouble(y_squared));
  double sign = 2. * static_cast<double>(bosonic) - 1.;
  double factor = sign;
  double sum = factor * D_real_K2(y, 1);

  for (int n = 2; n <= max_n; n += 1) {
    factor *= sign;
    const double term = factor * D_real_K2(static_cast<double>(n) * y, 1);
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
  
  sum *= 0.25 / std::abs(y_squared);

  return sum;
}

double D1_J_bessel(double y_squared, double abs_error, double rel_error, int max_n, bool fast, bool bosonic) {
  // Derivative of thermal function with respect to y^2

  const double d = - y_squared * D1_bessel_sum(y_squared, abs_error, rel_error, max_n, bosonic) 
                    + bessel_sum(y_squared, abs_error, rel_error, max_n, fast, bosonic) / y_squared;
  
  return d;
}

double D2_J_bessel(double y_squared, double abs_error, double rel_error, int max_n, bool fast, bool bosonic) {
  // Second derivative of thermal function with respect to y^2  

  const double d = - y_squared * D2_bessel_sum(y_squared, abs_error, rel_error, max_n, bosonic) 
                    - 0.5 * D1_bessel_sum(y_squared, abs_error, rel_error, max_n, fast, bosonic);
  
  return d;
}






