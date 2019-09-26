/**
    @file
    @brief Thermal functions in finite-temperature field theory.

    I refer throughout to <a href="https://arxiv.org/pdf/1612.00466.pdf">
    Thermal Resummation and Phase Transitions by Curtin et al.</a>
    The thermal functions are eq. (2.12):
    \f[
    J_{B/F}(y^2)=\Re\int_0^{\infty} dx\,x^2 \ln\left(1\mp\exp\left(-\sqrt{x^2 + y^2}\right)\right).
    \f]
    We calculate only the real part.

    I also refer to <a href="http://escholarship.org/uc/item/2r84s7h9">
    Phase Transitions in the Early Universe by Wainwright</a>
    for eq. (2.18) and (2.19).
    
    @example fortran_example.f
    @example example.cpp
*/


#include <thermal_funcs.h>
#include <zeta.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_fermi_dirac.h>
#include <gsl/gsl_sf_pow_int.h>
#include <stdexcept>
#include <complex>
#include <algorithm>
#include <limits>

/*! \f$y^2\f$ that is considered \f$y^2 \ll 0\f$ */
constexpr double neg_y_squared = -1.E3;
/*! \f$y^2\f$ that is considered \f$y^2 \gg 1\f$ */
constexpr double pos_y_squared = 1.E3;
// 
/*! Term in bosonic sum defined below eq. (2.13) */
constexpr double a_b = pow(M_PI, 2) * exp(1.5 - 2. * M_EULER);
/*! Term in fermionic sum defined below eq. (2.13) */
constexpr double a_f = 16. * a_b;
/*! \f$\pi^2\f$ */
constexpr double M_PI_POW_2 = pow(M_PI, 2);
/*! \f$1/\pi^2\f$ */
constexpr double M_PI_POW_M2 = pow(M_PI, -2);
/*! \f$\pi^4\f$ */
constexpr double M_PI_POW_4 = pow(M_PI, 4);

/*!
    Bosonic thermal function at \f$y^2 \to 0\f$. Found in Mathematica:
    `-Sum[1/n^2 * Limit[x^2 * BesselK[2, x *n], x -> 0], {n, 1, Infinity}]`
    or from Taylor expansions.
*/
constexpr double J_B_0 = - pow(M_PI, 4) / 45.;
/*!
    Fermionic thermal function at \f$y^2 \to 0\f$. Found in Mathematica:
    `-Sum[(-1)^n/n^2 * Limit[x^2 * BesselK[2, x *n], x -> 0], {n, 1, Infinity}]`
    or from Taylor expansions.
*/
constexpr double J_F_0 = 7. * pow(M_PI, 4) / 360.;


// Thermal functions by numerical integration


double J_integrand(double x, double y_squared, bool bosonic) {
  /**
      @returns Integrand for thermal function

      Integrand in Curtin eq. (2.12)
      
      \f[x^2 \log\left(1 \mp \exp\left(-\sqrt{x^2 + y^2}\right)\right)\f]
      
      If \f$r^2 = x^2 + y^2 < 0\f$, this can be written
      
      \f[\frac12 x^2 \log\left(2 + 2\cos r\right)\f]

      @param x Argument of integrand
      @param y_squared Argument of thermal function
      @param bosonic Whether bosonic (or fermionic) function required
  */
  const double sign = 1. - 2. * static_cast<double>(bosonic);
  const double x_squared = gsl_sf_pow_int(x, 2);
  const double r_squared = x_squared + y_squared;
  const double abs_r = sqrt(std::abs(r_squared));
  if (r_squared >= 0.) {
    return x_squared * gsl_log1p(sign * exp(-abs_r));
  } else {
    return 0.5 * x_squared * (gsl_log1p(sign * cos(abs_r)) + M_LN2);
  }
}

/** Parameters for quadrature contained in a struct */
struct J_integrand_params {
  /*! Argument of integrand */
  double y_squared;
  /*! Whether bosonic (or fermionic) function required */
  bool bosonic;
};

double J_integrand_wrapper(double x, void *p) {
  /**
      @brief Wrapper for J_integrand with signature required by
      numerical integration.

      @returns Integrand for thermal function

      @param x Argument of integrand
      @param p Additional arguments including whether bosonic function required
  */
  const struct J_integrand_params *params = (struct J_integrand_params *)p;
  const double y_squared = params->y_squared;
  const bool bosonic = params->bosonic;
  return J_integrand(x, y_squared, bosonic);
}

int n_integrand_points(double y_squared, const bool bosonic) {
  /**
      @returns Number of integrand points, i.e. number of singularities + 2 for
      endpoints of integration

      @warning We don't include an endpoint twice if it is singular.
      This means we exclude \f$n = 0\f$.

      @param y_squared Argument of thermal function
      @param bosonic Whether bosonic (or fermionic) function required
  */
  #ifdef THROW
    if (y_squared >= 0.) {
      throw std::invalid_argument("|y_squared| >= 0. - no singularities possible");
    }
  #endif

  const double y = sqrt(std::abs(y_squared));
  const int max_n = floor(y * M_1_PI);
  int n_singularity;

  if (bosonic) {
    // Even numbers <= max_n (excluding 0)
    n_singularity = static_cast<int>(floor(0.5 * max_n));
  } else {
    // Odd numbers <= max_n
    n_singularity = static_cast<int>(ceil(0.5 * max_n));
  }

  return n_singularity + 2;
}

double *integrand_points(double y_squared, bool bosonic) {
  /**
      @returns Singularities in integrand, present if \f$y^2 < 0\f$,
      and boundaries of integration, in ascending order

      @warning We don't include an endpoint twice if it is singular.

      Singularities occur at
      \f[0 \le x = \sqrt{-n^2 \pi^2 - y^2}\le |y|\f]
      for \f$n\f$ even/odd for bosonic/fermionic.

      @param y_squared Argument of thermal function
      @param bosonic Whether bosonic (or fermionic) function required
  */
  #ifdef THROW
    if (y_squared >= 0.) {
      throw std::invalid_argument("|y_squared| >= 0."
                                   " - no singularities possible");
    }
  #endif

  const double y = sqrt(std::abs(y_squared));
  const int n_singularity = n_integrand_points(y_squared, bosonic) - 2;
  double *p = reinterpret_cast<double *>(malloc(sizeof(double)
    * (n_singularity + 2)));

  for (int i = 1; i <= n_singularity; i += 1) {
    // Insure result is in ascending order
    const int reverse_i = n_singularity - i + 1;

    if (bosonic) {
      p[reverse_i] = sqrt(-gsl_sf_pow_int(2 * i, 2) * M_PI_POW_2
        - y_squared);
    } else {
      p[reverse_i] = sqrt(-gsl_sf_pow_int(2 * i - 1, 2) * M_PI_POW_2
        - y_squared);
    }

    #ifdef DEBUG
      printf("i = %d. p[%d] = %e\n", i, reverse_i, p[reverse_i]);
    #endif
  }

  // Add limits of integration
  p[0] = 0.;
  p[n_singularity + 1] = y;

  #ifdef DEBUG
    printf("integrand_points:\n");
    for (int i=0; i <= n_singularity + 1; i+=1) {
      printf("p[%d] = %e\n", i, p[i]);
    }
  #endif

  return p;
}

double J_quad(double y_squared, double abs_error, double rel_error,
              int max_n, bool bosonic) {
  /**
      @returns Numerical integration for Curtin eq. (2.12)

      Break integral into two domains: \f$0\f$ to \f$\Im y \f$, \f$\Im y\f$ to infinity.
      The first domain (if non empty) may contain singularities.
      The locations of the singularities are calculated and supplied to the
      integration routine.

      @warning For bosonic case, there is a discontinuity at the boundary,
      \f$x = \Im y \f$ if \f$y^2 \le 0\f$. This isn't treated explicitly in
      either domain, though doesn't appear to be problematic.

      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of subdivisions for memory allocation
      @param bosonic Whether bosonic (or fermionic) function required
  */
  gsl_set_error_handler_off();  // Default handler aborts
  gsl_function integrand;
  integrand.function = &J_integrand_wrapper;
  struct J_integrand_params params = {y_squared, bosonic};
  integrand.params = &params;
  double d1_result;
  double d2_result;
  double error;
  gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(max_n);

  if (y_squared >= 0.) {
    d1_result = 0.;
  } else {
    gsl_integration_qagp(&integrand,
                         integrand_points(y_squared, bosonic),
                         n_integrand_points(y_squared, bosonic),
                         abs_error,
                         rel_error,
                         max_n,
                         workspace,
                         &d1_result,
                         &error);
  }

  const double imag_y = imag(sqrt(cdouble(y_squared)));
  gsl_integration_qagiu(&integrand,
                        imag_y,
                        abs_error,
                        rel_error,
                        max_n,
                        workspace,
                        &d2_result,
                        &error);

  gsl_integration_workspace_free(workspace);

  #ifdef DEBUG
    printf("d1 = %e\nd2 = %e:\n", d1_result, d2_result);
  #endif

  return d1_result + d2_result;
}

double J_B_quad(double y_squared, double abs_error, double rel_error,
                int max_n) {
  /**
      @returns Bosonic thermal function found by numerical integration
      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of subdivisions for memory allocation
  */
  const double integral = J_quad(y_squared, abs_error, rel_error, max_n, true);
  #ifdef THROW
    const double bound = J_B_lim(y_squared);
    if (y_squared < neg_y_squared && std::abs(integral) > bound) {
      printf("quad = %e > bound = %e\n", integral, bound);
      throw std::runtime_error("quad exceeds upper bound");
    }
  #endif
  return integral;
}

double J_F_quad(double y_squared, double abs_error, double rel_error,
                int max_n) {
  /**
      @returns Fermionic thermal function found by numerical integration
      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of subdivisions for memory allocation
  */
  const double integral = J_quad(y_squared, abs_error, rel_error, max_n, false);
  #ifdef THROW
    const double bound = J_F_lim(y_squared);
    if (y_squared < neg_y_squared && std::abs(integral) > bound) {
      printf("quad = %e > bound = %e\n", integral, bound);
      throw std::runtime_error("quad exceeds upper bound");
    }
  #endif
  return integral;
}


// Thermal functions by Taylor expansion

/*! Factor that appears in terms in sum */
const double prefactor = 2. * gsl_sf_gamma(1.5) * gsl_sf_pow_int(4., -3) /
  (gsl_sf_fact(3) * M_PI_POW_2 * M_SQRTPI);
/*! \f$\zeta(3)\f$ */
const double zeta_3 = gsl_sf_zeta_int(3);

double gamma_sum(double y_squared, double abs_error, double rel_error,
                 int max_n, const bool bosonic, double sum = 0.) {
  /**
      @returns Sum of gamma functions in Wainwright eq. (2.18) and eq. (2.19)

      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of terms in sum
      @param bosonic Whether bosonic (or fermionic) function required
      @param sum Taylor expansion without gamma functions
  */
  #ifdef THROW
    if (std::abs(y_squared) >= 1.) {
      throw std::invalid_argument("|y_squared| >= 1."
                                   " - Taylor expansion invalid");
    }
  #endif

  double factor = gsl_sf_pow_int(y_squared, 3) * prefactor;
  double zeta = zeta_3;
  double term = factor * zeta;
  sum += term;

  for (int n = 2; n <= max_n; n += 1) {
    factor *= - y_squared * n / (n + 2.) * 0.25 * M_PI_POW_M2;
    zeta = gsl_sf_zeta_int(2 * n + 1);

    if (bosonic) {
      term = zeta * factor;
    } else {
      term = zeta * factor * (-1. + 0.125 * gsl_sf_pow_int(0.25, n + 2));
    }

    sum += term;

    if (std::abs(term) < std::max(abs_error, rel_error * std::abs(sum))) {
      #ifdef DEBUG
        printf("number of terms in sum = %d\n", n);
      #endif
      return sum;
    } else if (std::isinf(sum)) {
      #ifdef THROW
        throw std::runtime_error("sum diverging");
      #endif
      return sum;
    }

    #ifdef DEBUG
    if (n == max_n) {
      printf("reached max number of terms in sum = %d\n", n);
    }
    #endif
  }

  #ifdef DEBUG
    printf("gamma sum = %e\n", sum);
  #endif

  return sum;
}

double J_B_taylor(double y_squared, double abs_error, double rel_error,
                  int max_n) {
  /**
      @returns Bosonic thermal function from a Taylor expansion of
      Curtin eq. (2.12) as in Curtin eq. (2.13)

      @warning Valid for \f$|y^2| \ll 1\f$
      
      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of terms in sum
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_B_0;
  } else if (std::abs(y_squared) >= 1.) {
    #ifdef THROW
      throw std::invalid_argument("|y_squared| >= 1."
                                   "- Taylor expansion invalid");
    #endif
  }

  double real_y_cubed;

  if (y_squared >= 0.) {
    real_y_cubed = y_squared * sqrt(y_squared);
  } else {
    real_y_cubed = std::abs(y_squared) * sqrt(std::abs(y_squared)) * M_SQRT1_2;
  }

  double taylor_sum = - M_PI_POW_4 / 45.
                      + M_PI_POW_2 / 12. * y_squared
                      - M_PI / 6. * real_y_cubed
                      - gsl_sf_pow_int(y_squared, 2)
                        * log(std::abs(y_squared) / a_b) / 32.;

  const double sum = gamma_sum(y_squared, abs_error, rel_error, max_n, true,
    taylor_sum);

  #ifdef DEBUG
    const double gamma_sum_ = sum - taylor_sum;
    printf("gamma = %e and taylor = %e\n", gamma_sum_, taylor_sum);
  #endif

  return sum;
}

double J_F_taylor(double y_squared, double abs_error, double rel_error,
                  int max_n) {
  /**
      @returns Fermionic thermal function from a Taylor expansion of
      Curtin eq. (2.12) as in Curtin eq. (2.13)

      @warning Valid for \f$|y^2| \ll 1\f$
      
      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of terms in sum
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_F_0;
  } else if (std::abs(y_squared) >= 1.) {
    #ifdef THROW
      throw std::invalid_argument("|y_squared| >= 1."
                                   " - Taylor expansion invalid");
    #endif
  }

  double taylor_sum = 7. / 360. * M_PI_POW_4
                      - M_PI_POW_2 / 24. * y_squared
                      - gsl_sf_pow_int(y_squared, 2)
                        * log(std::abs(y_squared) / a_f) / 32.;

  const double sum = gamma_sum(y_squared, abs_error, rel_error, max_n, false,
    taylor_sum);

  #ifdef DEBUG
    const double gamma_sum_ = sum - taylor_sum;
    printf("gamma = %e and taylor = %e\n", gamma_sum_, taylor_sum);
  #endif

  return sum;
}


// Thermal functions by infinite sum of Bessel functions


double K2(cdouble x, bool fast = false) {
  /**
      @returns \f$K_2\f$ Bessel function
      
      Utilize fact that
      \f[\Re K_2(x i) = \frac{\pi}{2} Y_2(x)\f]
      to define \f$K_2\f$ for imaginary arguments.
      
      @param x
      @param fast Whether to use fast approximations
  */
  #ifdef THROW
  if (real(x) != 0. && imag(x) != 0.) {
    throw std::invalid_argument("K2 only implemented for x real or imaginary");
  } else  if (std::abs(x) == 0.) {
    throw std::invalid_argument("K2 diverges for |x| = 0");
  } else
  #endif
  if (real(x) != 0.) {
    if (fast) {
      // This is an asymptotic approximation for K_2 Bessel function
      return M_SQRTPI / M_SQRT2 / sqrt(real(x)) * exp(-real(x));
    } else {
      return gsl_sf_bessel_Kn(2, real(x));
    }
  } else if (imag(x) != 0.) {
      if (fast) {
        // This is an asymptotic approximation for Y_2 Bessel function
        return M_SQRTPI / M_SQRT2 / sqrt(imag(x)) * sin(M_PI_4 - imag(x));
      } else {
        return M_PI_2 * gsl_sf_bessel_Yn(2, imag(x));
      }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

double bessel_sum(double y_squared, double abs_error, double rel_error,
                  int max_n, bool fast, bool bosonic) {
  /**
      @returns Thermal functions from sum of Bessel functions

      Curtin eq. (2.14). 
      
      @warning Converges fastest (i.e. fewer terms required)
      if \f$|y^2| \gg 1\f$, though valid for all \f$y^2\f$ except
      \f$y^2 = 0\f$.
      
      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of terms in sum
      @param fast Whether to use fast approximations
      @param bosonic Whether bosonic (or fermionic) function required
  */
  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  const cdouble y = sqrt(cdouble(y_squared));
  const double sign = 2. * static_cast<double>(bosonic) - 1.;
  double factor = - y_squared * sign;

  gsl_set_error_handler_off();  // Default handler aborts
  double sum = factor * K2(y, fast);

  for (int n = 2; n <= max_n; n += 1) {
    const double n_double = static_cast<double>(n);
    factor *= sign * gsl_sf_pow_int((n_double - 1.) / n_double, 2);
    const double term = factor * K2(n_double * y, fast);
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

double J_F_bessel(double y_squared, double abs_error, double rel_error,
                  int max_n, bool fast) {
  /**
      @returns Fermionic thermal functions from sum of Bessel functions

      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of terms in sum
      @param fast Whether to use fast approximations
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_F_0;
  }

  double sum = bessel_sum(y_squared, abs_error, rel_error, max_n, fast, false);
  #ifdef THROW
    if (y_squared < neg_y_squared && std::abs(sum) > J_F_lim(y_squared)) {
      throw std::runtime_error("bessel sum exceeds upper bound");
    }
  #endif

  return sum;
}

double J_B_bessel(double y_squared, double abs_error, double rel_error,
                  int max_n, bool fast) {
  /**
      @returns Bosonic thermal functions from sum of Bessel functions
      
      @param y_squared Argument of thermal function     
      @param abs_error Maximum absolute error
      @param rel_error Maximum relative error
      @param max_n Maximum number of terms in sum
      @param fast Whether to use fast approximations
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_B_0;
  }

  double sum = bessel_sum(y_squared, abs_error, rel_error, max_n, fast, true);
  #ifdef THROW
    if (y_squared < neg_y_squared && std::abs(sum) > J_B_lim(y_squared)) {
      throw std::runtime_error("bessel sum exceeds upper bound");
    }
  #endif

  return sum;
}


// Limits and approximations


/*! Maximum of \f$\zeta_{-3/2}(x)\f$ from Mathematica */
constexpr double zeta_maximum = 0.024145376807240715;
/*! Minimum of \f$\zeta_{-3/2}(x)\f$ from Mathematica */
constexpr double zeta_minimum = -0.03154228985099239;

double J_B_lim(double y_squared, bool upper) {
  /**
      @returns Bound for bosonic thermal function
      
      @warning Valid for \f$y^2 \ll 0\f$
      
      @param y_squared Argument of thermal function     
      @param upper Whether to return upper (or lower) bound
  */
  #ifdef DEBUG
    if (y_squared > neg_y_squared) {
      printf("limit applicable for y_squared << 0. only\n");
    }
  #endif

  const double y = sqrt(std::abs(y_squared));

  if (upper) {
    return -y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta_minimum;
  } else {
    return -y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta_maximum;
  }
}

double J_F_lim(double y_squared, bool upper) {
  /**
      @returns Bound for fermionic thermal function
      
      @warning Valid for \f$y^2 \ll 0\f$

      @param y_squared Argument of thermal function     
      @param upper Whether to return upper (or lower) bound
  */
  return J_B_lim(y_squared, upper);
}

double J_F_approx(double y_squared) {
  /**
      @returns Approximation for fermionic thermal function
      
      @warning Valid for \f$y^2 \ll 0\f$
      
      @param y_squared Argument of thermal function     
  */
  double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for y_squared << 0. only\n");
      }
    #endif
    return M_SQRTPI / M_SQRT2 * y * sqrt(y) * sin(y - M_PI_4);
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for y_squared >> 0. only\n");
      }
    #endif
    return M_SQRTPI / M_SQRT2 * y * sqrt(y) * exp(-y);
  }
}

double J_B_approx(double y_squared) {
  /**
      @returns Approximation for bosonic thermal function
      
      @warning Valid for \f$y^2 \ll 0\f$
      
      @param y_squared Argument of thermal function  
  */
  return -J_F_approx(y_squared);
}

double shift_F(double y) {
  /**
      @returns Argument shifted into required domain
      @param y
  */
  double y_shift = fmod(y, 2. * M_PI);
  if (y_shift > M_PI) {
    y_shift -= 2. * M_PI;
  }
  return y_shift;
}

double shift_B(double y) {
  /**
      @returns Argument shifted into required domain
      @param y
  */
  return fmod(y, 2. * M_PI) - 2. * M_PI;
}

double J_F_zeta(double y_squared, int max_n) {
  /**
      @returns Fermionic thermal function from zeta function
      
      @warning Valid for \f$|y^2| \gg 1\f$
      
      @param y_squared Argument of thermal function 
      @param max_n Maximum number of terms in sum 
  */
  const double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for |y_squared| >> 1. only\n");
      }
    #endif
    const double zeta = std::real(hurwitz_zeta(-1.5, 0.5 - shift_F(y) * 0.5 * M_1_PI, max_n));
    return - y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta;
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for |y_squared| >> 1. only\n");
      }
    #endif
    const double poly = -gsl_sf_fermi_dirac_3half(-y);
    return - M_SQRTPI / M_SQRT2 * y * sqrt(y) * poly;
  }
}

double J_B_zeta(double y_squared, int max_n) {
  /**
      @returns Bosonic thermal function from zeta function

      @warning Valid for \f$|y^2| \gg 1\f$

      @param y_squared Argument of thermal function 
      @param max_n Maximum number of terms in sum 
  */
  const double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for |y_squared| >> 1. only\n");
      }
    #endif
    const double zeta = std::real(hurwitz_zeta(-1.5, - shift_B(y)
      / (2. * M_PI), max_n));
    return - y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta;
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for |y_squared| >> 1. only\n");
      }
    #endif
    const double poly = std::real(polylog(2.5, exp(-y), max_n));
    return - M_SQRTPI / M_SQRT2 * y * sqrt(y) * poly;
  }
}
