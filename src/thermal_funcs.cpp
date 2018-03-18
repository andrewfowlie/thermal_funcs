/*

  Implementation of thermal functions in finite-temperature field theory. I
  refer throughout to Thermal Resummation and Phase Transitions by Curtin et al:

  https://arxiv.org/pdf/1612.00466.pdf

  The thermal functions are eq. (2.12):

  J_{B/F}(y^2) = \int_0^\infnty dx x^2 \log(1 \mp \exp(-\sqrt{x^2 + y^2}))

  We calculate only the real part.

  I also refer to Phase Transitions in the Early Universe by Wainwright:

  http://escholarship.org/uc/item/2r84s7h9

  for eq. (2.18) and (2.19).

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


constexpr double neg_y_squared = -1.E3;
constexpr double pos_y_squared = 1.E3;
constexpr double a_b = pow(M_PI, 2) * exp(1.5 - 2. * M_EULER);  // Below eq. (2.13)
constexpr double a_f = 16. * a_b;  // Below eq. (2.13)
constexpr double M_PI_POW_2 = pow(M_PI, 2);
constexpr double M_PI_POW_4 = pow(M_PI, 4);


// Thermal functions at y_squared -> 0. Found in Mathematica:
//
// -Sum[1/n^2 * Limit[x^2 * BesselK[2, x *n], x -> 0], {n, 1, Infinity}]
// -Sum[(-1)^n/n^2 * Limit[x^2 * BesselK[2, x *n], x -> 0], {n, 1, Infinity}]
//
// or from Taylor expansions.

constexpr double J_B_0 = - pow(M_PI, 4) / 45.;
constexpr double J_F_0 = 7. * pow(M_PI, 4) / 360.;


// Thermal functions by numerical integration


double J_integrand(double x, double y_squared, bool bosonic) {
  /**
     @returns Integrand in Curtin eq. (2.12)
     x^2 Log[1 -+ Exp[-Sqrt[x^2 + y^2]]]
     If r^2 = x^2 + y^2 < 0, this can be written
     x^2 * 0.5 * Log[2. + 2. * Cos[r]]
  */
  const double sign = 1. - 2. * static_cast<double>(bosonic);
  const double r_squared = gsl_sf_pow_int(x, 2) + y_squared;
  const double abs_r = sqrt(std::abs(r_squared));
  if (r_squared >= 0.) {
    return gsl_sf_pow_int(x, 2) * gsl_log1p(sign * exp(-abs_r));
  } else {
    return 0.5 * gsl_sf_pow_int(x, 2) * (gsl_log1p(sign * cos(abs_r)) + M_LN2);
  }
}

struct J_integrand_params {
  double y_squared;
  bool bosonic;
};

double J_integrand_wrapper(double x, void *p) {
  /**
     Wrapper for J_integrand with signature required by numerical integration.
  */
  const struct J_integrand_params *params = (struct J_integrand_params *)p;
  const double y_squared = params->y_squared;
  const bool bosonic = params->bosonic;
  return J_integrand(x, y_squared, bosonic);
}

int n_integrand_points(double y_squared, const bool bosonic) {
  /**
     @returns Number of integrand points, i.e. number of singularities + 2 for
     endpoints of integration.

     NB don't include an endpoint twice if it is singular. This means we exclude
     n = 0.
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
     @returns Singularities in integrand, present if y_squared < 0.,
     and boundaries of integration, in ascending order.

     NB don't include an endpoint twice if it is singular.

     Singularities occur at
     0. <= x = Sqrt[-n^2 Pi^2 - y_squared] <= |y|
     for n even/odd for bosonic/fermionic.
  */
  #ifdef THROW
    if (y_squared >= 0.) {
      throw std::invalid_argument("|y_squared| >= 0. - no singularities possible");
    }
  #endif

  const double y = sqrt(std::abs(y_squared));
  const int n_singularity = n_integrand_points(y_squared, bosonic) - 2;
  double *p = reinterpret_cast<double *>(malloc(sizeof(double) * (n_singularity + 2)));

  for (int i = 1; i <= n_singularity; i += 1) {
    // Insure result is in ascending order
    const int reverse_i = n_singularity - i + 1;

    if (bosonic) {
      p[reverse_i] = sqrt(-gsl_sf_pow_int(2 * i, 2) * M_PI_POW_2 - y_squared);
    } else {
      p[reverse_i] = sqrt(-gsl_sf_pow_int(2 * i - 1, 2) * M_PI_POW_2 - y_squared);
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

double J_quad(double y_squared, double abs_error, double rel_error, int max_n,
              bool bosonic) {
  /**
      @returns Numerical integration for Curtin eq. (2.12)

      Break integral into two domains: 0 to Im(y), Im(y) to infinity.
      The first domain (if non empty) may contain singularities.
      The locations of the singularities are calculated and supplied to the
      integration routine.

      NB For bosonic case, there is a discontinuity at the boundary, x = Im(y)
      if y_squared <= 0. This isn't treated explicity in either domain, though
      doesn't appear to be problematic.
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
     @returns Bosonic thermal function found by numerical integration.
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
     @returns Fermionic thermal function found by numerical integration.
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


double gamma_sum(double y_squared, double abs_error, double rel_error,
                 int max_n, const bool bosonic, double sum = 0.) {
  /**
     @returns Sum of Gamma functions in Wainwright eq. (2.18) and eq. (2.19).
  */
  #ifdef THROW
    if (std::abs(y_squared) >= 1.) {
      throw std::invalid_argument("|y_squared| >= 1. - Taylor expansion invalid");
    }
  #endif

  double factor = 2. * gsl_sf_gamma(1.5) * gsl_sf_pow_int(y_squared, 3) * gsl_sf_pow_int(4., -3) /
                  (gsl_sf_fact(3) * M_PI_POW_2 * M_SQRTPI);
  double zeta = gsl_sf_zeta_int(3);
  double term = factor * zeta;
  sum += term;

  for (int n = 2; n <= max_n; n += 1) {
    factor *= - y_squared * n / (n + 2.) * 0.25 / M_PI_POW_2;
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
      Curtin eq. (2.12) as in Curtin eq. (2.13).

      Valid for |y_squared| << 1.
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_B_0;
  } else if (std::abs(y_squared) >= 1.) {
    #ifdef THROW
      throw std::invalid_argument("|y_squared| >= 1. - Taylor expansion invalid");
    #endif
  }

  const double real_y_cubed = std::abs(real(pow(cdouble(y_squared), 1.5)));

  double taylor_sum = - M_PI_POW_4 / 45.
                      + M_PI_POW_2 / 12. * y_squared
                      - M_PI / 6. * real_y_cubed
                      - gsl_sf_pow_int(y_squared, 2) * log(std::abs(y_squared) / a_b) / 32.;

  const double sum = gamma_sum(y_squared, abs_error, rel_error, max_n, true, taylor_sum);

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
      Curtin eq. (2.12) as in Curtin eq. (2.13).

      Valid for |y_squared| << 1.
  */

  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_F_0;
  } else if (std::abs(y_squared) >= 1.) {
    #ifdef THROW
      throw std::invalid_argument("|y_squared| >= 1. - Taylor expansion invalid");
    #endif
  }

  double taylor_sum = 7. / 360. * M_PI_POW_4
                      - M_PI_POW_2 / 24. * y_squared
                      - gsl_sf_pow_int(y_squared, 2) * log(std::abs(y_squared) / a_f) / 32.;

  const double sum = gamma_sum(y_squared, abs_error, rel_error, max_n, false, taylor_sum);

  #ifdef DEBUG
    const double gamma_sum_ = sum - taylor_sum;
    printf("gamma = %e and taylor = %e\n", gamma_sum_, taylor_sum);
  #endif

  return sum;
}


// Thermal functions by infinite sum of Bessel functions


double K2(cdouble x, bool fast = false) {
  /**
      @returns K2 Bessel function.
      Utilize fact that
      Re[BesselK[2, x * I]] = 0.5 * Pi BesselY[2, x]
      to define K2 for imaginary arguments.
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
}

double bessel_sum(double y_squared, double abs_error, double rel_error,
                  int max_n, bool fast, bool bosonic) {
  /**
      @returns Thermal functions from sum of Bessel functions.

      Curtin eq. (2.14). Converges fastest (i.e. fewer terms required)
      if |y_squared| >> 1, though valid for all |y_squared| except
      y_squared = 0.
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
      @returns Fermionic thermal functions from sum of Bessel functions.
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
      @returns Bosonic thermal functions from sum of Bessel functions.
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


// Maxima and minima of Zeta[-3/2, x] from Mathematica

constexpr double zeta_maxima = 0.024145376807240715;
constexpr double zeta_minima = -0.03154228985099239;

double J_B_lim(double y_squared, bool upper) {
  /**
      @returns Bound for bosonic thermal function.
      Applicable for y_squared << 0.
  */
  #ifdef DEBUG
    if (y_squared > neg_y_squared) {
      printf("limit applicable for y_squared << 0. only\n");
    }
  #endif

  const double y = sqrt(std::abs(y_squared));

  if (upper) {
    return -y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta_minima;
  } else {
    return -y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta_maxima;
  }
}

double J_F_lim(double y_squared, bool upper) {
  /**
      @returns Bound for fermionic thermal function.
      Applicable for y_squared << 0.
  */
  return J_B_lim(y_squared, upper);
}

double J_F_approx(double y_squared) {
  /**
      @returns Approximation for fermionic thermal function.
      Applicable for y_squared << 0.
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
      @returns Approximation for bosonic thermal function.
      Applicable for y_squared << 0.
  */
  return -J_F_approx(y_squared);
}

double shift_F(double y) {
  /**
      @returns Argument shifted into required domain.
  */
  double y_shift = fmod(y, 2. * M_PI);
  if (y_shift > M_PI) {
    y_shift -= 2. * M_PI;
  }
  return y_shift;
}

double shift_B(double y) {
  /**
      @returns Argument shifted into required domain.
  */
  return fmod(y, 2. * M_PI) - 2. * M_PI;
}

double J_F_zeta(double y_squared, int max_n) {
  /**
      @returns Fermionic thermal function from Zeta function.
  */
  const double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for y_squared << 0. only\n");
      }
    #endif
    const double zeta = std::real(hurwitz_zeta(-1.5, 0.5 - shift_F(y) * 0.5 * M_1_PI, max_n));
    return - y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta;
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for y_squared >> 0. only\n");
      }
    #endif
    const double poly = -gsl_sf_fermi_dirac_3half(-y);
    return - M_SQRTPI / M_SQRT2 * y * sqrt(y) * poly;
  }
}

double J_B_zeta(double y_squared, int max_n) {
  /**
      @returns Bosonic thermal function from Zeta function.
  */
  const double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for y_squared << 0. only\n");
      }
    #endif
    const double zeta = std::real(hurwitz_zeta(-1.5, - shift_B(y) / (2. * M_PI), max_n));
    return - y * sqrt(y) * 8. / 3. * M_PI_POW_2 * M_SQRTPI * zeta;
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for y_squared >> 0. only\n");
      }
    #endif
    const double poly = std::real(polylog(2.5, exp(-y), max_n));
    return - M_SQRTPI / M_SQRT2 * y * sqrt(y) * poly;
  }
}
