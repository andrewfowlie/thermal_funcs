/*

  Implementation of thermal functions in finite-temperature field theory. I refer
  throughout to Thermal Resummation and Phase Transitions by Curtin et al:

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
#include <stdexcept>
#include <complex>
#include <algorithm>

const float neg_y_squared = -1.E3;
const float pos_y_squared = 1.E3;
const float a_b = pow(M_PI, 2) * exp(1.5 - 2. * M_EULER);  // Below eq. (2.13)
const float a_f = 16. * a_b;  // Below eq. (2.13)

// Thermal functions at y_squared -> 0. Found in Mathematica:
//
// -Sum[1/n^2 * Limit[x^2 * BesselK[2, x *n], x -> 0], {n, 1, Infinity}]
// -Sum[(-1)^n/n^2 * Limit[x^2 * BesselK[2, x *n], x -> 0], {n, 1, Infinity}]
//
// or from Taylor expansions.

const float J_B_0 = - pow(M_PI, 4) / 45.;
const float J_F_0 = 7. * pow(M_PI, 4) / 360.;

// Thermal functions by numerical integration


double J_integrand(double x, double y_squared, bool bosonic) {
  // Integrand in Curtin eq. (2.12)
  // x^2 Log[1 -+ Exp[-Sqrt[x^2 + y^2]]]
  // If r^2 = x^2 + y^2 < 0, this can be written
  // x^2 * 0.5 * Log[2. + 2. * Cos[r]]
  double sign = 1. - 2. * static_cast<double>(bosonic);
  double r_squared = pow(x, 2) + y_squared;
  double abs_r = sqrt(std::abs(r_squared));
  if (r_squared >= 0.) {
    return pow(x, 2) * gsl_log1p(sign * exp(-abs_r));
  } else {
    return 0.5 * pow(x, 2) * (gsl_log1p(sign * cos(abs_r)) + log(2.));
  }
}

struct J_integrand_params {double y_squared; bool bosonic;};

double J_integrand_wrapper(double x, void *p) {
  // Wrapper for J_integrand with signature required by numerical integration
  struct J_integrand_params *params = (struct J_integrand_params *)p;
  double y_squared = params->y_squared;
  bool bosonic = params->bosonic;
  return J_integrand(x, y_squared, bosonic);
}

int n_integrand_points(double y_squared, bool bosonic) {
  // Number of integrand points, i.e. number of singularities + 2 for endpoints
  // of integration. NB don't include an endpoint twice if it is singular. This
  // means we exclude n = 0.

  #ifdef THROW
    if (y_squared >= 0.) {
      throw std::invalid_argument("|y_squared| >= 0. - no singularities possible");
    }
  #endif

  double y = sqrt(std::abs(y_squared));
  int n_singularity;
  int max_n = floor(y / M_PI);

  if (bosonic) {
    n_singularity = static_cast<int>(floor(0.5 * max_n));  // Even numbers <= max_n (excluding 0)
  } else {
    n_singularity = static_cast<int>(ceil(0.5 * max_n));  // Odd numbers <= max_n
  }

  return n_singularity + 2;
}

double *integrand_points(double y_squared, bool bosonic) {
  // Singularities in integrand, present if y_squared < 0., and boundaries of
  // integration, in ascending order.
  // NB don't include an endpoint twice if it is singular.
  //
  // Singularities occur at
  // 0. <= x = Sqrt[-n^2 Pi^2 - y_squared] <= |y|
  // for n even/odd for bosonic/fermionic.

  #ifdef THROW
    if (y_squared >= 0.) {
      throw std::invalid_argument("|y_squared| >= 0. - no singularities possible");
    }
  #endif

  double y = sqrt(std::abs(y_squared));
  int n_singularity = n_integrand_points(y_squared, bosonic) - 2;
  double *p = reinterpret_cast<double *>(malloc(sizeof(double) * (n_singularity + 2)));
  int reverse_i;

  for (int i=1; i <= n_singularity; i+=1) {
    reverse_i = n_singularity - i + 1;  // Insure result is in ascending order
    if (bosonic) {
      p[reverse_i] = sqrt(-pow(2 * i, 2) * pow(M_PI, 2) - y_squared);
    } else {
      p[reverse_i] = sqrt(-pow(2 * i - 1, 2) * pow(M_PI, 2) - y_squared);
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

double J_quad(double y_squared, double abs_error, double rel_error, int max_n, bool bosonic) {
  // Numerical integration for Curtin eq. (2.12)
  // Break integral into two domains: 0 to Im(y), Im(y) to infinity.
  // The first domain (if non empty) may contain singularities.
  // The locations of the singularities are calculated and supplied to the
  // integration routine.
  //
  // NB For bosonic case, there is a discontinuity at the boundary, x = Im(y)
  // if y_squared <= 0. This isn't treated explicity in either domain, though doesn't
  // appear to be problematic.

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

  double imag_y = imag(sqrt(std::complex<double>(y_squared)));
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

double J_B_quad(double y_squared, double abs_error, double rel_error, int max_n) {
  double integral = J_quad(y_squared, abs_error, rel_error, max_n, true);
  #ifdef THROW
    double bound = J_B_lim(y_squared);
    if (y_squared < neg_y_squared && std::abs(integral) > bound) {
      printf("quad = %e > bound = %e\n", integral, bound);
      throw std::runtime_error("quad exceeds upper bound");
    }
  #endif
  return integral;
}

double J_F_quad(double y_squared, double abs_error, double rel_error, int max_n) {
  double integral = J_quad(y_squared, abs_error, rel_error, max_n, false);
  #ifdef THROW
    double bound = J_F_lim(y_squared);
    if (y_squared < neg_y_squared && std::abs(integral) > bound) {
      printf("quad = %e > bound = %e\n", integral, bound);
      throw std::runtime_error("quad exceeds upper bound");
    }
  #endif
  return integral;
}


// Thermal functions by Taylor expansion


double gamma_sum(double y_squared, double abs_error, double rel_error, int max_n, bool bosonic, double sum = 0.) {
  // Sum of Gamma functions in Wainwright eq. (2.18) and eq. (2.19).

  #ifdef THROW
    if (std::abs(y_squared) >= 1.) {
      throw std::invalid_argument("|y_squared| >= 1. - Taylor expansion invalid");
    }
  #endif

  double factor = 2. * pow(M_PI, -2.5) * pow(y_squared, 3) * pow(4., -3) / gsl_sf_fact(3) * gsl_sf_gamma(1.5);
  double zeta = gsl_sf_zeta_int(3);
  double term = factor * zeta;
  sum += term;

  for (int n = 2; n <= max_n; n += 1) {
    factor *= - y_squared * n / (n + 2.) * 0.25 / pow(M_PI, 2);
    zeta = gsl_sf_zeta_int(2 * n + 1);

    if (bosonic) {
      term = zeta * factor;
    } else {
      term = zeta * factor * (-1. + 0.125 * pow(0.25, n + 2));
    }

    sum += term;

    if (std::abs(term) < std::max(abs_error, rel_error * sum)) {
      #ifdef DEBUG
        printf("number of terms in sum = %d\n", n);
      #endif
      return sum;
    } else if (isinf(sum)) {
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

double J_B_taylor(double y_squared, double abs_error, double rel_error, int max_n) {
  // Taylor expansion of Curtin eq. (2.12) as in Curtin eq. (2.13). Valid for |y_squared| << 1.
  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_B_0;
  } else if (std::abs(y_squared) >= 1.) {
    #ifdef THROW
      throw std::invalid_argument("|y_squared| >= 1. - Taylor expansion invalid");
    #endif
  }

  double real_y_cubed = std::abs(real(pow(std::complex<double>(y_squared), 1.5)));

  double taylor_sum = - pow(M_PI, 4) / 45.
                      + pow(M_PI, 2) / 12. * y_squared
                      - M_PI / 6. * real_y_cubed
                      - pow(y_squared, 2) * log(std::abs(y_squared) / a_b) / 32.;

  double sum = gamma_sum(y_squared, abs_error, rel_error, max_n, true, taylor_sum);

  #ifdef DEBUG
    double gamma_sum_ = sum - taylor_sum;
    printf("gamma = %e and taylor = %e\n", gamma_sum_, taylor_sum);
  #endif

  return sum;
}

double J_F_taylor(double y_squared, double abs_error, double rel_error, int max_n) {
  // Taylor expansion of Curtin eq. (2.12) as in Curtin eq. (2.13). Valid for |y_squared| << 1.
  // If y_squared = 0, known limit returned.
  if (y_squared == 0.) {
    return J_F_0;
  } else if (std::abs(y_squared) >= 1.) {
    #ifdef THROW
      throw std::invalid_argument("|y_squared| >= 1. - Taylor expansion invalid");
    #endif
  }

  double taylor_sum = 7. / 360. * pow(M_PI, 4)
                      - pow(M_PI, 2) / 24. * y_squared
                      - pow(y_squared, 2) * log(std::abs(y_squared) / a_f) / 32.;

  double sum = gamma_sum(y_squared, abs_error, rel_error, max_n, false, taylor_sum);

  #ifdef DEBUG
    double gamma_sum_ = sum - taylor_sum;
    printf("gamma = %e and taylor = %e\n", gamma_sum_, taylor_sum);
  #endif

  return sum;
}


// Thermal functions by infinite sum of Bessel functions


double K2(std::complex<double> x, bool fast = false) {
  // Utilize fact that
  // Re[BesselK[2, x * I]] = 0.5 * Pi BesselY[2, x]
  // to define K2 for imaginary arguments.

  gsl_set_error_handler_off();  // Default handler aborts

  if (real(x) != 0. && imag(x) != 0.) {
    #ifdef THROW
      throw std::invalid_argument("K2 only implemented for x real or imaginary");
    #endif
  } else  if (std::abs(x) == 0.) {
    #ifdef THROW
      throw std::invalid_argument("K2 diverges for |x| = 0");
    #endif
  } else if (real(x) != 0.) {
    return gsl_sf_bessel_Kn(2, real(x));
  } else if (imag(x) != 0.) {
      if (fast) {
        // This is an asymptotic approximation for Y_2 Bessel function
        return sqrt(0.5 * M_PI / imag(x)) * sin(M_PI * 0.25 - imag(x));
      } else {
        return 0.5 * M_PI * gsl_sf_bessel_Yn(2, imag(x));
      }
  }
}

double bessel_sum(double y_squared, double abs_error, double rel_error, int max_n, bool fast, bool bosonic) {
  // Curtin eq. (2.14). Converges fastest (i.e. fewer terms required) if |y_squared| >> 1, though
  // valid for all |y_squared| except y_squared = 0.

  #ifdef THROW
    if (y_squared == 0.) {
      throw std::invalid_argument("y_squared == 0 invalid");
    }
  #endif

  std::complex<double> y = sqrt(std::complex<double>(y_squared));
  double sign = 2. * static_cast<double>(bosonic) - 1.;
  double factor = sign;
  double sum = factor * K2(y, fast);

  for (int n = 2; n <= max_n; n += 1) {
    factor *= sign * pow((static_cast<double>(n) - 1.) / static_cast<double>(n), 2);
    const double term = factor * K2(static_cast<double>(n) * y, fast);
    sum += term;

    #ifdef DEBUG
      printf("term_%d = %e and sum = %e\n", n, term, sum);
    #endif

    if (std::abs(term) < std::max(abs_error, rel_error * sum)) {
      #ifdef DEBUG
        printf("number of terms in sum = %d\n", n);
      #endif
      break;
    } else if (isinf(sum)) {
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

  return -y_squared * sum;
}

double J_F_bessel(double y_squared, double abs_error, double rel_error, int max_n, bool fast) {
  // If y_squared = 0, known limit returned, otherwise, J_F from sum of Bessel functions.
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

double J_B_bessel(double y_squared, double abs_error, double rel_error, int max_n, bool fast) {
  // If y_squared = 0, known limit returned, otherwise, J_B from sum of Bessel functions.
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
const double zeta_maxima = 0.024145376807240715;
const double zeta_minima = -0.03154228985099239;

double J_B_lim(double y_squared, bool upper) {
  // Found from Eq. (14) in http://mathworld.wolfram.com/BesselFunctionoftheSecondKind.html
  #ifdef DEBUG
    if (y_squared > neg_y_squared) {
      printf("limit applicable for y_squared << 0. only\n");
    }
  #endif

  double y = sqrt(std::abs(y_squared));
  double zeta;
  if (upper) {
    return -pow(y, 1.5) * 8. / 3. * pow(M_PI, 2.5) * zeta_minima;
  } else {
    return pow(y, 1.5) * 8. / 3. * pow(M_PI, 2.5) * zeta_maxima;
  }
}

double J_F_lim(double y_squared, bool upper) {
  return J_B_lim(y_squared, upper);
}

double J_F_approx(double y_squared) {
  double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for y_squared << 0. only\n");
      }
    #endif
    return sqrt(0.5 * M_PI) * pow(y, 1.5) * sin(y - M_PI * 0.25);
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for y_squared >> 0. only\n");
      }
    #endif
    return sqrt(0.5 * M_PI) * pow(y, 1.5) * exp(-y);
  }
}

double J_B_approx(double y_squared) {
  return -J_F_approx(y_squared);
}

double shift_F(double y) {
  double y_shift = fmod(y, 2. * M_PI);
  if (y_shift > M_PI) {
    y_shift -= 2. * M_PI;
  }
  return y_shift;
}

double shift_B(double y) {
  return fmod(y, 2. * M_PI) - 2. * M_PI;
}

double J_F_zeta(double y_squared, int max_n) {
  double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for y_squared << 0. only\n");
      }
    #endif
    const double zeta = std::real(hurwitz_zeta(-1.5, 0.5 - shift_F(y) / (2. * M_PI), max_n));
    return - pow(y, 1.5) * 8. / 3. * pow(M_PI, 2.5) * zeta;
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for y_squared >> 0. only\n");
      }
    #endif
    double poly = -gsl_sf_fermi_dirac_3half(-y);
    return - sqrt(0.5 * M_PI) * pow(y, 1.5) * poly;
  }
}

double J_B_zeta(double y_squared, int max_n) {
  double y = sqrt(std::abs(y_squared));
  if (y_squared < 0.) {
    #ifdef DEBUG
      if (y_squared > neg_y_squared) {
        printf("approx applicable for y_squared << 0. only\n");
      }
    #endif
    const double zeta = std::real(hurwitz_zeta(-1.5, - shift_B(y) / (2. * M_PI), max_n));
    return - pow(y, 1.5) * 8. / 3. * pow(M_PI, 2.5) * zeta;
  } else {
    #ifdef DEBUG
      if (y_squared < pos_y_squared) {
        printf("approx applicable for y_squared >> 0. only\n");
      }
    #endif
    const double poly = std::real(polylog(2.5, exp(-y), max_n));
    return -sqrt(0.5 * M_PI) * pow(y, 1.5) * poly;
  }
}
