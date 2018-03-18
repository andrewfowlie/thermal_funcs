/*

  Hurwitz Zeta function and polylogarithm.

  http://fredrikj.net/math/hurwitz_zeta.pdf

*/

#include <zeta.h>
#include <bernoulli.h>

#include <gsl/gsl_sf.h>
#include <complex>
#include <stdexcept>


cdouble S(double s, cdouble a, int N) {
  /**
      @returns S part of zeta function.
  */
  cdouble sum = 0.;

  for (int k = 0; k <= N - 1; k += 1) {
    sum += pow(a + static_cast<cdouble>(k), -s);
  }

  return sum;
}

cdouble I(double s, cdouble a, int N) {
  /**
      @returns I part of zeta function.
  */
  return pow(a + static_cast<cdouble>(N), 1. - s) / (s - 1.);
}

cdouble T(double s, cdouble a, int N, int M) {
  /**
      @returns T part of zeta function.
  */
  const cdouble d = a + static_cast<cdouble>(N);
  const cdouble factor = pow(d, -s);

  if (M > B_2n_fact_size) {
    #ifdef DEBUG
      printf("M = %d > B_2n_fact_size = %d. Using M = B_2n_fact_size\n", M, B_2n_fact_size);
    #endif

    #ifdef THROW
      throw std::invalid_argument("Bernoulli numbers out of bounds");
    #endif

    M = B_2n_fact_size;
  }

  cdouble sum = 0.;

  for (int k = 1; k <= M; k += 1) {
    sum += B_2n_fact[k] * gsl_sf_poch(s, 2. * k - 1.) / pow(d, 2 * k - 1);
  }

  return factor * (0.5 + sum);
}

cdouble hurwitz_zeta(double s, cdouble a, int N) {
  /**
      @returns Hurwitz zeta function.
  */
  if (N > B_2n_fact_size) {
    #ifdef DEBUG
      printf("N = %d > B_2n_fact_size = %d. Using N = B_2n_fact_size\n", N, B_2n_fact_size);
    #endif

    #ifdef THROW
      throw std::invalid_argument("Bernoulli numbers out of bounds");
    #endif

    N = B_2n_fact_size;
  }

  return S(s, a, N) + I(s, a, N) + T(s, a, N, N);
}

cdouble polylog(double s, cdouble a, int N) {
  /**
      @returns Polylogarithm function.
  */
  if (s > 1 && std::abs(a) <= 1) {
    cdouble term = a;
    cdouble sum = term;
    for (int i = 2; i <= N; i += 1) {
      term *= a * pow((i - 1) / i, 2.5);
      sum += term;
    }
    return sum;
  }

  const cdouble j(0., 1.);
  const cdouble factor = pow(0.5 * j * M_1_PI, 1. - s) * gsl_sf_gamma(1. - s);
  const cdouble arg = 0.5 * j * log(a) * M_1_PI;
  const cdouble zeta = - pow(j, 2. * s) * hurwitz_zeta(1. - s, arg, N) +
                       hurwitz_zeta(1. - s, 1. - arg, N);
  return factor * zeta;
}
