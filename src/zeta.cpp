/*

  Hurwitz Zeta function.

  http://fredrikj.net/math/hurwitz_zeta.pdf

*/

#include <zeta.h>
#include <bernoulli.h>

#include <gsl/gsl_sf.h>
#include <complex>


cdouble S(double s, cdouble a, int N) {
  cdouble sum = 0.;

  for (int k = 0; k <= N - 1; k += 1) {
    sum += pow(a + static_cast<cdouble>(k), -s);
  }

  return sum;
}

cdouble I(double s, cdouble a, int N) {
  return pow(a + static_cast<cdouble>(N), 1. - s) / (s - 1.);
}

cdouble T(double s, cdouble a, int N, int M) {
  const cdouble factor = pow(a + static_cast<cdouble>(N), -s);

  if (M > B_2n_fact_size) {
    #ifdef DEBUG
      printf("M = %d > B_2n_fact_size = %d. Using M = B_2n_fact_size\n");
    #endif

    #ifdef THROW
      throw std::invalid_argument("Bernoulli numbers out of bounds");
    #endif

    M = B_2n_fact_size;
  }

  cdouble sum = 0.;

  for (int k = 1; k <= M; k += 1) {
    sum += B_2n_fact[k] * gsl_sf_gamma(s + 2. * k - 1.) / gsl_sf_gamma(s) / pow(a + static_cast<cdouble>(N), 2. * k - 1.);
  }

  return factor * (0.5 + sum);
}

cdouble hurwitz_zeta(double s, cdouble a, int N) {
  if (N > B_2n_fact_size) {
    #ifdef DEBUG
      printf("N = %d > B_2n_fact_size = %d. Using N = B_2n_fact_size\n");
    #endif

    #ifdef THROW
      throw std::invalid_argument("Bernoulli numbers out of bounds");
    #endif

    N = B_2n_fact_size;
  }

  return S(s, a, N) + I(s, a, N) + T(s, a, N, N);
}

cdouble polylog(double s, cdouble a, int N) {
  const cdouble j(0., 1.);
  const cdouble factor = pow(0.5 * j / M_PI, 1. - s) * gsl_sf_gamma(1. - s);
  const cdouble arg = 0.5 * j * log(a) / M_PI;
  const cdouble zeta = - pow(j, 2. * s) * hurwitz_zeta(1. - s, arg, N) + hurwitz_zeta(1. - s, 1. - arg, N);
  return factor * zeta;
}
