/*

  Hurwitz Zeta function.

  http://fredrikj.net/math/hurwitz_zeta.pdf

*/

#include <zeta.h>
#include <bernoulli.h>

#include <gsl/gsl_sf.h>
#include <complex>


double S(double s, double a, int N) {
  double sum = 0.;

  for (int k = 0; k <= N - 1; k += 1) {
    sum += pow(a + k, -s);
  }

  return sum;
}

double I(double s, double a, int N) {
  return pow(a + N, 1. - s) / (s - 1.);
}

double T(double s, double a, int N, int M) {
  const double k = pow(a + N, -s);

  double sum = 0.;

  for (int k = 1; k <= M; k += 1) {
    sum += B_2n_fact[k] * gsl_sf_gamma(s + 2. * k - 1.) / gsl_sf_gamma(s) / pow(a + N, 2. * k - 1.);
  }

  return k * (0.5 + sum);
}

double hurwitz_zeta(double s, double a, int N) {
  return S(s, a, N) + I(s, a, N) + T(s, a, N, N);
}
