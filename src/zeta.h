/*

  Hurwitz Zeta function and polylogarithm.

  http://fredrikj.net/math/hurwitz_zeta.pdf

*/

#ifndef _ZETA_H_
#define _ZETA_H_

#include <complex>

typedef std::complex<double> cdouble;

#ifdef __cplusplus
extern "C" {
#endif

cdouble hurwitz_zeta(double s, cdouble a, int N = 50);
cdouble polylog(double s, cdouble a, int N = 50);

#ifdef __cplusplus
}
#endif

#endif  // _ZETA_H_
