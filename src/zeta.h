/*

*/

#ifndef _ZETA_H_
#define _ZETA_H_

#include <complex>

typedef std::complex<double> cdouble;

extern "C" cdouble hurwitz_zeta(double s, cdouble a, int N = 50);
extern "C" cdouble polylog(double s, double a, int N = 50);

#endif  // _ZETA_H_
