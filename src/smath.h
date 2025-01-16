#ifndef __SMATH_H__
#define __SMATH_H__

const double EPSILON = .1192e-06;
const double INF = 1e32;

#ifdef __cplusplus
extern "C" {
#endif

// ln(gamma(z)), for z > 0
double lngamma(const double z);

// Tail area under the normal curve.
double alnorm(const double x, const bool upper);

// Incomplete Gamma integral
double gammad(const double x, const double p);

// Chi^2 distribution function
double chi_squared(const int ndf, const double chi2);


#ifdef __cplusplus
}
#endif
/*
template<typename T>
T min(const T a, const T b) { return a < b ? a : b; }

template<typename T>
T max(const T a, const T b) { return a > b ? a : b; }
*/
#endif // __SMATH_H__
