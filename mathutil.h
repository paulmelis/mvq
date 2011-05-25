#ifndef MVQ_UTIL_H
#define MVQ_UTIL_H

#ifdef WIN32
// For M_PI and friends
#define _USE_MATH_DEFINES
#endif
#include <cmath>

#ifndef M_PI
#define M_PI 3.141592653589793238462643350279
#endif

namespace mvq {

const double epsilon=1.0e-8;
const double infinity=1.0e8;

inline bool
equal(double d1, double d2, double eps=epsilon)
{
	return (fabs(d1-d2) < eps);
}

inline bool
zero(double d, double eps=epsilon)
{
	return (fabs(d) < eps);
}

inline double
deg2rad(double d)
{
    return d/180.0*M_PI;
}

inline double
rad2deg(double r)
{
    return r/M_PI*180.0;
}

template<typename T>
inline T
min(T a, T b)
{
    if (a <= b)
        return a;
    else
        return b;
}
template<typename T>
inline T
max(T a, T b)
{
    if (a >= b)
        return a;
    else
        return b;
}

}

#endif
