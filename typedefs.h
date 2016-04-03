#ifndef TYPEDEFS
#define TYPEDEFS

#include <iostream>
#include <cmath>
#include <limits>
#include <iomanip>
#include <cstdlib>

using namespace std;
//  define precision by commenting out one of the two lines:
typedef float reals;       //  defines reals as double (standard for scientific calculations)
const reals One=1.0f,Two=2.0f,Three=3.0f,Four=4.0f,Five=5.0f,Six=6.0f,Ten=10.0f;
const reals Pi=3.141592653589793238462643383f;
const reals REAL_MAX=numeric_limits<reals>::max();
const reals REAL_MIN=numeric_limits<reals>::min();
const reals REAL_EPSILON=numeric_limits<reals>::epsilon();


#endif // TYPEDEFS

