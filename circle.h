#ifndef CIRCLE_H
#define CIRCLE_H

#include "typedefs.h"

// A circle has 7 fields:
//     a, b, r (of type reals), the circle parameters
//     s (of type reals), the estimate of sigma (standard deviation)
//     g (of type reals), the norm of the gradient of the objective function
//     i and j (of type int), the iteration counters (outer and inner, respectively)

class Circle
{
public:

    // The fields of a Circle
    reals a, b, r, s, g, Gx, Gy;
    int i, j;

    // constructors
    Circle();
    Circle(reals aa, reals bb, reals rr);

    // routines
    void print(void);

    // no destructor we didn't allocate memory by hand.
};


#endif // CIRCLE_H
