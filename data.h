#ifndef DATA_H
#define DATA_H

// A data has 5 fields:
//       n (of type int), the number of data points
//       X and Y (arrays of type reals), arrays of x- and y-coordinates
//       meanX and meanY (of type reals), coordinates of the centroid (x and y sample means)

#include "typedefs.h"

class Data
{
public:

    int n;
    reals *X;		//space is allocated in the constructors
    reals *Y;		//space is allocated in the constructors
    reals meanX, meanY;

    // constructors
    Data();
    Data(int N);
    Data(int N, reals *X, reals *Y);
    // destructors
    ~Data();

    // routines
    void means(void);
    void center(void);
    void scale(void);
    void print(void);

};


#endif // DATA_H
