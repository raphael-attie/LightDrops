#include "typedefs.h"
#include "data.h"
#include "circle.h"

//****************** Sigma ************************************
//
//   estimate of Sigma = square root of RSS divided by N
//   gives the root-mean-square error of the geometric circle fit

reals Sigma (Data& data, Circle& circle)
{
    reals sum=0.,dx,dy;

    for (int i=0; i<data.n; i++)
    {
        dx = data.X[i] - circle.a;
        dy = data.Y[i] - circle.b;
        sum += pow(sqrt(dx*dx+dy*dy) - circle.r, 2);
    }
    return sqrt(sum/data.n);
}

//****************** SigmaReduced ************************************
//
//   estimate of Sigma = square root of RSS divided by N
//   gives the root-mean-square error of the geometric circle fit
//
//   uses only the center of the circle (a,b), not the radius
//   the function computes the optimal radius here

reals SigmaReduced (Data& data, Circle& circle)
{
    int i,n=data.n;
    reals sum=0.,dx,dy,r,D[n];

    for (i=0; i<n; i++)
    {
        dx = data.X[i] - circle.a;
        dy = data.Y[i] - circle.b;
        D[i] = sqrt(dx*dx+dy*dy);
        sum += D[i];
    }
    r = sum/n;

    for (sum=0., i=0; i<n; i++)  sum += pow(D[i] - r, 2.0);

    return sqrt(sum/n);
}

//****************** SigmaReducedNearLinearCase ****************
//
//   estimate of Sigma = square root of RSS divided by N

reals SigmaReducedNearLinearCase (Data& data, Circle& circle)
{
    int i,n=data.n;
    reals a0,b0,del,s,c,x,y,z,p,t,g,W,Z;

    a0 = circle.a-data.meanX;  b0 = circle.b-data.meanY;
    del = One/sqrt(a0*a0 + b0*b0);
    s = b0*del;  c = a0*del;

    for (W=Z=0.,i=0; i<n; i++)
    {
        x = data.X[i] - data.meanX;
        y = data.Y[i] - data.meanY;
        z = x*x + y*y;
        p = x*c + y*s;
        t = del*z - Two*p;
        g = t/(One+sqrt(One+del*t));
        W += (z+p*g)/(Two+del*g);
        Z += z;
    }
    W /= n;
    Z /= n;

    return sqrt(Z-W*(Two+del*del*W));
}

//****************** SigmaReducedForCenteredScaled ****************
//
//   estimate of Sigma = square root of RSS divided by N

reals SigmaReducedForCenteredScaled (Data& data, Circle& circle)
{
    int i,n=data.n;
    reals sum=0.,dx,dy,r;

    for (i=0; i<n; i++)
    {
        dx = data.X[i] - circle.a;
        dy = data.Y[i] - circle.b;
        sum += sqrt(dx*dx+dy*dy);
    }
    r = sum/n;

    return sqrt(pow(circle.a, 2.0)+pow(circle.b, 2.0)-r*r+Two);
}

//****************** OptimalRadius ******************************
//
//     compute the optimal radius of a circle, given its center (a,b)

reals OptimalRadius (Data& data, Circle& circle)
{
    reals Mr=0.,dx,dy;

    for (int i=0; i<data.n; i++)
    {
        dx = data.X[i] - circle.a;
        dy = data.Y[i] - circle.b;
        Mr += sqrt(dx*dx + dy*dy);
    }
    return Mr/data.n;
}

int CircleFitByLevenbergMarquardtFull (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle)
/*                                     <------------------ Input ------------------->  <-- Output -->

       Geometric circle fit to a given set of data points (in 2D)

       Input:  data     - the class of data (contains the given points):

           data.n   - the number of data points
           data.X[] - the array of X-coordinates
           data.Y[] - the array of Y-coordinates

               circleIni - parameters of the initial circle ("initial guess")

           circleIni.a - the X-coordinate of the center of the initial circle
           circleIni.b - the Y-coordinate of the center of the initial circle
           circleIni.r - the radius of the initial circle

           LambdaIni - the initial value of the control parameter "lambda"
                       for the Levenberg-Marquardt procedure
                       (common choice is a small positive number, e.g. 0.001)

       Output:
           integer function value is a code:
                      0:  normal termination, the best fitting circle is
                          successfully found
                      1:  the number of outer iterations exceeds the limit (99)
                          (indicator of a possible divergence)
                      2:  the number of inner iterations exceeds the limit (99)
                          (another indicator of a possible divergence)
                      3:  the coordinates of the center are too large
                          (a strong indicator of divergence)

           circle - parameters of the fitting circle ("best fit")

           circle.a - the X-coordinate of the center of the fitting circle
           circle.b - the Y-coordinate of the center of the fitting circle
           circle.r - the radius of the fitting circle
           circle.s - the root mean square error (the estimate of sigma)
           circle.i - the total number of outer iterations (updating the parameters)
           circle.j - the total number of inner iterations (adjusting lambda)

       Algorithm:  Levenberg-Marquardt running over the full parameter space (a,b,r)

       See a detailed description in Section 4.5 of the book by Nikolai Chernov:
       "Circular and linear regression: Fitting circles and lines by least squares"
       Chapman & Hall/CRC, Monographs on Statistics and Applied Probability, volume 117, 2010.

        Nikolai Chernov,  February 2014
*/
{
    int code,i,iter,inner,IterMAX=99;

    reals factorUp=10.,factorDown=0.04,lambda,ParLimit=1.e+6;
    reals dx,dy,ri,u,v;
    reals Mu,Mv,Muu,Mvv,Muv,Mr,UUl,VVl,Nl,F1,F2,F3,dX,dY,dR;
    reals epsilon=3.e-8;
    reals G11,G22,G33,G12,G13,G23,D1,D2,D3;

    Circle Old,New;

//       starting with the given initial circle (initial guess)

    New = circleIni;

//       compute the root-mean-square error via function Sigma; see Utilities.cpp

    New.s = Sigma(data,New);

//       initializing lambda, iteration counters, and the exit code

    lambda = LambdaIni;
    iter = inner = code = 0;

NextIteration:

    Old = New;
    if (++iter > IterMAX) {code = 1;  goto enough;}

//       computing moments

    Mu=Mv=Muu=Mvv=Muv=Mr=0.;

    for (i=0; i<data.n; i++)
    {
        dx = data.X[i] - Old.a;
        dy = data.Y[i] - Old.b;
        ri = sqrt(dx*dx + dy*dy);
        u = dx/ri;
        v = dy/ri;
        Mu += u;
        Mv += v;
        Muu += u*u;
        Mvv += v*v;
        Muv += u*v;
        Mr += ri;
    }
    Mu /= data.n;
    Mv /= data.n;
    Muu /= data.n;
    Mvv /= data.n;
    Muv /= data.n;
    Mr /= data.n;

//       computing matrices

    F1 = Old.a + Old.r*Mu - data.meanX;
    F2 = Old.b + Old.r*Mv - data.meanY;
    F3 = Old.r - Mr;

    Old.g = New.g = sqrt(F1*F1 + F2*F2 + F3*F3);

try_again:

    UUl = Muu + lambda;
    VVl = Mvv + lambda;
    Nl = One + lambda;

//         Cholesly decomposition

    G11 = sqrt(UUl);
    G12 = Muv/G11;
    G13 = Mu/G11;
    G22 = sqrt(VVl - G12*G12);
    G23 = (Mv - G12*G13)/G22;
    G33 = sqrt(Nl - G13*G13 - G23*G23);

    D1 = F1/G11;
    D2 = (F2 - G12*D1)/G22;
    D3 = (F3 - G13*D1 - G23*D2)/G33;

    dR = D3/G33;
    dY = (D2 - G23*dR)/G22;
    dX = (D1 - G12*dY - G13*dR)/G11;

    if ((abs(dR)+abs(dX)+abs(dY))/(One+Old.r) < epsilon) goto enough;

//       updating the parameters

    New.a = Old.a - dX;
    New.b = Old.b - dY;

    if (abs(New.a)>ParLimit || abs(New.b)>ParLimit) {code = 3; goto enough;}

    New.r = Old.r - dR;

    if (New.r <= 0.)
    {
        lambda *= factorUp;
        if (++inner > IterMAX) {code = 2;  goto enough;}
        goto try_again;
    }

//       compute the root-mean-square error via function Sigma; see Utilities.cpp

    New.s = Sigma(data,New);

//       check if improvement is gained

    if (New.s < Old.s)    //   yes, improvement
    {
        lambda *= factorDown;
        goto NextIteration;
    }
    else                       //   no improvement
    {
        if (++inner > IterMAX) {code = 2;  goto enough;}
        lambda *= factorUp;
        goto try_again;
    }

    //       exit

enough:

    Old.i = iter;    // total number of outer iterations (updating the parameters)
    Old.j = inner;   // total number of inner iterations (adjusting lambda)

    circle = Old;

    return code;
}

