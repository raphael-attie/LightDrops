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
        sum += pow(sqrt(dx*dx+dy*dy) - circle.r, 2.0);
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

int CircleFitByLevenbergMarquardtReduced (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle)
/*                                        <------------------ Input ------------------->  <-- Output -->

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

       Algorithm:  Levenberg-Marquardt running over the "reduced" parameter space (a,b)
                   (the radius r is always set to its optimal value)

       See a detailed description in Section 4.6 of the book by Nikolai Chernov:
       "Circular and linear regression: Fitting circles and lines by least squares"
       Chapman & Hall/CRC, Monographs on Statistics and Applied Probability, volume 117, 2010.

        Nikolai Chernov,  February 2014
*/
{
    int code,i,pivot,iter,inner,IterMAX=99;

    reals factorUp=10.,factorDown=0.04,lambda,ParLimit=1.e+6;
    reals dx,dy,ri,u,v;
    reals Mu,Mv,Muu,Mvv,Muv,Mr,A11,A12,A22,F1,F2,dX,dY;
    reals epsilon=3.e-8;
    reals G11,G12,G22,D1,D2;

    Circle Old,New;

    data.means();   // Compute x- and y-means (via a function in class "data")

//       starting with the given initial circle (initial guess)

    New = circleIni;

//     compute the root-mean-square error via function SigmaReduced; see Utilities.cpp

    New.s = SigmaReduced(data,New);

//     initializing lambda, iteration counters, and the exit code

    lambda = LambdaIni;
    iter = inner = code = 0;

NextIteration:

    Old = New;
    if (++iter > IterMAX) {code = 1;  goto enough;}

//     computing moments

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

//    computing matrices

    F1 = Old.a + Mu*Mr - data.meanX;
    F2 = Old.b + Mv*Mr - data.meanY;

try_again:

    A11 = (Muu - Mu*Mu) + lambda;
    A22 = (Mvv - Mv*Mv) + lambda;
    A12 = Muv - Mu*Mv;

    if (A11<epsilon || A22<epsilon || A11*A22-A12*A12<epsilon)
    {
        lambda *= factorUp;
        goto try_again;
    }

//      Cholesky decomposition with pivoting

    pivot=0;
    if (A11 < A22)
    {
        swap(A11,A22);
        swap(F1,F2);
        pivot=1;
    }

    G11 = sqrt(A11);
    G12 = A12/G11;
    G22 = sqrt(A22 - G12*G12);

    D1 = F1/G11;
    D2 = (F2 - G12*D1)/G22;

//    updating the parameters

    dY = D2/G22;
    dX = (D1 - G12*dY)/G11;

    if ((abs(dX)+abs(dY))/(One+abs(Old.a)+abs(Old.b)) < epsilon) goto enough;

    if (pivot != 0) swap(dX,dY);

    New.a = Old.a - dX;
    New.b = Old.b - dY;

    if (abs(New.a)>ParLimit || abs(New.b)>ParLimit) {code = 3; goto enough;}

//    compute the root-mean-square error via function SigmaReduced; see Utilities.cpp

    New.s = SigmaReduced(data,New);

//      check if improvement is gained

    if (New.s < Old.s)    //    yes, improvement
    {
        lambda *= factorDown;
        goto NextIteration;
    }
    else                  //  no improvement
    {
        if (++inner > IterMAX) {code = 2;  goto enough;}
        lambda *= factorUp;
        goto try_again;
    }

//       exit

enough:

//    compute the optimal radius via a function in Utilities.cpp

    Old.r = OptimalRadius(data,Old);
    Old.i = iter;
    Old.j = inner;

    circle = Old;

    return code;
}

int CircleFitByChernovLesort (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle)
/*                            <------------------ Input ------------------->  <-- Output -->

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

       Algorithm by Nikolai Chernov and Claire Lesort

       See a detailed description in the journal paper:

       N. Chernov and C. Lesort, "Least squares fitting of circles"
          in J. Math. Imag. Vision, volume 23, (2005), pages 239-251.

       the algorithm is designed to converge from any initial guess,
       but it is complicated and generally very slow

        Nikolai Chernov,  February 2014
*/
{
    int code,i,iter,inner,IterMAX=99;

    reals factorUp=10.,factorDown=0.04,lambda;
    reals Aold,Fold,Told,Anew,Fnew,Tnew,DD,H,aabb;
    reals Xi,Yi,Zi,Ui,Vi,Gi,CT,ST,D,ADF,SQ,DEN,FACT,DGDAi,DGDFi,DGDTi;
    reals H11,H12,H13,H22,H23,H33,F1,F2,F3,dA,dF,dT;
    reals epsilon=3.e-8;
    reals G11,G22,G33,G12,G13,G23,D1,D2,D3;
    reals Xshift=0.,Yshift=0.,dX=One,dY=0.,aTemp,bTemp,rTemp;

    Circle Old,New;

//       starting with the given initial circle (initial guess)

    New = circleIni;

//       compute the root-mean-square error via function Sigma; see Utilities.cpp

    New.s = Sigma(data,New);

    Anew = One/Two/New.r;
    aabb = New.a*New.a + New.b*New.b;
    Fnew = (aabb - New.r*New.r)*Anew;
    Tnew = acos(-New.a/sqrt(aabb));
    if (New.b > 0.) Tnew = Two*Pi - Tnew;

    if (One+Four*Anew*Fnew < epsilon)
    {
        Xshift += dX;
        Yshift += dY;

        New.a += dX;
        New.b += dY;
        aabb = New.a*New.a + New.b*New.b;
        Fnew = (aabb - New.r*New.r)*Anew;
        Tnew = acos(-New.a/sqrt(aabb));
        if (New.b > 0.) Tnew = Two*Pi - Tnew;
    }

//       initializing lambda, iteration counters, and the exit code

    lambda = LambdaIni;
    iter = inner = code = 0;

NextIteration:

    Aold = Anew;
    Fold = Fnew;
    Told = Tnew;
    Old = New;

    if (++iter > IterMAX) {code = 1;  goto enough;}

//       computing moments

shiftXY:

    DD = One + Four*Aold*Fold;
    D = sqrt(DD);
    CT = cos(Told);
    ST = sin(Told);

    H11=H12=H13=H22=H23=H33=F1=F2=F3=0.;

    for (i=0; i<data.n; i++)
    {
        Xi = data.X[i] + Xshift;
        Yi = data.Y[i] + Yshift;
        Zi = Xi*Xi + Yi*Yi;
        Ui = Xi*CT + Yi*ST;
        Vi =-Xi*ST + Yi*CT;

        ADF = Aold*Zi + D*Ui + Fold;
        SQ = sqrt(Four*Aold*ADF + One);
        DEN = SQ + One;
        Gi = Two*ADF/DEN;
        FACT = Two/DEN*(One - Aold*Gi/SQ);
        DGDAi = FACT*(Zi + Two*Fold*Ui/D) - Gi*Gi/SQ;
        DGDFi = FACT*(Two*Aold*Ui/D + One);
        DGDTi = FACT*D*Vi;

        H11 += DGDAi*DGDAi;
        H12 += DGDAi*DGDFi;
        H13 += DGDAi*DGDTi;
        H22 += DGDFi*DGDFi;
        H23 += DGDFi*DGDTi;
        H33 += DGDTi*DGDTi;

        F1 += Gi*DGDAi;
        F2 += Gi*DGDFi;
        F3 += Gi*DGDTi;
    }
    Old.g = New.g = sqrt(F1*F1 + F2*F2 + F3*F3);

try_again:

//        Cholesky decomposition

    G11 = sqrt(H11 + lambda);
    G12 = H12/G11;
    G13 = H13/G11;
    G22 = sqrt(H22 + lambda - G12*G12);
    G23 = (H23 - G12*G13)/G22;
    G33 = sqrt(H33 + lambda - G13*G13 - G23*G23);

    D1 = F1/G11;
    D2 = (F2 - G12*D1)/G22;
    D3 = (F3 - G13*D1 - G23*D2)/G33;

    dT = D3/G33;
    dF = (D2 - G23*dT)/G22;
    dA = (D1 - G12*dF - G13*dT)/G11;

//       updating the parameters

    Anew = Aold - dA;
    Fnew = Fold - dF;
    Tnew = Told - dT;

    if (One+Four*Anew*Fnew < epsilon)
    {
        Xshift += dX;
        Yshift += dY;

        H = sqrt(One+Four*Aold*Fold);
        aTemp = -H*cos(Told)/(Aold+Aold) + dX;
        bTemp = -H*sin(Told)/(Aold+Aold) + dY;
        rTemp = One/abs(Aold+Aold);

        Aold = One/(rTemp + rTemp);
        aabb = aTemp*aTemp + bTemp*bTemp;
        Fold = (aabb - rTemp*rTemp)*Aold;
        Told = acos(-aTemp/sqrt(aabb));
        if (bTemp > 0.) Told = Two*Pi - Told;

        lambda *= factorUp;
        inner++;
        goto shiftXY;
    }

    H = sqrt(One+Four*Anew*Fnew);
    New.a = -H*cos(Tnew)/(Anew+Anew) - Xshift;
    New.b = -H*sin(Tnew)/(Anew+Anew) - Yshift;
    New.r = One/abs(Anew+Anew);
    New.s = Sigma(data,New);

    if ((abs(New.a-Old.a) + abs(New.b-Old.b) + abs(New.r-Old.r))/(One + Old.r) < epsilon) goto enough;

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

    Old.i = iter;
    Old.j = inner;

    circle = Old;

    return code;
}

//****************** Perturb *********************************

Circle Perturb (Circle& New, Circle& Old, reals range)
{
    Circle Perturbed;

    if (range==0.) return New;

    Perturbed.a = New.a + (New.a - Old.a)*(range*rand()/RAND_MAX - range/Two);
    Perturbed.b = New.b + (New.b - Old.b)*(range*rand()/RAND_MAX - range/Two);
    Perturbed.r = New.r + (New.r - Old.r)*(range*rand()/RAND_MAX - range/Two);

    return Perturbed;
}

Circle CircleFitByPratt (Data& data)
/*
      Circle fit to a given set of data points (in 2D)

      This is an algebraic fit, due to Pratt, based on the journal article

      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
      Computer Graphics, Vol. 21, pages 145-152 (1987)

      Input:  data     - the class of data (contains the given points):

          data.n   - the number of data points
          data.X[] - the array of X-coordinates
          data.Y[] - the array of Y-coordinates

     Output:
               circle - parameters of the fitting circle:

           circle.a - the X-coordinate of the center of the fitting circle
           circle.b - the Y-coordinate of the center of the fitting circle
           circle.r - the radius of the fitting circle
           circle.s - the root mean square error (the estimate of sigma)
           circle.j - the total number of iterations

     The method is based on the minimization of the function

                 F = sum [(x-a)^2 + (y-b)^2 - R^2]^2 / R^2

     This method is more balanced than the simple Kasa fit.

     It works well whether data points are sampled along an entire circle or
     along a small arc.

     It still has a small bias and its statistical accuracy is slightly
     lower than that of the geometric fit (minimizing geometric distances).

     It provides a good initial guess for a subsequent geometric fit.

       Nikolai Chernov  (September 2012)

*/
{
    int i,iter,IterMAX=99;

    reals Xi,Yi,Zi;
    reals Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
    reals A0,A1,A2,A22;
    reals Dy,xnew,x,ynew,y;
    reals DET,Xcenter,Ycenter;

    Circle circle;

    data.means();   // Compute x- and y- sample means (via a function in the class "data")

//     computing moments

    Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

    for (i=0; i<data.n; i++)
    {
        Xi = data.X[i] - data.meanX;   //  centered x-coordinates
        Yi = data.Y[i] - data.meanY;   //  centered y-coordinates
        Zi = Xi*Xi + Yi*Yi;

        Mxy += Xi*Yi;
        Mxx += Xi*Xi;
        Myy += Yi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
        Mzz += Zi*Zi;
    }
    Mxx /= data.n;
    Myy /= data.n;
    Mxy /= data.n;
    Mxz /= data.n;
    Myz /= data.n;
    Mzz /= data.n;

//    computing coefficients of the characteristic polynomial

    Mz = Mxx + Myy;
    Cov_xy = Mxx*Myy - Mxy*Mxy;
    Var_z = Mzz - Mz*Mz;

    A2 = Four*Cov_xy - Three*Mz*Mz - Mzz;
    A1 = Var_z*Mz + Four*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    A22 = A2 + A2;

//    finding the root of the characteristic polynomial
//    using Newton's method starting at x=0
//     (it is guaranteed to converge to the right root)

    for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
        Dy = A1 + x*(A22 + 16.*x*x);
        xnew = x - y/Dy;
        if ((xnew == x)||(!isfinite(xnew))) break;
        ynew = A0 + xnew*(A1 + xnew*(A2 + Four*xnew*xnew));
        if (abs(ynew)>=abs(y))  break;
        x = xnew;  y = ynew;
    }

//    computing paramters of the fitting circle

    DET = x*x - x*Mz + Cov_xy;
    Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/Two;
    Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/Two;

//       assembling the output

    circle.a = Xcenter + data.meanX;
    circle.b = Ycenter + data.meanY;
    circle.r = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz + x + x);
    circle.s = Sigma(data,circle);
    circle.i = 0;
    circle.j = iter;  //  return the number of iterations, too

    return circle;
}

Circle CircleFitByHyper (Data& data)
/*
      Circle fit to a given set of data points (in 2D)

      This is an algebraic fit based on the journal article

      A. Al-Sharadqah and N. Chernov, "Error analysis for circle fitting algorithms",
      Electronic Journal of Statistics, Vol. 3, pages 886-911, (2009)

      It is an algebraic circle fit with "hyperaccuracy" (with zero essential bias).
      The term "hyperaccuracy" first appeared in papers by Kenichi Kanatani around 2006

      Input:  data     - the class of data (contains the given points):

          data.n   - the number of data points
          data.X[] - the array of X-coordinates
          data.Y[] - the array of Y-coordinates

     Output:
               circle - parameters of the fitting circle:

           circle.a - the X-coordinate of the center of the fitting circle
           circle.b - the Y-coordinate of the center of the fitting circle
           circle.r - the radius of the fitting circle
           circle.s - the root mean square error (the estimate of sigma)
           circle.j - the total number of iterations

     This method combines the Pratt and Taubin fits to eliminate the essential bias.

     It works well whether data points are sampled along an entire circle or
     along a small arc.

     Its statistical accuracy is theoretically higher than that of the Pratt fit
     and Taubin fit, but practically they all return almost identical circles
     (unlike the Kasa fit that may be grossly inaccurate).

     It provides a very good initial guess for a subsequent geometric fit.

       Nikolai Chernov  (September 2012)

*/
{
    int i,iter,IterMAX=99;

    reals Xi,Yi,Zi;
    reals Mz,Mxy,Mxx,Myy,Mxz,Myz,Mzz,Cov_xy,Var_z;
    reals A0,A1,A2,A22;
    reals Dy,xnew,x,ynew,y;
    reals DET,Xcenter,Ycenter;

    Circle circle;

    data.means();   // Compute x- and y- sample means (via a function in the class "data")

//     computing moments

    Mxx=Myy=Mxy=Mxz=Myz=Mzz=0.;

    for (i=0; i<data.n; i++)
    {
        Xi = data.X[i] - data.meanX;   //  centered x-coordinates
        Yi = data.Y[i] - data.meanY;   //  centered y-coordinates
        Zi = Xi*Xi + Yi*Yi;

        Mxy += Xi*Yi;
        Mxx += Xi*Xi;
        Myy += Yi*Yi;
        Mxz += Xi*Zi;
        Myz += Yi*Zi;
        Mzz += Zi*Zi;
    }
    Mxx /= data.n;
    Myy /= data.n;
    Mxy /= data.n;
    Mxz /= data.n;
    Myz /= data.n;
    Mzz /= data.n;

//    computing the coefficients of the characteristic polynomial

    Mz = Mxx + Myy;
    Cov_xy = Mxx*Myy - Mxy*Mxy;
    Var_z = Mzz - Mz*Mz;

    A2 = Four*Cov_xy - Three*Mz*Mz - Mzz;
    A1 = Var_z*Mz + Four*Cov_xy*Mz - Mxz*Mxz - Myz*Myz;
    A0 = Mxz*(Mxz*Myy - Myz*Mxy) + Myz*(Myz*Mxx - Mxz*Mxy) - Var_z*Cov_xy;
    A22 = A2 + A2;

//    finding the root of the characteristic polynomial
//    using Newton's method starting at x=0
//     (it is guaranteed to converge to the right root)

    for (x=0.,y=A0,iter=0; iter<IterMAX; iter++)  // usually, 4-6 iterations are enough
    {
        Dy = A1 + x*(A22 + 16.*x*x);
        xnew = x - y/Dy;
        if ((xnew == x)||(!isfinite(xnew))) break;
        ynew = A0 + xnew*(A1 + xnew*(A2 + Four*xnew*xnew));
        if (abs(ynew)>=abs(y))  break;
        x = xnew;  y = ynew;
    }

//    computing paramters of the fitting circle

    DET = x*x - x*Mz + Cov_xy;
    Xcenter = (Mxz*(Myy - x) - Myz*Mxy)/DET/Two;
    Ycenter = (Myz*(Mxx - x) - Mxz*Mxy)/DET/Two;

//       assembling the output

    circle.a = Xcenter + data.meanX;
    circle.b = Ycenter + data.meanY;
    circle.r = sqrt(Xcenter*Xcenter + Ycenter*Ycenter + Mz - x - x);
    circle.s = Sigma(data,circle);
    circle.i = 0;
    circle.j = iter;  //  return the number of iterations, too

    return circle;
}

