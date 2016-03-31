#ifndef UTILITIES_H
#define UTILITIES_H

#include "typedefs.h"
#include "data.h"
#include "circle.h"

reals Sigma (Data& data, Circle& circle);
reals SigmaReduced (Data& data, Circle& circle);
reals SigmaReducedNearLinearCase (Data& data, Circle& circle);
reals SigmaReducedForCenteredScaled (Data& data, Circle& circle);
reals OptimalRadius (Data& data, Circle& circle);

int CircleFitByLevenbergMarquardtFull (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle);
int CircleFitByLevenbergMarquardtReduced (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle);
int CircleFitByChernovLesort (Data& data, Circle& circleIni, reals LambdaIni, Circle& circle);

Circle Perturb (Circle& New, Circle& Old, reals range);
Circle CircleFitByPratt (Data& data);
Circle CircleFitByHyper (Data& data);


#endif // UTILITIES_H

