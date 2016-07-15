#ifndef LIMB_H
#define LIMB_H

#include "utilities.h"
#include <opencv2/world.hpp>

int getMaxIncl(int *array, int length);
int getMaxIncl2(ushort *array, int length);
void FindLimb(int *data,Data *dat, int naxis1, int naxis2, int numDots);
void FindLimb2(ushort *data, Data *dat, int naxis1, int naxis2, int numDots);
void FindLimb3(cv::Mat matImage, Data *dat, int numDots);
Circle getSun(int *data, int naxis1, int naxis2, int numDots);

#endif // LIMB_H
