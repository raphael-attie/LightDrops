#include "RawImage2.h"

RawImage2::RawImage2()//rawProcess(NULL)
{
  //LibRaw rawProcess;
  matCFA.create(500, 500, CV_16U);
}

cv::Mat RawImage2::getMatCFA() const
{
    return matCFA;
}
