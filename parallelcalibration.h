#ifndef PARALLELCALIBRATION_H
#define PARALLELCALIBRATION_H

#include "winsockwrapper.h"
#include <QtCore>

//opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/video.hpp>
#include <opencv2/core/ocl.hpp>

#include "imagemanager.h"
#include "rmat.h"


class ParallelCalibration : public cv::ParallelLoopBody
{

public:
    ParallelCalibration(QList<QUrl> lightFiles, RMat* masterDark, RMat* masterFlat, QList<RMat*> rMatList);

    virtual void operator()(const cv::Range& range) const;

private:

    QList<QUrl> files;
    RMat dark;
    RMat flat;
    QList<RMat*> resultList;
};

#endif // PARALLELCALIBRATION_H
