#ifndef PARALLELCALIBRATION_H
#define PARALLELCALIBRATION_H


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
private:

    QString exportDir;
    QList<QUrl> files;
    RMat* dark;
    RMat* flat;
    QList<RMat*> resultList;

public:
    ParallelCalibration(QString dir, QList<QUrl> lightFiles, RMat* masterDark, RMat* masterFlat, QList<RMat*> rMatList);

    virtual void operator()(const cv::Range& range) const;

};

#endif // PARALLELCALIBRATION_H
