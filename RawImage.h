#ifndef RAWWIDGET_H
#define RAWWIDGET_H

#include "winsockwrapper.h"
#include <iostream>
#include <QString>
#include <QColor>

//opencv
#include <opencv2/core.hpp>
#include <opencv2/core/ocl.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

// libraw
#include <libraw/libraw.h>

class RawImage
{
public:
    RawImage();
    RawImage(QString filePath);
    ~RawImage();

    LibRaw getRawProcess() const;

    cv::Mat matCFA;

    cv::Mat getImRed() const;
    cv::Mat getImGreen() const;
    cv::Mat getImBlue() const;

    qint32 getNaxis1() const;
    qint32 getNaxis2() const;

    int getNPixels() const;

    float getDataMinRed() const;
    float getDataMaxRed() const;

    float getDataMinGreen() const;
    float getDataMaxGreen() const;

    float getDataMinBlue() const;
    float getDataMaxBlue() const;

    float getWbRed() const;
    float getWbGreen() const;
    float getWbBlue() const;


private:

    LibRaw rawProcess;
    libraw_processed_image_t* imageProcessed;


    // dimensions
    int naxis1;
    int naxis2;
    int nPixels;

    float dataMinRed;
    float dataMaxRed;

    float dataMinGreen;
    float dataMaxGreen;

    float dataMinBlue;
    float dataMaxBlue;

    // White balance values from camera
    float wbRed;
    float wbGreen;
    float wbBlue;

    cv::Mat imRed;
    cv::Mat imGreen;
    cv::Mat imBlue;
    cv::Mat imGreen2;




};

#endif // RAWWIDGET_H
