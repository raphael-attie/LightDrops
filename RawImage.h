#ifndef RAWWIDGET_H
#define RAWWIDGET_H

#include "winsockwrapper.h"
#include <iostream>
#include <QString>
#include <QVector>
#include <QColor>

//opencv
#include <opencv2/core.hpp>
#include <opencv2/core/ocl.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

// libraw
#include <libraw/libraw.h>

// Exiv2
//#include <exiv2/exiv2.hpp>
#include "exiv2/exiv2.hpp"

class RawImage
{
public:
    RawImage();
    RawImage(QString filePath);
    ~RawImage();

    void extractExif();

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

    QVector<QString> getKeyNames() const;
    QVector<QString> getKeyValues() const;
    QVector<QString> getKeyComments() const;

    void dispatchMetaDatum(const char *keyName, QString keyValue, const char *comment = "");


private:

    QString filePath_;
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

    QVector<QString> keyNames;
    QVector<QString> keyValues;
    QVector<QString> keyComments;

    Exiv2::Image::AutoPtr image;




};

#endif // RAWWIDGET_H
