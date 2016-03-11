#ifndef IMAGEMANAGER_H
#define IMAGEMANAGER_H

#include <QtCore>
#include <QTableWidget>

//opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "MyFitsImage.h"
#include "RawImage.h"
#include "rmat.h"

class ImageManager
{
public:
    ImageManager(QString filePathQStre);
    ~ImageManager();

    RMat rMatImage;

    // getters
    MyFitsImage* getNewFitsImage();
    QString getFileName();
    QTableWidget *getTableWidget() const;

    float getWbRed() const;
    float getWbGreen() const;
    float getWbBlue() const;

private:

    void loadFits();
    void loadRaw();
    void createTableWidget();

    QList<QString> fitsExtList;
    QList<QString> rawExtList;
    QString filePathQStr;
    QString fileName;
    QString fileExt;
    QTableWidget* tableWidget;

    MyFitsImage *newFitsImage;
    QVector<MyFitsImage*> fitsSeries;
    RawImage *newRawImage;

    // white balance
    float wbRed;
    float wbGreen;
    float wbBlue;
    double dataMin;
    double dataMax;


};

#endif // IMAGEMANAGER_H
