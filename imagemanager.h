#ifndef IMAGEMANAGER_H
#define IMAGEMANAGER_H

#include "winsockwrapper.h"
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


    // getters
    bool getError();
    MyFitsImage* getNewFitsImage();
    QString getFileName();
    QTableWidget *getTableWidget() const;
    RMat* getRMatImage();

    float getWbRed() const;
    float getWbGreen() const;
    float getWbBlue() const;

private:

    bool error;
    void loadFits();
    void loadRaw();
    void createTableWidget();

    RMat* rMatImage;
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

};

#endif // IMAGEMANAGER_H
