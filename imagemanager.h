#ifndef IMAGEMANAGER_H
#define IMAGEMANAGER_H

#include "winsockwrapper.h"
#include <QtCore>
#include <QTableWidget>

//opencv
//#include <opencv2/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/world.hpp>

#include "MyFitsImage.h"
#include "RawImage.h"
#include "rmat.h"


class ImageManager
{
public:
    ImageManager(QUrl url);
    ~ImageManager();


    // getters
    bool getError() const;
    MyFitsImage* getNewFitsImage() const;
    QUrl getUrl() const;
    QString getFileName() const;
    QTableWidget *getTableWidget() const;
    RMat* getRMatImage() const;
    QString getDate_obs() const;
    QString getTime_obs() const;

    float getWbRed() const;
    float getWbGreen() const;
    float getWbBlue() const;

private:

    bool error;
    void loadFits();
    void loadRaw();
    void createTableWidget();
    void setupFitsTableWidget();
    void setupRawTableWidget();
    void fixUset();

    RMat* rMatImage;
    QList<QString> fitsExtList;
    QList<QString> rawExtList;
    QUrl url;
    QString filePathQStr;
    QString fileName;
    QString fileExt;
    QTableWidget* tableWidget;

    MyFitsImage *newFitsImage;
    QVector<MyFitsImage*> fitsSeries;
    RawImage *newRawImage;

    // Header info
    QString date_obs;
    QString time_obs;
    uint nKeys;

    // white balance
    float wbRed;
    float wbGreen;
    float wbBlue;

};

#endif // IMAGEMANAGER_H
