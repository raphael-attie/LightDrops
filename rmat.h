#ifndef RMAT_H
#define RMAT_H

#include "winsockwrapper.h"
#include <QString>
#include <QTreeWidgetItem>
#include <QFileInfo>
//opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//using namespace cv;


enum class instruments {USET, MAG, DSLR, generic};

class RMat
{
public:
    RMat();
    RMat(RMat const& rMat);
    RMat(cv::Mat mat);
    RMat(cv::Mat mat, bool bayer);
    RMat(cv::Mat mat, bool bayer, instruments instrument);
    ~RMat();

    cv::Mat matImage;

    // Methods for getting some statistics
    void computeHist(int nBins, float minRange, float maxRange);
    float calcMedian();
    float calcThreshold(float cutOff, double histWidth, float minRange);

    void calcStats();
    void calcMinMax();


    // getters
    bool isBayer() const;
    float getBscale() const;
    int getBzero() const;
    double getDataMin() const;
    double getDataMax() const;
    float getExpTime() const;
    float getWbRed() const;
    float getWbGreen() const;
    float getWbBlue() const;
    instruments getInstrument() const;
    QString getImageTitle() const;
    cv::Mat getMatHist() const;
    float getMean() const;
    float getStdDev() const;
    float getMedian() const;
    float getIntensityLow() const;
    float getIntensityHigh() const;
    double getHistWidth() const;
    uint getNPixels() const;
    float getMinHistRange() const;
    float getMaxHistRange() const;
    float getDataRange() const;
    float getNormalizeRange() const;
    QTreeWidgetItem* getItem() const;
    QFileInfo getFileInfo() const;
    QString getDate_obs() const;
    QString getTime_obs() const;
    QString getDate_time() const;



    // setters
    void setBayer(bool bayer);
    void setBscale(float bscale);
    void setBzero(int bzero);
    void setDataMin(float dataMin);
    void setDataMax(float dataMax);
    void setExpTime(float expTime);
    void setWbRed(float wbRed);
    void setWbGreen(float wbGreen);
    void setWbBlue(float wbBlue);
    void setImageTitle(QString title);
    void setInstrument(instruments instrument);
    void setItem(QTreeWidgetItem* item);
    void setFileInfo(QFileInfo fileInfo);
    void setDate_obs(QString date_obs);
    void setTime_obs(QString time_obs);
    void setDate_time(QString date_time);



private:

   bool bayer;
   float bscale;
   int bzero;
   double dataMin;
   double dataMax;
   float expTime;
   float wbRed;
   float wbGreen;
   float wbBlue;

   instruments instrument;
   QString imageTitle;
   QFileInfo fileInfo;
   QString date_obs;
   QString time_obs;
   QString date_time;

   cv::Mat matImageGray;

   // Image statistics
   cv::Mat matHist;
   float mean;
   float stdDev;
   float median;
   float intensityLow, intensityHigh;
   double histWidth;
   float minHistRange;
   float maxHistRange;
   uint nPixels;
   float dataRange;

   // Normalization range
   float normalizeRange;

   // for QTreeWidget
   QTreeWidgetItem* item;


};

#endif // RMAT_H
