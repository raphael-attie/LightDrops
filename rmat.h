#ifndef RMAT_H
#define RMAT_H

#include "winsockwrapper.h"
#include <QUrl>
#include <QString>
#include <QTreeWidgetItem>
#include <QFileInfo>
//opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

//using namespace cv;


enum class instruments {generic, MAG, DSLR, USET, TIFF};

class RMat
{
public:
    RMat();
    RMat(cv::Mat mat);
    RMat(RMat const& rMat);
    RMat(cv::Mat mat, RMat* rMat);
    RMat(cv::Mat mat, bool bayer);
    RMat(cv::Mat mat, bool bayer, instruments instrument);
    RMat(cv::Mat mat, bool bayer, instruments instrument, float XPOSURE, float TEMP);
    ~RMat();

    cv::Mat matImage;
    cv::Mat matImageGray;
    cv::Mat matImageRGB;

    void initialize();
    void prepImages();
    // Methods for getting some statistics
    void computeHist(int nBins, float minRange, float maxRange);
    float calcMedian(double histWidth, float minRange);
    float calcThreshold(float cutOff, double histWidth, float minRange);

    void calcStats();
    void calcMinMax();

    // Extract channels
    cv::Mat extractChannel(unsigned int channel);

    // Used to know when to flip image upside down
    // Yeah I know, it's disgusting... but i'm tired!
    bool flipUD;


    // getters
    bool isBayer() const;
    float getBscale() const;
    int getBzero() const;
    double getDataMin() const;
    double getDataMax() const;
    float getExpTime() const;
    float getXPOSURE() const;
    float getTEMP() const;
    float getSOLAR_R() const;
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
    QUrl getUrl() const;



    // setters
    void setBayer(bool bayer);
    void setBscale(float bscale);
    void setBzero(int bzero);
    void setDataMin(float dataMin);
    void setDataMax(float dataMax);
    void setExpTime(float expTime);
    void setSOLAR_R(float SOLAR_R);
    void setXPOSURE(float XPOSURE);
    void setTEMP(float temperature);
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
    void setUrl(QUrl url);

private:

   bool bayer;
   float bscale;
   int bzero;
   double dataMin;
   double dataMax;
   float expTime;
   float XPOSURE;
   float TEMP;
   float SOLAR_R;
   float wbRed;
   float wbGreen;
   float wbBlue;

   instruments instrument;
   QString imageTitle;
   QFileInfo fileInfo;
   QUrl url;
   QString date_obs;
   QString time_obs;
   QString date_time;

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
