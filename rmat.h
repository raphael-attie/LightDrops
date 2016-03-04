#ifndef RMAT_H
#define RMAT_H

#include <QString>
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
    RMat(cv::Mat &matImage);
    RMat(cv::Mat &matImage, bool bayer);
    RMat(cv::Mat &matImage, bool bayer, instruments instrument);
    ~RMat();

    cv::Mat matImage;

    bool empty();

    // Methods for getting some statistics
    void computeHist(int nBins, float minRange, float maxRange);
    float calcMedian();
    float calcThreshold(float cutOff, double histWidth, float minRange);

    void calcStats();

    // getters
    bool isBayer() const;
    int getBscale() const;
    int getBzero() const;
    double getDataMin() const;
    double getDataMax() const;
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


    // setters
    void setBayer(bool bayer);
    void setBscale(int bscale);
    void setBzero(int bzero);
    void setDataMin(float dataMin);
    void setDataMax(float dataMax);
    void setWbRed(float wbRed);
    void setWbGreen(float wbGreen);
    void setWbBlue(float wbBlue);
    void setImageTitle(QString title);
    void setInstrument(instruments instrument);

private:

   bool bayer;
   int bscale;
   int bzero;
   double dataMin;
   double dataMax;
   float wbRed;
   float wbGreen;
   float wbBlue;

   instruments instrument;
   QString imageTitle;

   cv::Mat matImageGray;

   // Image statistics
   cv::Mat matHist;
   float mean;
   float stdDev;
   float median;
   float intensityLow, intensityHigh;
   double histWidth;
   uint nPixels;


};

#endif // RMAT_H
