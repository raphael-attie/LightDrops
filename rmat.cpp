#include "rmat.h"
#include <qdebug.h>

RMat::RMat()
{
    cv::Mat emptyMat;
    this->matImage = emptyMat;
    this->bayer = false;
    this->bscale = 1;
    this->bzero = 0;

}

RMat::RMat(const RMat &rMat)
    : bayer(rMat.bayer), bscale(rMat.bscale), bzero(rMat.bzero), expTime(rMat.expTime),
       instrument(rMat.instrument), item(NULL)
{
    rMat.matImage.copyTo(this->matImage);
    this->imageTitle = QString("Image #");


    if (matImage.channels() > 1)
    {
        cv::cvtColor(matImage, matImageGray, CV_RGB2GRAY);
    }
    else
    {
        matImageGray = matImage;
    }

    calcStats();
}



RMat::RMat(cv::Mat mat) : dataMin(0), dataMax(0), bscale(1), bzero(0), expTime(0), item(NULL)
{
    //matImage = cv::Mat(cv::Size(mat.cols, mat.rows), mat.type(), mat.data, mat.step);
    mat.copyTo(this->matImage);
    this->bayer = false;
    this->imageTitle = QString("Image #");
    this->instrument = instruments::generic;

    if (matImage.channels() > 1)
    {
        cv::cvtColor(matImage, matImageGray, CV_RGB2GRAY);
    }
    else
    {
        matImageGray = matImage;
    }

    calcStats();
}

RMat::RMat(cv::Mat mat, bool bayer) : dataMin(0), dataMax(0), bscale(1), bzero(0), expTime(0), item(NULL)
{
    //matImage = cv::Mat(cv::Size(mat.cols, mat.rows), mat.type(), mat.data, mat.step);
    mat.copyTo(this->matImage);
    this->bayer = bayer;
    this->imageTitle = QString("Image #");
    this->instrument = instruments::generic;

    if (matImage.channels() > 1)
    {
        cv::cvtColor(matImage, matImageGray, CV_RGB2GRAY);
    }
    else
    {
        matImageGray = matImage;
    }

    calcStats();
}

RMat::RMat(cv::Mat mat, bool bayer, instruments instrument) : dataMin(0), dataMax(0), bscale(1), bzero(0), expTime(0), item(NULL)
{
    //matImage = cv::Mat(cv::Size(mat.cols, mat.rows), mat.type(), mat.data, mat.step);
    mat.copyTo(this->matImage);
    this->bayer = bayer;
    this->imageTitle = QString("Image #");
    this->instrument = instrument;

    if (matImage.channels() > 1)
    {
        cv::cvtColor(matImage, matImageGray, CV_RGB2GRAY);
    }
    else
    {
        matImageGray = matImage;
    }


    calcStats();
}



RMat::~RMat()
{
    if (item != NULL)
    {
        qDebug("RMat:: deleting item from QTreeWidget");
        item->~QTreeWidgetItem();
    }
}

void RMat::computeHist(int nBins, float minRange, float maxRange)
{
    float range[2];
    range[0] = minRange;
    range[1] = maxRange;
    const float* histRange = { range };
    bool uniform = true;
    bool accumulate = false;

    /// Compute the histograms:

    if (minRange == maxRange)
    {
       matHist = cv::Mat::zeros(nBins, 1, CV_32F);
    }
    else
    {
        cv::calcHist( &matImageGray, 1, 0, cv::Mat(), matHist, 1, &nBins, &histRange, uniform, accumulate);   
    }

    float cdf = 0;
    for (int i = 0 ; i < matHist.rows ; ++i)
    {
        cdf += matHist.at<float>(i);
    }
}

void RMat::calcStats()
{

    calcMinMax();
    qDebug("RMat::calcStats():: [dataMin , dataMax] = [%f , %f]", (float) dataMin, (float) dataMax);


    /// Determine a range of the histogram
    /// These params depend on the instrument's sensor (different bit depths).

    if (instrument == instruments::USET)
    {
        maxHistRange = 4095.0f;
    }
    else if (instrument == instruments::DSLR || matImage.type() == CV_16U)
    {
        maxHistRange = 65535.0f;
    }
    else if (matImage.type() == CV_8U || matImage.type() == CV_8UC3)
    {
        maxHistRange = 255.0f;
    }
    else
    {
        maxHistRange = (float) dataMax;
    }

    minHistRange = (float) dataMin;

    cv::Scalar meanScalar;
    cv::Scalar stdDevScalar;
    cv::meanStdDev(matImageGray, meanScalar, stdDevScalar);
    mean = (float) meanScalar.val[0];
    stdDev = (float) stdDevScalar.val[0];
    nPixels = (uint) matImage.cols * matImage.rows;

    /// Histogram bin counts
    int nBins = 256;


    // Get histogram
    computeHist(nBins, minHistRange, maxHistRange + 1);
   // Calculate median
    median = calcMedian();

    histWidth = (maxHistRange - minHistRange)/(double) (nBins -1);

    /// Define percentiles for the low and high thresholds
    /// i.e, pertcentage of pixels below the lowest intensity
    /// and above the highest intensity
    float cutOffLow;
    float cutOffHigh;


    if (instrument == instruments::USET)
    {
        cutOffLow = 0;
        cutOffHigh = 99.97f;
    }
    else
    {
        cutOffLow = 5.0f;
        cutOffHigh = 99.85f;
    }

    intensityLow = calcThreshold(cutOffLow, histWidth, minHistRange);
    intensityHigh = calcThreshold(cutOffHigh, histWidth, minHistRange);

    int matType = matImage.type();
    if (instrument == instruments::USET)
    {
        dataRange = 4096.0f;
        normalizeRange = 4096.0f;
    }
    else if (matType == CV_16U || matType == CV_16UC3)
    {
        dataRange = 65536.0f;
        normalizeRange = 65536.0f;
    }
    else if (matType == CV_32F || matType == CV_32FC3)
    {
        dataRange = (float) (dataMax - dataMin) + 1.0f;
        /// For normalization / contrast stretching
        /// we assume that no CCD will go beyond 16 bits
        /// thus keeping a maximum 16-bit data range
        normalizeRange = 65536.0f;
    }
    else if (matType == CV_8U || matType == CV_8UC3)
    {
        dataRange = 256.0f;
        normalizeRange = 256.0f;
    }
    else
    {
        qDebug("RMat::calcStats():: ERROR. Unknown dataRange");
        exit(1);
    }
}

void RMat::calcMinMax()
{
    cv::minMaxLoc(matImageGray, &dataMin, &dataMax);
}



float RMat::calcMedian()
{
    float cdf = 0;
    int nBins = matHist.rows;
    float medianLimit = std::round((float)nBins/2.0f);

    float medianVal = -1;
    for (int i = 1; i < nBins && medianVal < 0.0; i++){
        cdf += matHist.at<float>(i);
        if (cdf > medianLimit) { medianVal = i;}
    }

    return medianVal;
}

float RMat::calcThreshold(float cutOff, double histWidth, float minRange)
{
    float cdf = 0.0f;
    bool checkThreshold = false;
    int nBins = matHist.rows;

    int index = 0;
    int i = 0;
    while (i < nBins && checkThreshold == false){
        cdf += matHist.at<float>(i);
        checkThreshold = (100.0f * cdf/nPixels > cutOff);
        if (checkThreshold) { index = i;}
        i++;
    }


    /// Sometimes the percentile isn't reach (for example after normalizing)
    /// This makes sure that when the maximum percentile is not reached even
    /// when the last bin has been taken into the cdf, the index is set to the last bin.
    if (!checkThreshold)
    {
        index = nBins - 1;
    }

   float intensity = minRange + (float) histWidth * (float) index;

   return intensity;
}

// Getters
bool RMat::isBayer() const
{
    return bayer;
}

float RMat::getBscale() const
{
    return bscale;
}

int RMat::getBzero() const
{
    return bzero;
}

double RMat::getDataMin() const
{
    return dataMin;
}

double RMat::getDataMax() const
{
    return dataMax;
}

float RMat::getExpTime() const
{
    return expTime;
}

float RMat::getWbRed() const
{
    return wbRed;
}

float RMat::getWbGreen() const
{
    return wbGreen;
}

float RMat::getWbBlue() const
{
    return wbBlue;
}

instruments RMat::getInstrument() const
{
    return instrument;
}

QString RMat::getImageTitle() const
{
    return imageTitle;
}

cv::Mat RMat::getMatHist() const
{
    return matHist;
}

float RMat::getMean() const
{
    return mean;
}

float RMat::getStdDev() const
{
    return stdDev;
}

float RMat::getMedian() const
{
    return median;
}

float RMat::getIntensityLow() const
{
    return intensityLow;
}

float RMat::getIntensityHigh() const
{
    return intensityHigh;
}

double RMat::getHistWidth() const
{
    return histWidth;
}

uint RMat::getNPixels() const
{
    return nPixels;
}

float RMat::getMinHistRange() const
{
    return minHistRange;
}

float RMat::getMaxHistRange() const
{
    return maxHistRange;
}

float RMat::getDataRange() const
{
    return dataRange;
}

float RMat::getNormalizeRange() const
{
    return normalizeRange;
}

QTreeWidgetItem* RMat::getItem() const
{
    return item;
}

QFileInfo RMat::getFileInfo() const
{
    return fileInfo;
}

QString RMat::getDate_obs() const
{
    return date_obs;
}

QString RMat::getTime_obs() const
{
    return time_obs;
}

QString RMat::getDate_time() const
{
    return date_time;
}


// Setters
void RMat::setBayer(bool bayer)
{
    this->bayer = bayer;
}

void RMat::setBscale(float bscale)
{
    this->bscale = bscale;
}

void RMat::setBzero(int bzero)
{
    this->bzero = bzero;
}

void RMat::setDataMin(float dataMin)
{
    this->dataMin = dataMin;
}

void RMat::setDataMax(float dataMax)
{
    this->dataMax = dataMax;
}

void RMat::setExpTime(float expTime)
{
    this->expTime = expTime;
}

void RMat::setWbRed(float wbRed)
{
    this->wbRed = wbRed;
}

void RMat::setWbGreen(float wbGreen)
{
    this->wbGreen = wbGreen;
}

void RMat::setWbBlue(float wbBlue)
{
    this->wbBlue = wbBlue;
}

void RMat::setImageTitle(QString title)
{
    this->imageTitle = title;
}

void RMat::setInstrument(instruments instrument)
{
    this->instrument = instrument;
}

void RMat::setItem(QTreeWidgetItem *item)
{
    this->item = item;
}

void RMat::setFileInfo(QFileInfo fileInfo)
{
    this->fileInfo = fileInfo;
}

void RMat::setDate_obs(QString date_obs)
{
    this->date_obs = date_obs;
}

void RMat::setTime_obs(QString time_obs)
{
    this->time_obs = time_obs;
}

void RMat::setDate_time(QString date_time)
{
    this->date_time = date_time;
}


