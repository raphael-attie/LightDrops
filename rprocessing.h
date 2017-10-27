#ifndef RPROCESSING_H
#define RPROCESSING_H

#include "winsockwrapper.h"
#include <QApplication>
#include <QtCore>


#include <arrayfire.h>


#include "rmat.h"
#include "rlistimagemanager.h"
#include "rtreewidget.h"
#include "rlineedit.h"
#include "ropenglwidget.h"

#include "typedefs.h"
#include "data.h"
#include "circle.h"
#include "utilities.h"


namespace Ui {
class RMainWindow;
}

class RProcessing: public QObject
{
    Q_OBJECT

public:
    RProcessing(QObject *parent = NULL);
    ~RProcessing();

    void loadRMatLightList(QList<QUrl> urls);
    void loadRMatBiasList(QList<QUrl> urls);
    void loadRMatDarkList(QList<QUrl> urls);
    void loadRMatFlatList(QList<QUrl> urls);

    /// Calibration (bias, dark, ...)
    bool makeMasterBias();
    bool makeMasterDark();
    bool makeMasterFlat();
    void stack(QList<RMat*> rMatImageList);

    RMat *average(QList<RMat*> rMatList);
    RMat *sigmaClipAverage(QList<RMat *> rMatImageList);

    /// Fourier filtering
    double pixelDistance(double u, double v);
    cv::Mat make2DGaussian(int matSize, double sigma);
    cv::Mat fftshift(cv::Mat matFourier);
    cv::Mat makePowerSpectrumFFT(cv::Mat matImage);
    cv::Mat makeImageHPF(cv::Mat matImage, double sigma);

    /// Sharpenning
    QList<RMat*> sharpenSeries(QList<RMat*> rMatImageList, float weight1, float weight2);
    RMat* sharpenCurrentImage(RMat* rMatImage, float weight1, float weight2);

    /// Lucky imaging
    /// Block Processing

    // void blockProcessing(QList<RMat*> rMatImageList, void (*functionPtr)(af::array array3D));

    void blockProcessingLocal(QList<RMat*> rMatImageList);

    void sortBestBlocks(af::array & bestBlk, af::array & templateBlk, af::array & arfSeries, af::array & arrayBinnedSeries,
                        int blkSize, int blkSizeL, QPoint blkPos, int nFrames, int binning);

    void extractLuckySample(QList<RMat *> rMatImageList, QList<RMat *> &blockList, int x, int y, bool isCentered = false);
    void extractLuckySample2(QList<RMat *> rMatImageList, QList<RMat *> &blockList, int x, int y, bool isCentered = false);

    void blockProcessingGradient(QList<RMat*> rMatImageList, af::array &arfSeries, af::array &qualityBinnedSeries);
    void blockProcessingSobel(QList<RMat*> rMatImageList, af::array &arfSeries, af::array &qualitySeries);
    void blockProcessingLaplace(QList<RMat*> rMatImageList, af::array &arfSeries, af::array &qualitySeries);

    void blockProcessingGlobalGradients(QList<RMat*> rMatImageList);
    void blockProcessingGlobalLaplace(QList<RMat*> rMatImageList);

    void extractBestBlock(af::array & bestBlk, af::array & arfSeries, af::array & arrayBinnedSeries,
                           const int & blkSize, const int &binnedBlkSize, const int & x, const int & y,
                           const int & nFrames, const int & binning);
    void makeAlignedStack(af::array & bestBlks, const af::array & arfSeries, const af::array & qualityBinnedSeries,
                          const int nBest, const int & blkSize, const int &binnedBlkSize, int & x, int & y,
                          const int bufferSpace, const int & binning);
    void makeAlignedStack2(af::array & stackedBlks, const af::array & arfSeries, const af::array & qualityBinnedSeries,
                          const af::array & arDim, const af::array &xRange, const af::array &yRange,
                           const int nBest, const int & blkSize, const int &binnedBlkSize, const int & binning, int & x, int & y);
    void makeAlignedStackGradient(af::array & stackedBlks, const af::array & arfSeries, const af::array & qualityBinnedSeries,
                           const af::array & arDim, const af::array &xRange, const af::array &yRange,
                           int & x, int & y);
    void makeAlignedStackLaplace(af::array & stackedBlks, const af::array & arfSeries, const af::array & qualitySeries,
                           const af::array & arDim, const af::array &xRange, const af::array &yRange,
                           int & x, int & y);

    void populateResultListWithAr(QList<RMat*> rMatImageList, af::array &canvas, QString title);

    void matchTemplate2(af::array &res, af::array &searchImg, af::array &img);
    void matchTemplate3(af::array &res, af::array &A, af::array &k, af::array &arrOnes);
    void phaseCorrelate(af::array &refBlk, af::array &shiftedArray, const af::array & arDim, af::array &shifts);
    void phaseCorrelate2(af::array &stackedBlks, af::array &shifts, const af::array &arDim);


    void findMinLoc(af::array & ar, int &dx, int &dy);
    void findMinLoc(af::array & a, af::array &dx, af::array &dy);
    void findMaxLoc(af::array & a, af::array &locxy);
    void findMaxLoc2(af::array & a, af::array &locxy);

    void fetchTMatch2Shifts(af::array & a, int &dx, int &dy, af::dim4 dims);
    void fetchTMatch2Shifts(af::array & a, af::array &dx, af::array &dy, af::dim4 dims);

    void fetchTMatchShiftsGfor(af::array & ar, af::array &dxAr, af::array &dyAr);


    /// export methods
    void exportMastersToFits();
    void exportFramesToFits(QList<RMat*> rMatImageList, QDir exportDir, bool useBasename);
    void exportToFits(RMat *rMatImage, QString QStrFilename);
    void batchExportToFits(QList<QUrl> urls, QString exportDir);
    cv::Mat rescaleForExport8Bits(cv::Mat matImage, float alpha, float beta);
    void exportToTiff(RMat *rMatImage, QString QStrFilename);
    void exportToJpeg(RMat *rMatImage, QString QStrFilename);
    QString setupFileName(QFileInfo fileInfo);
    void loadMasterBias();
    void loadMasterDark();
    void loadMasterFlat();

    /// Statistics
    void showMinMax(const cv::Mat & matImage);
    cv::Mat histogram(cv::Mat matVector, int &nBins, float &width);
    float calcMedian(std::vector<float> data, float width);

    /// General matrix manipulation
    // ArrayFire matrix multiplication to be used only with af::batchFunc
    static af::array rMatMul(const af::array &lhs, const af::array &rhs);
    // Rebin
    void rebin(const af::array &a, af::array & b, int binning);
    void rebin2(const af::array &a, af::array & b2);
    void rebin4(const af::array &a, af::array & b4);
    void rebin(RMat *rMatImage, RMat *binnedRMatImage, int binning);
    // Masking & ROI
    cv::Mat circleMaskMat(cv::Mat matImage, int circleX, int circleY, int radius);
    cv::Mat circleMask(cv::Mat matImage, int circleX, int circleY, int radius);


    /// setters
    void setTreeWidget(RTreeWidget *treeWidget);
    void setCurrentROpenGLWidget(ROpenGLWidget *rOpenGLWidget);
    void setShowContours(bool status);
    void setShowLimb(bool status);
    void setUseXCorr(bool useXCorr);
    void setCvRectROI(cv::Rect cvRect);
    void setUseROI(bool status);
    void setBlurSigma(double sigma);
    void setUseHPF(bool status);
    void setHPFSigma(double sigma);
    void setSharpenLiveStatus(bool status);
    void setStackWithMean(bool status);
    void setStackWithSigmaClip(bool status);
    void setBinning(int binning);
    void setBlkSize(int blkSize);
    void setNBest(int nBest);
    void setQualityMetric(QString qualityMetric);
    void setApplyMask(bool status);
    // ROI
    void setCvRectROIList(QList<cv::Rect> cvRectList);
    void setMaskCircleX(int circleX);
    void setMaskCircleY(int circleY);
    void setMaskCircleRadius(int circleRadius);
    // TreeWidget
    void setUseUrlsFromTreeWidget(bool status);

    /// getters
    QString getExportMastersDir();
    QString getExportCalibrateDir();

    RMat* getMasterBias();
    RMat* getMasterDark();
    RMat* getMasterFlat();
    RMat* getCannyRMat();
    RMat* getContoursRMat();
    RMat* getEllipseRMat();
    QImage* getCannyQImage();
    QList<RMat*> getContoursRMatList();
    QList<RMat*> getResultList();
    QList<RMat*> getResultList2();
    QList<RMat*> getLimbFitResultList1();
    QList<RMat*> getLimbFitResultList2();
    QList<RMat*> getLuckyBlkList();
    QVector<Circle> getCircleOutList();
    float getMeanRadius();


    // public properties (for "easier" referencing)
    QList<RMat*> rMatLightList;
    QList<RMat*> rMatBiasList;
    QList<RMat*> rMatDarkList;
    QList<RMat*> rMatFlatList;

    static bool compareContourAreas(std::vector<cv::Point> contour1, std::vector<cv::Point> contour2);

    // Display - lookup table
    void red_tab(int* red, int* green ,int* blue);
    cv::Mat scalePreviewImage(float sunX,float sunY,float sunR, cv::Mat matImage, char filter);
    void printArDims(af::array &ar);

    /// Colorize
    cv::Mat wSolarColorize(cv::Mat matImage, char filter);
    QList<RMat*> wSolarColorizeSeries(QList<RMat*> rMatImageList, char filter);

signals:

   void resultSignal(RMat* rMatResult);
   void resultSignal(cv::Mat matImage, bool bayer, instruments instrument);
   void resultSignal(cv::Mat matImage, bool bayer, instruments instrument, QString imageTitle);
   void listResultSignal(QList<RMat*> rMatListResult);
   void ellipseSignal(cv::RotatedRect ellRect);
   void resultQImageSignal(QImage &image);
   void processingQImageSignal(QImage &image);
   void cannyResultsReady();
   void cannyUpdatesReady();
   void messageSignal(QString message);
   void tempMessageSignal(QString message, int = 3000);

public slots:
   // The following are setters and slots as well.
   // They are used as slots by the RMainWindow for sending the
   // urls from the treeWidget and a drag and drop of files.
   void createMasters();
   // Preprocessing
   void calibrateOffScreen();
   //void calibrateOnScreen();
   void setupMasterWithSigmaClip(bool enabled);
   void setupMasterWithMean(bool enabled);

   void setExportMastersDir(QString dir);
   void setExportCalibrateDir(QString dir);

   /// Registration
   // Ultimately I need to use function pointers. It will get messy otherwise.
   bool prepRegistration();
   void registerSeries();
   void registerSeriesXCorrPropagate();
   void registerSeriesOnLimbFit();
   void registerSeriesByPhaseCorrelation();
   void registerSeriesCustom();
   void registerSeriesCustomPropagate();
   cv::Point calculateSADShift(cv::Mat refMat, cv::Mat matImage, cv::Rect fov, int maxLength);
   cv::Point calculateSADShift(cv::Mat refMat, cv::Mat matImage, QList<cv::Rect> fovList, int maxLength);
   cv::Mat calculateXCorrShift(cv::Mat refMat, cv::Mat matImage, cv::Rect fov);
   cv::Mat calculateXCorrShift(cv::Mat refMat, cv::Mat matImage, QList<cv::Rect> fovList);
   cv::Mat shiftToWarp(cv::Point shift);

   // Template Matching
   void registerSeriesByTemplateMatching();
   void registerSeriesByTemplateMatchingPropagate();
   cv::Point templateMatch(cv::Mat img, cv::Mat templ, int matchMethod);
   cv::Mat calculateTemplateMatchShift(cv::Mat refMat, cv::Mat matImage, cv::Rect fov);

   cv::Mat shiftImage(RMat* rMatImage, cv::Mat warpMat);
   cv::Mat shiftImage(RMat* rMatImage, cv::Point shift);
   void cannyEdgeDetectionOffScreen(int thresh);
   bool cannyEdgeDetection(int thresh);
   void setupCannyDetection(int i);
   void cannyDetect(int thresh);
   bool limbFit(int i);
   bool wernerLimbFit(QList<RMat*> rMatImageList, bool smooth, int smoothSize = 5);
   bool solarLimbRegisterSeries(QList<RMat*> rMatImageList);
   void raphFindLimb(cv::Mat matImage, Data *dat, int numDots, bool smooth, int smoothSize);

   void blurRMat(RMat* rMat);
   QList<RMat*> normalizeSeriesByStats(QList<RMat*> rMatImageList);
   RMat* normalizeByStats(RMat* rMat);
   RMat* normalizeToXposure(RMat* rMat);
   QList<RMat*> normalizeSeriesToXposure(QList<RMat*> rMatImageList);
   void normalizeByStatsInPlace(RMat* rMat);
   cv::Mat normalizeByThresh(cv::Mat matImage, float oldMin, float oldMax, float newRange);
   cv::Mat normalizeClipByThresh(cv::Mat matImage, float newMin, float newMax, float dataRange);
   void fixUset(cv::Mat matImage);
    // ROI
   void setupMaskingCircle(int circleX, int circleY, int radius);
   void appendROIList(QRect qRect);
   void clearROIs();


private:

    void normalizeFlat();
    void calibrate();

    //int circleFitLM(Data& data, Circle& circleIni, reals LambdaIni, Circle& circle);

    RListImageManager *listImageManager;
    RMat *masterBias;
    RMat *masterDark;
    RMat *masterFlat;
    RMat *masterFlatN;
    RMat *stackedRMat;
    QList<RMat*> resultList;
    QList<RMat*> resultList2;
    QList<RMat*> limbFitResultList1;
    QList<RMat*> limbFitResultList2;

    cv::Mat limbFitWarpMat;
    cv::Mat sampleMat8, sampleMatN, contoursMat;
    QImage *cannyQImage;

    RMat *cannyRMat;
    RMat *contoursRMat;
    RMat *ellipseRMat;


    QList<ImageManager*> imageManagerList;
    QList<RMat*> contoursRMatList;
    QVector<cv::RotatedRect> ellRectList;
    QVector<cv::Point2f> centers;
    QVector< std::vector< std::vector<cv::Point> > > selectedContoursList;

    RTreeWidget *treeWidget;
    bool useUrlsFromTreeWidget;
    ROpenGLWidget *currentROpenGLWidget;

    QString exportMastersDir;
    QString exportCalibrateDir;
    QDir exportQDir;

    QUrl masterBiasUrl, masterDarkUrl, masterFlatUrl;

    bool biasSuccess, darkSuccess, flatSuccess;
    bool showContours, showLimb;
    bool useXCorr;
    bool masterWithMean, masterWithSigmaClip;
    bool stackWithMean, stackWithSigmaClip;

    float radius, radius1, radius2, radius3, meanRadius;
    // ROI
    bool useROI;
    bool applyMask;
    QList<cv::Rect> cvRectROIList;
    cv::Rect cvRectROI;
    int maskCircleX, maskCircleY, maskCircleRadius;


    Circle circleOut;
    QVector<Circle> circleOutList;

    // plots
    QCustomPlot *limbFitPlot;

    /// Canny parameters
    /// Blur
    double blurSigma;

    // Fourier Filters
    bool useHPF;
    double hpfSigma;

    // Sharpenning
    bool sharpenLiveStatus;

    // Block Processing (lucky imaging)
    int blkSize;
    int binning;
    int nBest;
    QList<RMat*> luckyBlkList;
    QString qualityMetric;

    // Normalization of the images
    double normFactor;


};

#endif // RPROCESSING_H
