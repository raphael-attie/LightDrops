#ifndef RPROCESSING_H
#define RPROCESSING_H

#include "winsockwrapper.h"
#include <QApplication>
#include <QtCore>

#include "rmat.h"
#include "rlistimagemanager.h"
#include "rtreewidget.h"
#include "rlineedit.h"
#include "ropenglwidget.h"

#include "typedefs.h"
#include "data.h"
#include "circle.h"
#include "utilities.h"
#include "werner/limb.h"

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

    /// export methods
    void exportMastersToFits();
    void exportToFits(RMat *rMatImage, QString QStrFilename);
    QString setupFileName(QFileInfo fileInfo, QString format);
    void loadMasterDark();
    void loadMasterFlat();

    /// Statistics
    void showMinMax(const cv::Mat & matImage);
    cv::Mat histogram(cv::Mat matVector, int &nBins, float &width);
    float calcMedian(std::vector<float> data, float width);


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
   void registerSeries();
   void registerSeriesOnLimbFit();
   void registerSeriesByPhaseCorrelation();
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
   void normalizeByStatsInPlace(RMat* rMat);
   cv::Mat normalizeByThresh(cv::Mat matImage, float oldMin, float oldMax, float newRange);
   cv::Mat normalizeClipByThresh(cv::Mat matImage, float newMin, float newMax, float dataRange);
   void fixUset(cv::Mat matImage);

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
    ROpenGLWidget *currentROpenGLWidget;

    QString exportMastersDir;
    QString exportCalibrateDir;

    QString masterBiasPath, masterDarkPath, masterFlatPath;

    bool biasSuccess, darkSuccess, flatSuccess;
    bool showContours, showLimb;
    bool useXCorr;
    bool masterWithSigmaClip;
    bool masterWithMean;
    bool stackWithMean, stackWithSigmaClip;

    float radius, radius1, radius2, radius3, meanRadius;
    cv::Rect cvRectROI;
    bool useROI;
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

};

#endif // RPROCESSING_H
