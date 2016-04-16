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


class RProcessing: public QObject
{
    Q_OBJECT

public:
    RProcessing(QObject *parent = 0);
    ~RProcessing();

    void loadRMatLightList(QList<QUrl> urls);
    void loadRMatBiasList(QList<QUrl> urls);
    void loadRMatDarkList(QList<QUrl> urls);
    void loadRMatFlatList(QList<QUrl> urls);

    // Calibration (bias, dark, ...)
    bool makeMasterBias();
    bool makeMasterDark();
    bool makeMasterFlat();

    RMat *average(QList<RMat*> rMatList);

    // export methods
    void exportMastersToFits();
    void exportToFits(RMat *rMatImage, QString QStrFilename);
    QString setupFileName(QFileInfo fileInfo, QString format);
    void loadMasterDark();
    void loadMasterFlat();

    void showMinMax(const cv::Mat & matImage);
    // setters
    void setTreeWidget(RTreeWidget *treeWidget);
    void setCurrentROpenGLWidget(ROpenGLWidget *rOpenGLWidget);
    void setShowContours(bool status);
    void setShowLimb(bool status);
    void setUseXCorr(bool useXCorr);
    void setCvRectROI(cv::Rect cvRect);
    void setUseROI(bool status);

    // getters
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
    QList<RMat*> getLimbFitPreviewList();
    QVector<Circle> getCircleOutList();

    // public properties (for "easier" referencing)
    QList<RMat*> rMatLightList;
    QList<RMat*> rMatBiasList;
    QList<RMat*> rMatDarkList;
    QList<RMat*> rMatFlatList;

    static bool compareContourAreas(std::vector<cv::Point> contour1, std::vector<cv::Point> contour2);

signals:

   void resultSignal(RMat* rMatResult);
   void resultSignal(cv::Mat matImage, bool bayer, instruments instrument);
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


   void setExportMastersDir(QString dir);
   void setExportCalibrateDir(QString dir);
   void registerSeries();
   void registerSeries(QList<RMat*> rMatList);
   void cannyEdgeDetectionOffScreen(int thresh);
   void cannyEdgeDetection(int thresh);
   void setupCannyDetection(int i);
   void cannyDetect(int thresh);
   void limbFit(int i);
   void cannyRegisterSeries();
   void blurRMat(RMat* rMat);
   QList<RMat*> normalizeSeriesByStats(QList<RMat*> rMatImageList);
   RMat* normalizeByStats(RMat* rMat);
   void normalizeByStatsInPlace(RMat* rMat);
   cv::Mat normalizeClipByThresh(RMat* rMat, float newMin, float newMax);

private:

    void normalizeFlat();
    void calibrate();
    //int circleFitLM(Data& data, Circle& circleIni, reals LambdaIni, Circle& circle);

    RListImageManager *listImageManager;
    RMat *masterBias;
    RMat *masterDark;
    RMat *masterFlat;
    RMat *masterFlatN;
    QList<RMat*> resultList;
    QList<RMat*> resultList2;
    cv::Mat sampleMat8, sampleMatN, contoursMat;
    QImage *cannyQImage;

    RMat *cannyRMat;
    RMat *contoursRMat;
    RMat *ellipseRMat;

    QList<ImageManager*> imageManagerList;
    QList<RMat*> contoursRMatList;
    QList<RMat*> limbFitPreviewList;
    QVector<cv::RotatedRect> ellRectList;
    QVector<cv::Point2f> centers;
    QVector< std::vector< std::vector<cv::Point> > > biggestContoursList;

    RTreeWidget *treeWidget;
    ROpenGLWidget *currentROpenGLWidget;

    QString exportMastersDir;
    QString exportCalibrateDir;

    QString masterBiasPath, masterDarkPath, masterFlatPath;

    bool biasSuccess, darkSuccess, flatSuccess;
    bool showContours, showLimb;
    bool useXCorr;

    float radius, radius1, radius2, radius3;
    cv::Rect cvRectROI;
    bool useROI;
    Circle circleOut;
    QVector<Circle> circleOutList;

    // plots
    QCustomPlot *limbFitPlot;

};

#endif // RPROCESSING_H
