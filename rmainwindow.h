#ifndef RMAINWINDOW_H
#define RMAINWINDOW_H

#include "winsockwrapper.h"
#include <QtCore>
#include <QMainWindow>
#include <QTreeWidgetItem>

#include "rmat.h"
#include "imagemanager.h"
#include "rlistimagemanager.h"
#include "ropenglwidget.h"
#include "rprocessing.h"
#include "RFrame.h"
#include "rgraphicsscene.h"

//QCustomPlot
#include <qcustomplot/qcustomplot.h>

namespace Ui {
class RMainWindow;
}

class RMainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit RMainWindow(QWidget *parent = 0);
    ~RMainWindow();

signals:
    void tempMessageSignal(QString message, int = 3000);
    void plotSignal(QCustomPlot *customPlot);

    void radioLightImages(QList<RMat*> rMatImageList);
    void radioBiasImages(QList<RMat*> rMatImageList);
    void radioDarkImages(QList<RMat*> rMatImageList);
    void radioFlatImages(QList<RMat*> rMatImageList);

protected:
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);
    void closeEvent(QCloseEvent *event);


private slots:

    void createNewImage(RListImageManager *newRListImageManager);
    void createNewImage(QList<RMat*> newRMatImageList);
    void createNewImage(RMat *rMatImage);
    void createNewImage(cv::Mat cvImage, bool bayer = false, instruments instrument = instruments::generic, QString imageTitle = QString("Result"));
    void createNewImage(QImage &image, ROpenGLWidget *rOpenGLWidget = NULL, bool inverted = true);

    void displayQImage(QImage &image, RGraphicsScene *scene, QMdiSubWindow *subWindow, QString windowTitle = QString("Processing Window"));
    void initPreviewQImage(bool status);

    void selectROI(bool isSquare = false, int blkSize = 0);
    void setRect(QRect rect);
    void extractNewImageROI();
    void disableROIaction();

    //TreeWidget
    //void removeTreeWidgetItems(QTreeWidgetItem* item);

    //void addEllipseToScene(cv::RotatedRect rect);

    //void createNewImage(QImage *image);
    void addImage(ROpenGLWidget *rOpenGLWidget);
    void updateDoubleSpinBox(int);
    void scaleImageSlot(int value);
    void gammaScaleImageSlot(int value);
    void updateSliderValueSlot();
    void updateWB(int value);
    void autoScale();
    void autoScale(ROpenGLWidget *rOpenGLWidget);
    void minMaxScale();
    void rangeScale();
    void changeROpenGLWidget(ROpenGLWidget *rOpenGLWidget);
    void updateFrameInSeries(int frameIndex);
    void updateSubFrame(QImage *image, float intensity, int x, int y);
    void displayPlotWidget(ROpenGLWidget *rOpenGLWidget);
    void addHeaderWidget();
    void setupExportCalibrateDir();
    void exportMastersToFits();
    void exportFrames();
    QString makeFilePath(QString basename, int frameNumber);
    void exportFramesToFits();
    void exportFramesToJpeg();
    void exportFramesToTiff();
    void convertTo8Bit();
    cv::Mat convertTo8Bit(RMat *rMatImage);
    void convertToNeg();
    void solarColorizeSeriesSlot();
    void initSharpenImageSlot(bool status);
    void sharpenImageSlot();
    void sharpenSliderSlot();

    ///  Graphs
    void changeZoomAxisSlot();
    QVector<double> extractTemperatureFromSeries();
    QVector<double> extractMeanFromSeries();
    void displayTemperatureSeries();
    void displayMeanSeries();
    void plotData(QVector<double> data, QString xLabel, QString yLabel);

    // Processing
    void calibrateOffScreenSlot();
    void registerSlot();
    void solarLimbFit();
    void solarLimbRegisterSlot();
    void updateLimbResults();
    void normalizeCurrentSeries();
    void convert14to16bitSeriesSlot();
    void previewMatImageHPFSlot();
    void stackSlot();
    void binningSlot();
    void HDRSlot();

    // Lucky image
    void blockProcessingSlot();
    void extractLuckySampleSlot();
    void luckyROISlot();
    void setupLuckyImaging();


    // Plots
    void showLimbFitStats();

    // ToneMapping
    void setupToneMappingCurve();
    void updateToneMappingSlot();
    void applyToneMappingCurve();
    void hAlphaToneMappingSlot();
    void applyScaleLimbSlot();

    //sliderFrame playback buttons
    void stepForward();
    void stepBackward();
    void playForward();
    void increaseFps();
    void decreaseFps();

    //Statistics
    void updateStats(int frameNumber);



    void stopButtonPressed();

    void radioRMatSlot();

    void on_actionHeader_toggled(bool arg1);

    void uncheckActionHeaderState();

    void on_actionTileView_triggered();

    void on_actionROISelect_triggered();

    void on_actionROIExtract_triggered();

    void on_actionHeader_triggered();

    void on_actionTemperature_triggered();

    void on_actionTemperature_toggled(bool arg1);

    void on_actionROISelect_toggled(bool arg1);

    void on_actionCircle_toggled(bool arg1);

    void on_actionMultiROI_toggled(bool arg1);

    void on_actionSave_ROI_triggered();

    void on_actionClear_ROIs_triggered();

    void on_actionplot_metadata_triggered();

private:

    // functions
    void updateCurrentROpenGLWidget();
    void resizeScrollArea(ROpenGLWidget *rOpenGLWidget, QScrollArea *scrollArea);
    void loadSubWindow(QScrollArea *scrollArea);
    void dispatchRMatImages(QList<RMat*> rMatList);
    float convertSliderToScale(int value);
    float convertLimbSliderToScale(int value);
    int convertScaleToSlider(float value);
    float convertSliderToGamma(int value);
    int convertGammaToSlider(float gamma);
    void setupSliders(ROpenGLWidget *rOpenGLWidget);
    void setupSubImage();
    void updateCustomPlotLineItems();
    void updateInvGaussianParams();



    // properties
    //cv::Mat ellMat;
//    QImage image2;
    //QGraphicsPixmapItem *item;
    RGraphicsScene* scene;
    QGraphicsView* graphicsView;
    QGraphicsPixmapItem *pixMapItem;

    Ui::RMainWindow *ui;
    RProcessing *processing;
    QImage *subQImage;
    QCustomPlot *customPlot;
    QCPItemLine *vertLineHigh;
    QCPItemLine *vertLineLow;

    QCustomPlot *toneMappingPlot;
    QCPGraph *toneMappingGraph;

    QString checkExistingDir();

    void tileView();

    ROpenGLWidget *currentROpenGLWidget;
    ROpenGLWidget *lastROpenGLWidget;
    ROpenGLWidget *resultROpenGLWidget;
    ROpenGLWidget *cannyContoursROpenGLWidget;
    ROpenGLWidget *limbFittingROpenGLWidget;

    QMdiSubWindow *limbRegisterSubWindow;
    QMdiSubWindow *currentSubWindow;
    QMdiSubWindow *plotSubWindow;
    QMdiSubWindow *tempSubWindow;
    QMdiSubWindow *previewSubWindow;
    QScrollArea *currentScrollArea;
    QScrollArea *cannyScrollArea;
    QScrollArea *processingScrollArea;

    QSize oglSize;
    QSize defaultWindowSize;
    float sliderScale, sliderRange, sliderScaleWB, sliderToScaleMinimum;
    float limbSliderScale, limbSliderRange, limbSliderToScaleMinimum;

    float gammaScale, alpha, beta, gammaMin, gamma;
    double iMax, lambda, mu;

    int doubleSpinBoxDecimals;
    int sliderValueHigh, sliderValueLow, sliderValueGamma;
    int scrollBarHeight;
    int frameIndex;

    bool stopButtonStatus;

    QRect rect;
    //QCustomPlot *limbFitPlot;

    /// image objects used for previews
    cv::Mat matImageHPFPreview;
    QImage previewQImage;

    /// Limb fitting temporary results
    QVector<Circle> fittedLimbList;  // With smooth
    QVector<Circle> fittedLimbList2; // Without smooth

    QCustomPlot *limbFitPlot;

};

#endif // RMAINWINDOW_H
