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
    void messageSignal(QString message);
    void plotSignal(QCustomPlot *customPlot);

    void radioLightImages(QList<RMat*> rMatImageList);
    void radioBiasImages(QList<RMat*> rMatImageList);
    void radioDarkImages(QList<RMat*> rMatImageList);
    void radioFlatImages(QList<RMat*> rMatImageList);

protected:
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);


private slots:

    void createNewImage(RListImageManager *newRListImageManager);
    void createNewImage(QList<RMat*> newRMatImageList);
    void createNewImage(RMat *rMatImage);
    void addImage(ROpenGLWidget *rOpenGLWidget);
    void updateDoubleSpinBox(int);
    void scaleImageSlot(int value);
    void gammaScaleImageSlot(int value);
    void updateSliderValueSlot(double valueD);
    void updateWB(int value);
    void autoScale();
    void minMaxScale();
    void updateCurrentData(ROpenGLWidget *rOpenGLWidget);
    void updateFrameInSeries(int frameIndex);
    void updateSubFrame(QImage *image, float intensity, int x, int y);
    void addPlotWidget(QCustomPlot *plotWidget);
    void addHeaderWidget();
    void setupExportCalibrateDir();
    void exportMastersToFits();
    void exportFrames();


    //sliderFrame playback buttons
    void stepForward();
    void stepBackward();
    void playForward();
    void increaseFps();
    void decreaseFps();

    void stopButtonPressed();

    void radioRMatSlot();

    void on_actionHeader_triggered();

    void on_actionHeader_toggled(bool arg1);

    void uncheckActionHeaderState();

private:

    // functions
    void resizeOpenGLWidget(QScrollArea *scrollArea);
    void loadSubWindow(QScrollArea *scrollArea);
    void dispatchRMatImages(QList<RMat*> rMatList);
    float convertSliderToScale(int value);
    int convertScaleToSlider(float value);
    float convertSliderToGamma(int value);
    int convertGammaToSlider(float gamma);
    void setupSliders();
    void setupSubImage();
    void drawHist(cv::Mat matHist);
    void updateCustomPlotLineItems();

    // properties
    Ui::RMainWindow *ui;
    RProcessing *processing;
    QImage *subQImage;
    QCustomPlot *customPlot;
    QCPItemLine *vertLineHigh;
    QCPItemLine *vertLineLow;
    QString checkExistingDir();

    void tileView();

    ROpenGLWidget *currentROpenGLWidget;

    QSize oglSize;
    QSize defaultWindowSize;
    float sliderScale, sliderRange, sliderScaleWB;

    float gammaScale, newMax, newMin, alpha, beta, gammaMin, gamma;

    int doubleSpinBoxDecimals;
    int sliderValueHigh, sliderValueLow, sliderValueGamma;
    int scrollBarHeight;

    bool stopButtonStatus;


};

#endif // RMAINWINDOW_H
