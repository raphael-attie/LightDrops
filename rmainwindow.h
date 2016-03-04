#ifndef RMAINWINDOW_H
#define RMAINWINDOW_H

#include <QtCore>
#include <QMainWindow>
#include <QTreeWidgetItem>

#include "scrollarea.h"
#include "rmat.h"
#include "imagemanager.h"
#include "rlistimagemanager.h"
#include "ropenglwidget.h"
#include "rprocessing.h"
#include "RFrame.h"
#include "rtableworker.h"

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

protected:
    void dropEvent(QDropEvent* event);
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);


private slots:

    void createNewImage(RListImageManager *newRListImageManager);
    void createNewImage(QList<RMat> &newRMatImageList);
    void createNewImage(RMat &newRMatImage, QString windowTitle = QString(" "));
    void addImage(ROpenGLWidget *rOpenGLWidget);
    void updateDoubleSpinBox(int);
    void scaleImageSlot(int value);
    void gammaScaleImageSlot(int value);
    void updateSliderValueSlot(double valueD);
    void updateWB(int value);
    void autoScale();
    void updateCurrentData(ROpenGLWidget *rOpenGLWidget);
    void updateFrameInSeries(int frameIndex);
    void updateSubFrame(QImage *image, float intensity, int x, int y);
    void updateTableWidget(int value);
    void addPlotWidget(QCustomPlot *plotWidget);
    void addHeaderWidget();

    void setupExportMastersDir();
    void setupExportCalibrateDir();

    //sliderFrame playback buttons
    void stepForward();
    void stepBackward();
    void playForward();

    void stopButtonPressed();

private:

    // functions
    void loadSubWindow(QScrollArea *scrollArea);
    float convertSliderToScale(int value);
    int convertScaleToSlider(float value);
    float convertSliderToGamma(int value);
    int convertGammaToSlider(float gamma);
    void setupSliders();
    void setupSubImage();
    void drawHist(cv::Mat matHist);


    // properties
    Ui::RMainWindow *ui;
    RProcessing *processing;
    QImage *subQImage;
    QCustomPlot *customPlot;
    QCPItemLine *vertLineHigh;
    QCPItemLine *vertLineLow;
    QTableWidget *headerWidget;
    QThread *tableThread;
    RTableWorker *rTableWorker;
    QString checkExistingDir();

    void tileView();

    ROpenGLWidget *currentROpenGLWidget;
    QList<RMat> currentRMatList;

    QSize oglSize;
    float sliderScale, sliderRange, sliderScaleWB;

    float gammaScale, newMax, newMin, alpha, beta, gammaMin, gamma;

    int doubleSpinBoxDecimals;
    int sliderValueHigh, sliderValueLow, sliderValueGamma;

    bool stopButtonStatus;


};

#endif // RMAINWINDOW_H
