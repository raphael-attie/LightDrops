#include <QMdiSubWindow>
#include <QDropEvent>
#include <QFileDialog>

#include "rmainwindow.h"
#include "ui_rmainwindow.h"
#include "scrollarea.h"

//QCustomPlot
#include <qcustomplot/qcustomplot.h>


RMainWindow::RMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::RMainWindow),
    currentROpenGLWidget(NULL)
{
    ui->setupUi(this);
    setCentralWidget(ui->mdiArea);
    processing = new RProcessing(this);
    processing->setTreeWidget(ui->treeWidget);
    vertLineHigh = NULL;
    headerWidget = NULL;

    sliderScale = 1.0;
    gammaScale = 0.1;
    gammaMin = ui->doubleSpinBoxGamma->minimum();
    gamma = 1.0f;
    sliderScaleWB = 0.01;

    setupSubImage();

//    tabifyDockWidget(ui->currentFrameDock, ui->settingsDock);
//    this->setTabPosition(Qt::DockWidgetArea::BottomDockWidgetArea, QTabWidget::TabPosition::South);

    //this->setTabShape(QTabWidget::Triangular);

    // Connect the change of slider value to the number displayed in the spinBox
    connect(ui->sliderHigh, SIGNAL(valueChanged(int)), this, SLOT(updateDoubleSpinBox(int)));
    connect(ui->sliderLow, SIGNAL(valueChanged(int)), this, SLOT(updateDoubleSpinBox(int)));
    connect(ui->sliderGamma, SIGNAL(valueChanged(int)), this, SLOT(updateDoubleSpinBox(int)));

    // Connect the change of the number in the spinBox to the change of slider value
    connect(ui->doubleSpinBoxHigh, SIGNAL(valueChanged(double)), this, SLOT(updateSliderValueSlot(double)));
    connect(ui->doubleSpinBoxLow, SIGNAL(valueChanged(double)), this, SLOT(updateSliderValueSlot(double)));
    connect(ui->doubleSpinBoxGamma, SIGNAL(valueChanged(double)), this, SLOT(updateSliderValueSlot(double)));

    // Connect the change of slider value to the image linear scaling and gamma scaling.
    connect(ui->sliderHigh, SIGNAL(valueChanged(int)), this, SLOT(scaleImageSlot(int)));
    connect(ui->sliderLow, SIGNAL(valueChanged(int)), this, SLOT(scaleImageSlot(int)));
    connect(ui->sliderGamma, SIGNAL(valueChanged(int)), this, SLOT(gammaScaleImageSlot(int)));

    // Connect the white balance sliders
    connect(ui->redSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWB(int)));
    connect(ui->blueSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWB(int)));
    connect(ui->greenSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWB(int)));

    // Connect the pushButton and lineEdit of the ui for creating the masters in the processing class
    connect(ui->mastersDirLineEdit, SIGNAL(textChanged(QString)), processing, SLOT(setExportMastersDir(QString)));
    connect(ui->exportMastersDirButton, SIGNAL(released()), this, SLOT(setupExportMastersDir()));
    connect(ui->makeMasterButton, SIGNAL(pressed()), processing, SLOT(createMasters()));

    connect(ui->calibrateDirLineEdit, SIGNAL(textChanged(QString)), processing, SLOT(setExportCalibrateDir(QString)));

    // Connect the sliderFrame and playback buttons
    connect(ui->sliderFrame, SIGNAL(sliderMoved(int)), this, SLOT(updateFrameInSeries(int)));
    connect(ui->forwardButton, SIGNAL(released()), this, SLOT(stepForward()));
    connect(ui->backwardButton, SIGNAL(released()), this, SLOT(stepBackward()));
    connect(ui->playButton, SIGNAL(released()), this, SLOT(playForward()));
    connect(ui->stopButton, SIGNAL(released()), this, SLOT(stopButtonPressed()));

    // Connect the header button for header display
    connect(ui->headerPushButton, SIGNAL(released()), this, SLOT(addHeaderWidget()));
    connect(processing, SIGNAL(tempMessageSignal(QString,int)), this->statusBar(), SLOT(showMessage(QString,int)));
}

RMainWindow::~RMainWindow()
{
    delete ui;
    delete processing;
}

void RMainWindow::dropEvent(QDropEvent *event)
{
    const QMimeData* mimeData = event->mimeData();

    if (!mimeData->hasUrls())
        return;

    event->acceptProposedAction();

    RListImageManager *newRListImageManager = new RListImageManager(mimeData->urls());
    if (newRListImageManager->rMatImageList.empty())
    {
        qDebug("Nothing to load: wrong file extention?");
        return;
    }

    createNewImage(newRListImageManager);

}

void RMainWindow::dragEnterEvent(QDragEnterEvent *event)
{
    if (event->mimeData()->hasUrls())
    {
        event->acceptProposedAction();
    }
}

void RMainWindow::dragMoveEvent(QDragMoveEvent *event)
{
    if (event->mimeData()->hasUrls())
    {
        event->acceptProposedAction();
    }
}

void RMainWindow::createNewImage(RListImageManager *newRListImageManager)
{
    currentRMatList = newRListImageManager->rMatImageList;
    this->currentROpenGLWidget = new ROpenGLWidget(newRListImageManager, this);
    addImage(currentROpenGLWidget);
}

void RMainWindow::createNewImage(QList<RMat> &newRMatImageList)
{
    currentRMatList = newRMatImageList;
    this->currentROpenGLWidget = new ROpenGLWidget(newRMatImageList, this);
    addImage(currentROpenGLWidget);
}

void RMainWindow::createNewImage(RMat &newRMatImage, QString windowTitle)
{
    currentRMatList.clear();
    currentRMatList << newRMatImage;
    this->currentROpenGLWidget = new ROpenGLWidget(newRMatImage);
    addImage(currentROpenGLWidget);
}

void RMainWindow::addImage(ROpenGLWidget *rOpenGLWidget)
{
    ui->sliderFrame->setRange(0, currentRMatList.size()-1);
    ui->sliderFrame->setValue(0);
    ui->imageLabel->setText(QString("1") + QString("/") + QString::number(currentRMatList.size()));

    QScrollArea *scrollArea = new QScrollArea();
    scrollArea->setWidget(rOpenGLWidget);
    loadSubWindow(scrollArea);

    connect(rOpenGLWidget, SIGNAL(gotSelected(ROpenGLWidget*)), this, SLOT(updateCurrentData(ROpenGLWidget*)));
    connect(rOpenGLWidget, SIGNAL(sendSubQImage(QImage*,float,int,int)), this, SLOT(updateSubFrame(QImage*,float,int,int)));
    connect(this, SIGNAL(plotSignal(QCustomPlot*)), this, SLOT(addPlotWidget(QCustomPlot*)));

    newMax = (float) currentRMatList.at(0).getDataMax();
    newMin = (float) currentRMatList.at(0).getDataMin();

    drawHist(currentRMatList.at(0).getMatHist());
    setupSliders();

    autoScale();
}


void RMainWindow::loadSubWindow(QScrollArea *scrollArea)
{
    QMdiSubWindow *newSubWindow = new QMdiSubWindow;
    newSubWindow->setWidget(scrollArea);
    ui->mdiArea->addSubWindow(newSubWindow);

    oglSize = currentROpenGLWidget->getOglDefaultSize();
    currentROpenGLWidget->resize(oglSize);

    qDebug() << "RMainWindow::loadSubWindow() oglSize=" << oglSize;

    if (ui->mdiArea->subWindowList().size() > 2)
    {
        tileView();
    }
    else
    {
        // Resize the subWindow.
        QSize subWindowSize = oglSize;
        subWindowSize.rwidth() += 5;
        subWindowSize.rheight() += 25;
        newSubWindow->resize(subWindowSize);
    }
    newSubWindow->show();
}



void RMainWindow::updateDoubleSpinBox(int value)
{
    if (this->sender() == ui->sliderHigh)
    {
        float valueF = convertSliderToScale(value);
        ui->doubleSpinBoxHigh->setValue(valueF);
    }
    else if (this->sender() == ui->sliderLow)
    {
        float valueF = convertSliderToScale(value);
        ui->doubleSpinBoxLow->setValue(valueF);
    }
    if (this->sender() == ui->sliderGamma)
    {
        float valueF = convertSliderToGamma(value);
        ui->doubleSpinBoxGamma->setValue(valueF);
    }


}

void RMainWindow::scaleImageSlot(int value)
{
    if (currentROpenGLWidget == NULL || currentRMatList.empty())
        return;

    if (this->sender() == ui->sliderHigh)
    {
        newMax = convertSliderToScale(value);
    }
    else if (this->sender() == ui->sliderLow)
    {
        newMin = convertSliderToScale(value);
    }
    else
    {
        return;
    }

    if (newMin == newMax)
    {
        return;
    }

    currentROpenGLWidget->setNewMax(newMax);
    currentROpenGLWidget->setNewMin(newMin);
    float dataRange = (float) (newMax - newMin) ;
    alpha = 1.0f / dataRange;
    beta = (float) (-newMin / dataRange);
    qDebug("RMainWindow::scaleImageSlot() value = %i", value);
    qDebug("RMainWindow::scaleImageSlot() newMax = %f , newMin = %f", newMax, newMin);
    currentROpenGLWidget->setAlpha(alpha);
    currentROpenGLWidget->setBeta(beta);
    currentROpenGLWidget->updateSubQImage();
    currentROpenGLWidget->update();

    double nPixels = (double) currentRMatList.at(0).getNPixels();
//    vertLineHigh = new QCPItemLine(customPlot);
    vertLineHigh->start->setCoords(newMax, 0);
    vertLineHigh->end->setCoords(newMax, nPixels);
    vertLineLow->start->setCoords(newMin, 0);
    vertLineLow->end->setCoords(newMin, nPixels);

    customPlot->replot();

}

void RMainWindow::setupSliders()
{
    // The slider range of integer values must be consistent with the maximum data range.
    // For USET for example, it is 4096.

    if (currentRMatList.at(0).getInstrument() == instruments::USET)
    {
        qDebug("RMainWindow::setupSliders():: setting sliders for USET data.");
        ui->sliderHigh->setRange(1, 4096);
        ui->sliderLow->setRange(1, 4096);
    }
    else
    {
        ui->sliderHigh->setRange(1, 65536);
        ui->sliderLow->setRange(1, 65536);
    }
    sliderRange = (float) ui->sliderHigh->maximum() - ui->sliderHigh->minimum() + 1;

    if (currentRMatList.at(0).getBscale() == 1 && currentRMatList.at(0).getInstrument() != instruments::MAG)
    {
        sliderScale = 1.0;
    }
    else
    {
        float dataMax = (float) currentRMatList.at(0).getDataMax();
        float dataMin = (float) currentRMatList.at(0).getDataMin();
        float dataRange = dataMax - dataMin;
        sliderScale = (float) dataRange / sliderRange;
        qDebug() << "dataMax =" << dataMax;
        qDebug() << "dataMin =" << dataMin;
        qDebug() << "sliderScale =" << sliderScale;
    }
    ui->doubleSpinBoxHigh->setMaximum(convertSliderToScale(ui->sliderHigh->maximum()));
    ui->doubleSpinBoxHigh->setMinimum(convertSliderToScale(ui->sliderHigh->minimum()));
    ui->doubleSpinBoxLow->setMaximum(convertSliderToScale(ui->sliderLow->maximum()));
    ui->doubleSpinBoxLow->setMinimum(convertSliderToScale(ui->sliderLow->minimum()));
}

void RMainWindow::setupSubImage()
{
    int subNaxis = 50;
    uint subMagnification = 3;
    subQImage = new QImage(subNaxis, subNaxis, QImage::Format_ARGB32);
    subQImage->fill(QColor(125, 125, 125, 125));

    ui->subFrame->setImage(subQImage);
    ui->subFrame->setDrawCross(true);
    ui->subFrame->update();

    //subPaintWidget = new PaintWidget(subQImage, subMagnification);
    //subPaintWidget->setDrawCross(true);


}

void RMainWindow::drawHist(cv::Mat matHist)
{
    cv::Mat tempHist;
    matHist.convertTo(tempHist, CV_64F);

    int nBins = matHist.rows;
    double minHist;
    double maxHist;
    cv::minMaxLoc(matHist, &minHist, &maxHist);
    double histWidth = currentRMatList.at(ui->sliderFrame->value()).getHistWidth();

    double dataMax = currentRMatList.at(ui->sliderFrame->value()).getDataMax();
    double dataMin = currentRMatList.at(ui->sliderFrame->value()).getDataMin();

    QVector<double> x(nBins); // y must be the histogram values
    for (int i=0; i< nBins; ++i)
    {
      x[i] = dataMin + histWidth*i; // x goes from 0 to nBins-1
    }


    // Pointer to matHist data.
    const double* matHistPtr = tempHist.ptr<double>(0);
    std::vector<double> matHistStdVect(matHistPtr, matHistPtr + matHist.rows);
    QVector<double> y = QVector<double>::fromStdVector(matHistStdVect);

    // create graph and assign data to it:
    customPlot = new QCustomPlot();
    customPlot->addGraph();
    customPlot->plottable(0)->setPen(QPen(QColor(125, 125, 125, 50))); // line color gray
    customPlot->plottable(0)->setBrush(QBrush(QColor(125, 125, 125, 50)));


    customPlot->xAxis->setTicks(false);
    customPlot->yAxis->setTicks(false);

    // make left and bottom axes always transfer their ranges to right and top axes:
    connect(customPlot->xAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->xAxis2, SLOT(setRange(QCPRange)));
    connect(customPlot->yAxis, SIGNAL(rangeChanged(QCPRange)), customPlot->yAxis2, SLOT(setRange(QCPRange)));

    customPlot->yAxis->setScaleType(QCPAxis::stLogarithmic);
    customPlot->graph(0)->setData(x, y);
    customPlot->rescaleAxes();

    double minXRange;
    double maxXRange;

    if (currentRMatList.at(0).getInstrument() == instruments::USET)
    {
        minXRange = 0.0;
        maxXRange = 4095.0;
    }
    else
    {
        minXRange = std::min(dataMin, 0.0);
        maxXRange = dataMax;
    }

    qDebug("RMainWindow::drawHist():: maxXRange = %f", maxXRange);
    customPlot->xAxis->setRange(1.05*minXRange, 1.05*maxXRange);
    customPlot->yAxis->setRange(1, maxHist);

    customPlot->setInteraction(QCP::iRangeDrag, true);
    customPlot->setInteraction(QCP::iRangeZoom, true);
    customPlot->axisRect(0)->setRangeDrag(Qt::Horizontal);
    customPlot->axisRect(0)->setRangeZoom(Qt::Horizontal);

    double nPixels = (double) currentRMatList.at(0).getNPixels();
    vertLineHigh = new QCPItemLine(customPlot);
    vertLineLow = new QCPItemLine(customPlot);
    customPlot->addItem(vertLineHigh);
    customPlot->addItem(vertLineLow);

    vertLineHigh->start->setCoords(dataMax, 0);
    vertLineHigh->end->setCoords(dataMax, nPixels);
    vertLineLow->start->setCoords(dataMin, 0);
    vertLineLow->end->setCoords(dataMin, nPixels);

//    QCPItemText* itemText = new QCPItemText(customPlot);
//    customPlot->addItem(itemText);

    qDebug("customPlot ready for emit signal.");
    emit plotSignal(customPlot);
}

void RMainWindow::addPlotWidget(QCustomPlot *plotWidget)
{
    plotWidget->resize(200, 100);
    ui->histoDock->setWidget(plotWidget);


}

void RMainWindow::addHeaderWidget()
{
    int frameIndex = ui->sliderFrame->value();

    if (currentROpenGLWidget == NULL)
    {
        qDebug("RMainWindow::addHeaderWidget():: currentROpenGLWidget == NULL");
        return;
    }
    else if (currentROpenGLWidget->getRListImageManager() == NULL)
    {
      qDebug("RMainWindow::addHeaderWidget():: currentROpenGLWidget->getRListImageManager() == NULL");
      return;
    }
    else if (currentRMatList.at(frameIndex).isBayer())
    {
        qDebug("RMainWindow::addHeaderWidget():: current image is Bayer type");
        return;
    }


    qDebug("Creating headerWidget");
    tableThread = new QThread;
    rTableWorker = new RTableWorker(currentROpenGLWidget->getRListImageManager());
    rTableWorker->moveToThread(tableThread);

    headerWidget = rTableWorker->getNewTableWidget();

    QString windowTitle = currentROpenGLWidget->getWindowTitleList().at(frameIndex);
    headerWidget->setWindowTitle(windowTitle);
    headerWidget->setAttribute(Qt::WA_DeleteOnClose);

    QMdiSubWindow *newSubWindow = new QMdiSubWindow;
    newSubWindow->setWidget(headerWidget);
    newSubWindow->setAttribute(Qt::WA_DeleteOnClose);
    newSubWindow->setWindowTitle(windowTitle);

    ui->mdiArea->addSubWindow(newSubWindow);
    newSubWindow->show();

    // Connect the thread start/finish process to the worker.
    QObject::connect(tableThread, SIGNAL(started()), rTableWorker, SLOT(process()));
    QObject::connect(rTableWorker, SIGNAL(finished()), tableThread, SLOT(quit()));
    //Connect the slider value changed with the header dock change
    QObject::connect(ui->sliderFrame, SIGNAL(valueChanged(int)), this, SLOT(updateTableWidget(int)));
}

void RMainWindow::updateSubFrame(QImage *image, float intensity, int x, int y)
{
    char format = 'f';
    int precision = 0;

    if (std::abs(intensity) < 0.001)
    {
        format = 'g';
        precision = 2;
    }

    QString intensityQStr = QString::number(intensity, format, precision);

    ui->subFrame->setImage(image);
    ui->subFrame->setCursorText(intensityQStr);
    ui->subFrame->setImageCoordX(x);
    ui->subFrame->setImageCoordY(y);
    ui->subFrame->update();
}

void RMainWindow::updateTableWidget(int value)
{
    if (currentROpenGLWidget->rMatImageList.at(value).isBayer())
    {
      return;
    }

    rTableWorker->setFrameIndex(value);
    tableThread->start();

    QString windowTitle = currentROpenGLWidget->getWindowTitleList().at(value);
    headerWidget->setWindowTitle(windowTitle);
}


QString RMainWindow::checkExistingDir()
{
    QString openDir = QString();
    if (!ui->mastersDirLineEdit->text().isEmpty())
    {
        openDir = ui->mastersDirLineEdit->text();
    }
    else if (!ui->calibrateDirLineEdit->text().isEmpty())
    {
        openDir = ui->calibrateDirLineEdit->text();
    }
    else
    {
        openDir = QDir::homePath();
    }

    return openDir;
}

void RMainWindow::autoScale()
{
    if (currentRMatList.at(0).getInstrument() == instruments::MAG)
    {
        sliderValueHigh = convertScaleToSlider(100.0f);
        sliderValueLow = convertScaleToSlider(-100.0f);
    }
    else
    {
        sliderValueHigh = convertScaleToSlider(currentRMatList.at(0).getIntensityHigh());
        sliderValueLow = convertScaleToSlider(currentRMatList.at(0).getIntensityLow());
    }

    gamma = 1.0f;
    sliderValueGamma = convertGammaToSlider(gamma);
    // Below, the slider ignore the input if it is equal to the current state, but it prevents from going through
    // scaleImageSlot subsequently. So we have to tickle it a bit.
    ui->sliderHigh->setValue(sliderValueHigh-1);
    ui->sliderHigh->setValue(sliderValueHigh);
    ui->sliderLow->setValue(sliderValueLow +1);
    ui->sliderLow->setValue(sliderValueLow);
    ui->sliderGamma->setValue(sliderValueGamma);

    qDebug("RMainWindow::autoScale() sliderValueHigh = %i , sliderValueLow = %i", sliderValueHigh, sliderValueLow);

}

void RMainWindow::gammaScaleImageSlot(int value)
{
    if (currentROpenGLWidget == NULL || currentRMatList.empty())
        return;

    if (newMin == newMax)
        return;

    gamma = convertSliderToGamma(value);
    qDebug("RMainWindow::gammaScaleImageSlot() value =%i, gamma = %f", value, gamma);
    currentROpenGLWidget->setGamma(gamma);
    currentROpenGLWidget->update();
}

void RMainWindow::updateSliderValueSlot(double valueD)
{
    if (this->sender() == ui->doubleSpinBoxHigh)
    {
        sliderValueHigh = convertScaleToSlider((float) valueD);
        ui->sliderHigh->setValue(sliderValueHigh);
    }
    else if (this->sender() == ui->doubleSpinBoxLow)
    {
        sliderValueLow = convertScaleToSlider((float) valueD);
        ui->sliderLow->setValue(sliderValueLow);
    }
    else if (this->sender() == ui->doubleSpinBoxGamma)
    {
        sliderValueGamma = convertGammaToSlider((float) valueD);
        qDebug("RMainWindow::updateSliderValueSlot() valueD = %f , sliderValueGamma = %i", valueD, sliderValueGamma);
        ui->sliderGamma->setValue(sliderValueGamma);
    }
}

void RMainWindow::updateWB(int value)
{

    double colorFactor = (double) (value * sliderScaleWB);

    if (this->sender() == ui->redSlider)
    {
        ui->redSpinBox->setValue(colorFactor);
        if (currentROpenGLWidget == NULL)
        {
            return;
        }
        currentROpenGLWidget->setWbRed(colorFactor);
    }

    if (this->sender() == ui->greenSlider)
    {
        ui->greenSpinBox->setValue(colorFactor);
        if (currentROpenGLWidget == NULL)
        {
            return;
        }
        currentROpenGLWidget->setWbGreen(colorFactor);
    }

    if (this->sender() == ui->blueSlider)
    {
        ui->blueSpinBox->setValue(colorFactor);
        if (currentROpenGLWidget == NULL)
        {
            return;
        }
        currentROpenGLWidget->setWbBlue(colorFactor);
    }

    if (currentROpenGLWidget == NULL)
    {
        return;
    }

    currentROpenGLWidget->update();
}

int RMainWindow::convertGammaToSlider(float gamma)
{
    int sliderValueGamma = (int) std::round((gamma - gammaMin)/gammaScale + 1.0f);
    qDebug("RMainWindow::convertGammaToSlider()  gammaMin = %f, gammaScale = %f , gamma = %f ,sliderValueGamma = %i", gammaMin, gammaScale, gamma, sliderValueGamma);

    return sliderValueGamma;
}


float RMainWindow::convertSliderToScale(int value)
{

    float scaledValue = 0.0f;

    if (currentRMatList.empty())
    {
        qDebug("currentRMatList is empty.");
        return (float) scaledValue;
    }

    float dataMin = 0;

    if (currentRMatList.at(0).getDataMin() < 0)
    {
        dataMin = currentRMatList.at(0).getDataMin();
    }

    scaledValue = (float) ( dataMin + sliderScale * (value-1) );
    qDebug() << "RMainWindow::convertSliderToScale:: scaledValue =" << scaledValue;
    return scaledValue;
}


int RMainWindow::convertScaleToSlider(float valueF)
{


    int sliderValue = 1;

    if (currentRMatList.empty())
    {
        qDebug("currentRMatList is empty.");
        return sliderValue;
    }

    float dataMin = 0;
    if (currentRMatList.at(0).getDataMin() < 0)
    {
        dataMin = currentRMatList.at(0).getDataMin();
    }

    // If the slider Value is greater than the range of the slider, it is ignored instead of clipped.
    // So we need to clip to the max value of the slider which is the slider range.

    sliderValue = (int) std::min(std::round((valueF - dataMin)/sliderScale) + 1, sliderRange);

    return sliderValue;
}

float RMainWindow::convertSliderToGamma(int value)
{
    float valueF = (float) (gammaMin +  gammaScale * (value -1));
    return valueF;
}


void RMainWindow::updateCurrentData(ROpenGLWidget *rOpenGLWidget)
{
    currentROpenGLWidget = rOpenGLWidget;
    currentRMatList = currentROpenGLWidget->getRMatImageList();
    newMax = currentROpenGLWidget->getNewMax();
    newMin = currentROpenGLWidget->getNewMin();
    // Restore slider values
    // Scaling
    ui->sliderHigh->setValue(convertScaleToSlider(newMax));
    ui->sliderLow->setValue(convertScaleToSlider(newMin));
    ui->sliderGamma->setValue(convertGammaToSlider(currentROpenGLWidget->getGamma()));
    // White balance
    ui->redSlider->setValue((int) std::round(currentROpenGLWidget->getWbRed() / sliderScaleWB));
    ui->greenSlider->setValue((int) std::round(currentROpenGLWidget->getWbGreen() / sliderScaleWB));
    ui->blueSlider->setValue((int) std::round(currentROpenGLWidget->getWbBlue() / sliderScaleWB));
    // Time series
    ui->sliderFrame->setRange(0, currentRMatList.size()-1);
    ui->sliderFrame->setValue(currentROpenGLWidget->getFrameIndex());
    ui->imageLabel->setText(QString::number(currentROpenGLWidget->getFrameIndex()+1) + QString("/") + QString::number(currentRMatList.size()));

    drawHist(currentRMatList.at(0).getMatHist());


}

void RMainWindow::updateFrameInSeries(int frameIndex)
{
    // The sliderValue minimum is 1. The frameIndex minimum is 0;
    ui->imageLabel->setText(QString::number(frameIndex+1) + QString("/") + QString::number(currentRMatList.size()));

    currentROpenGLWidget->setFrameIndex(frameIndex);
    currentROpenGLWidget->initSubQImage();
    currentROpenGLWidget->updateSubQImage();
    currentROpenGLWidget->update();

    QString windowTitle = currentROpenGLWidget->getWindowTitleList().at(frameIndex);

    ui->mdiArea->activeSubWindow()->setWindowTitle(windowTitle);

}


void RMainWindow::setupExportMastersDir()
{
    // This is meant to be used when pressing the
    // exportMastersDirButton from the ui.

    QString dir = QFileDialog::getExistingDirectory(
                            0,
                            "Select directory to export Masters",
                            checkExistingDir());

    ui->mastersDirLineEdit->setText(dir);

}

void RMainWindow::setupExportCalibrateDir()
{
    QString dir = QFileDialog::getExistingDirectory(
                            0,
                            "Select directory to export Masters",
                            checkExistingDir());

    ui->calibrateDirLineEdit->setText(dir);

    emit messageSignal(QString(" Exporting calibrated files to: ") + dir);
}

void RMainWindow::stepForward()
{
    int newValue = ui->sliderFrame->value() + 1;
    if (newValue > ui->sliderFrame->maximum())
    {
        newValue = 0;
    }
    ui->sliderFrame->setValue(newValue);
    updateFrameInSeries(newValue);
}

void RMainWindow::stepBackward()
{
    int newValue = ui->sliderFrame->value() - 1;
    if (newValue < ui->sliderFrame->minimum())
    {
        newValue = ui->sliderFrame->maximum();
    }

    ui->sliderFrame->setValue(newValue);
    updateFrameInSeries(newValue);
}

void RMainWindow::playForward()
{
// Need to use QThread for using interface buttons. Otherwise intercept some keyboard key.
    stopButtonStatus = false;
    while (stopButtonStatus == false)
    {
        qApp->processEvents();
        stepForward();
        QThread::msleep(100);
    }

}

void RMainWindow::stopButtonPressed()
{
    stopButtonStatus = true;
    qDebug("Stop Button pressed");
}



void RMainWindow::tileView()
{
    ui->mdiArea->setViewMode(QMdiArea::SubWindowView);
    ui->mdiArea->tileSubWindows();
}

