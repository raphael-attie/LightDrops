#include <QMdiSubWindow>
#include <QDropEvent>
#include <QFileDialog>

#include "rmainwindow.h"
#include "ui_rmainwindow.h"
#include "rsubwindow.h"

//QCustomPlot
#include <qcustomplot/qcustomplot.h>


RMainWindow::RMainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::RMainWindow),
    currentROpenGLWidget(NULL),
    resultROpenGLWidget(NULL),
    cannyContoursROpenGLWidget(NULL),
    limbFittingROpenGLWidget(NULL)
{
    ui->setupUi(this);
    setCentralWidget(ui->mdiArea);
    processing = new RProcessing(this);
    processing->setTreeWidget(ui->treeWidget);

    vertLineHigh = NULL;

    sliderScale = 1.0;
    gammaScale = 0.1;
    gammaMin = ui->doubleSpinBoxGamma->minimum();
    gamma = 1.0f;
    sliderScaleWB = 0.01;

    setupSubImage();

//    tabifyDockWidget(ui->currentFrameDock, ui->settingsDock);
//    this->setTabPosition(Qt::DockWidgetArea::BottomDockWidgetArea, QTabWidget::TabPosition::South);

    //this->setTabShape(QTabWidget::Triangular);

    /// Connect the change of slider value to the number displayed in the spinBox
    connect(ui->sliderHigh, SIGNAL(valueChanged(int)), this, SLOT(updateDoubleSpinBox(int)));
    connect(ui->sliderLow, SIGNAL(valueChanged(int)), this, SLOT(updateDoubleSpinBox(int)));
    connect(ui->sliderGamma, SIGNAL(valueChanged(int)), this, SLOT(updateDoubleSpinBox(int)));

    /// Connect the change of the number in the spinBox to the change of slider value
    ///Use editingFinished() instead of valueChanged(double()) to avoid unwanted feedback on the slider.
    connect(ui->doubleSpinBoxHigh, SIGNAL(editingFinished()), this, SLOT(updateSliderValueSlot()));
    connect(ui->doubleSpinBoxLow, SIGNAL(editingFinished()), this, SLOT(updateSliderValueSlot()));
    connect(ui->doubleSpinBoxGamma, SIGNAL(editingFinished()), this, SLOT(updateSliderValueSlot()));

    /// Connect the change of slider value to the image linear scaling and gamma scaling.
    connect(ui->sliderHigh, SIGNAL(valueChanged(int)), this, SLOT(scaleImageSlot(int)));
    connect(ui->sliderLow, SIGNAL(valueChanged(int)), this, SLOT(scaleImageSlot(int)));
    connect(ui->sliderGamma, SIGNAL(valueChanged(int)), this, SLOT(gammaScaleImageSlot(int)));

    /// Connect the white balance sliders
    connect(ui->redSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWB(int)));
    connect(ui->blueSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWB(int)));
    connect(ui->greenSlider, SIGNAL(valueChanged(int)), this, SLOT(updateWB(int)));

    /// Connect the pushButton and lineEdit of the ui for creating the masters in the processing class
    connect(ui->makeMasterButton, SIGNAL(pressed()), processing, SLOT(createMasters()));

    connect(ui->calibrateDirLineEdit, SIGNAL(textChanged(QString)), processing, SLOT(setExportCalibrateDir(QString)));

    /// Connect the sliderFrame and playback buttons
    connect(ui->sliderFrame, SIGNAL(valueChanged(int)), this, SLOT(updateFrameInSeries(int)));
    connect(ui->forwardButton, SIGNAL(released()), this, SLOT(stepForward()));
    connect(ui->backwardButton, SIGNAL(released()), this, SLOT(stepBackward()));
    connect(ui->playButton, SIGNAL(released()), this, SLOT(playForward()));
    connect(ui->stopButton, SIGNAL(released()), this, SLOT(stopButtonPressed()));

    /// Connect processing outputs to the main ui.
    connect(processing, SIGNAL(tempMessageSignal(QString,int)), this->statusBar(), SLOT(showMessage(QString,int)));
    connect(processing, SIGNAL(resultSignal(RMat*)), this, SLOT(createNewImage(RMat*)));
    connect(processing, SIGNAL(listResultSignal(QList<RMat*>)), this, SLOT(createNewImage(QList<RMat*>)));

    /// Connect ui radio button to tree widget to send out window titles.
    connect(ui->lightRButton, SIGNAL(clicked(bool)), this, SLOT(radioRMatSlot()));
    connect(ui->biasRButton, SIGNAL(clicked(bool)), this, SLOT(radioRMatSlot()));
    connect(ui->darkRButton, SIGNAL(clicked(bool)), this, SLOT(radioRMatSlot()));
    connect(ui->flatRButton, SIGNAL(clicked(bool)), this, SLOT(radioRMatSlot()));

    // Connect the radioLight, radioBias,... signals to the treeWidget
    connect(this, SIGNAL(radioLightImages(QList<RMat*>)), ui->treeWidget, SLOT(rMatFromLightRButton(QList<RMat*>)));
    connect(this, SIGNAL(radioBiasImages(QList<RMat*>)), ui->treeWidget, SLOT(rMatFromBiasRButton(QList<RMat*>)));
    connect(this, SIGNAL(radioDarkImages(QList<RMat*>)), ui->treeWidget, SLOT(rMatFromDarkRButton(QList<RMat*>)));
    connect(this, SIGNAL(radioFlatImages(QList<RMat*>)), ui->treeWidget, SLOT(rMatFromFlatRButton(QList<RMat*>)));

    // Connect Canny detection
    connect(ui->cannySlider, SIGNAL(sliderReleased()), this, SLOT(updateCannyDetection()));
    connect(ui->cannyRegisterButton, SIGNAL(released()), this, SLOT(cannyRegisterSeries()));

    // Connect QImage scenes
    connect(processing, SIGNAL(resultQImageSignal(QImage&)), this, SLOT(createNewImage(QImage&)));
    //connect(processing, SIGNAL(ellipseSignal(cv::RotatedRect)), this, SLOT(addEllipseToScene(cv::RotatedRect)));

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
    if (newRListImageManager->getRMatImageList().empty())
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
    currentROpenGLWidget = new ROpenGLWidget(newRListImageManager, this);

    addImage(currentROpenGLWidget);
    dispatchRMatImages(currentROpenGLWidget->getRMatImageList());
}



void RMainWindow::createNewImage(QList<RMat*> newRMatImageList)
{
    currentROpenGLWidget = new ROpenGLWidget(newRMatImageList, this);
    addImage(currentROpenGLWidget);

}

void RMainWindow::createNewImage(RMat *rMatImage)
{
    currentROpenGLWidget = new ROpenGLWidget(rMatImage, this);
    addImage(currentROpenGLWidget);
}

void RMainWindow::createNewImage(QImage &image)
{
//    cv::Mat ellMat = cv::Mat::zeros(500,500, CV_8U);
//    //color = cv::Scalar( 0, 255, 0);
//    //cv::ellipse(ellMat, ellRect, color, 2, 8);

//    QImage image2(ellMat.data, ellMat.cols, ellMat.rows, QImage::Format_Grayscale8);

    QSize currentSize = ui->mdiArea->currentSubWindow()->size();
    QPixmap pixMap = QPixmap::fromImage(image);


    qDebug("RMainWindow::createNewImage(QImage &image):: pixMap.size() = [%i ; %i]", pixMap.size().width(), pixMap.size().height());
    //QPainter *painter = new QPainter(&image);
    scene = new QGraphicsScene;
    QGraphicsPixmapItem *pixMapItem = scene->addPixmap(pixMap);

    pixMapItem->setScale(currentROpenGLWidget->getResizeFac());

    graphicsView= new QGraphicsView(scene);
    graphicsView->scale(1, -1);

    QMdiSubWindow *newSubWindow = new QMdiSubWindow;
    newSubWindow->setWidget(graphicsView);

    ui->mdiArea->addSubWindow(newSubWindow);

    newSubWindow->show();
    newSubWindow->resize(currentSize);

}

//void RMainWindow::addEllipseToScene(cv::RotatedRect rect)
//{
////    QGraphicsEllipseItem * QGraphicsScene::addEllipse(qreal x, qreal y, qreal w, qreal h, const QPen & pen = QPen(), const QBrush & brush = QBrush())
//    QPen pen;
//    pen.setColor(Qt::green);
//    QGraphicsEllipseItem *ellipseItem = new QGraphicsEllipseItem(0, 0, rect.size.width, rect.size.height);
//    ellipseItem->setPen(pen);
//    //ellipseItem->setScale(-currentROpenGLWidget->getResizeFac());
//    scene->addItem(ellipseItem);
//    qDebug("RMainWindow::addEllipseToScene():: rect.size.width = %f", rect.size.width);
//}

void RMainWindow::dispatchRMatImages(QList<RMat*> rMatList)
{
    // Try to guess the calibration type (Light, Bias, Dark, or Flat)
    // with the file names, or directory name, etc... and "radio" it out as if we were pressing the corresponding radio button
    // to assign it to the corresponding QTreeWidgetItem in the treeWidget.

    bool noBias = !(rMatList.at(0)->getImageTitle().contains(QString("bias"), Qt::CaseInsensitive) |
                   rMatList.at(0)->getImageTitle().contains(QString("offset"), Qt::CaseInsensitive) |
                   rMatList.at(0)->getFileInfo().absolutePath().contains(QString("bias"), Qt::CaseInsensitive));

    bool noDark = !(rMatList.at(0)->getImageTitle().contains(QString("dark"), Qt::CaseInsensitive) |
                   rMatList.at(0)->getFileInfo().absolutePath().contains(QString("dark"), Qt::CaseInsensitive));

    bool noFlat = !(rMatList.at(0)->getImageTitle().contains(QString("flat"), Qt::CaseInsensitive) |
                   rMatList.at(0)->getFileInfo().absolutePath().contains(QString("flat"), Qt::CaseInsensitive));

    bool isAmbiguous = rMatList.at(0)->getFileInfo().absolutePath().contains(QString("light"), Qt::CaseInsensitive);

    bool isBias = (!noBias & noDark & noFlat & !isAmbiguous) | rMatList.at(0)->getImageTitle().contains(QString("bias"), Qt::CaseInsensitive) |  rMatList.at(0)->getImageTitle().contains(QString("offset"), Qt::CaseInsensitive);
    bool isDark = (!noDark & noBias & noFlat & !isAmbiguous) | rMatList.at(0)->getImageTitle().contains(QString("dark"), Qt::CaseInsensitive) ;
    bool isFlat = (!noFlat & noDark & noBias & !isAmbiguous) | rMatList.at(0)->getImageTitle().contains(QString("flat"), Qt::CaseInsensitive) ;


    if (isBias)
    {
        ui->biasRButton->click();
    }
    else if (isDark)
    {
        ui->darkRButton->click();
    }
    else if (isFlat)
    {
        ui->flatRButton->click();
    }
    else
    {
        emit radioLightImages(currentROpenGLWidget->getRMatImageList());
    }
}

void RMainWindow::addImage(ROpenGLWidget *rOpenGLWidget)
{

    ui->sliderFrame->setRange(0, rOpenGLWidget->getRMatImageList().size()-1);
    ui->sliderFrame->setValue(0);
    ui->imageLabel->setText(QString("1") + QString("/") + QString::number(rOpenGLWidget->getRMatImageList().size()));

    currentScrollArea = new QScrollArea();
    currentScrollArea->setWidget(rOpenGLWidget);
    /// The size of the scrollArea must be set to the size of whatever is needed for the ROpenGLWidget
    /// to display the whole FOV in the central widget when adding the height of the scroll bars.
    /// Note that the "height" of the scrollBar always represent the dimension perpendicular to the direction
    /// of the slider.
    resizeScrollArea(rOpenGLWidget, currentScrollArea);
    loadSubWindow(currentScrollArea);
    rOpenGLWidget->update();

    displayPlotWidget(rOpenGLWidget);
    setupSliders(rOpenGLWidget);

    autoScale(rOpenGLWidget);

    processing->setCurrentROpenGLWidget(rOpenGLWidget);

    connect(rOpenGLWidget, SIGNAL(gotSelected(ROpenGLWidget*)), this, SLOT(changeROpenGLWidget(ROpenGLWidget*)));
    connect(rOpenGLWidget, SIGNAL(sendSubQImage(QImage*,float,int,int)), this, SLOT(updateSubFrame(QImage*,float,int,int)));

}


void RMainWindow::loadSubWindow(QScrollArea *scrollArea)
{
    currentSubWindow = new QMdiSubWindow;
    currentSubWindow->setWidget(scrollArea);
    ui->mdiArea->addSubWindow(currentSubWindow);
    currentSubWindow->setWindowTitle(currentROpenGLWidget->getRMatImageList().at(ui->sliderFrame->value())->getImageTitle());
    /// ROpenGLWidget can only be resized once it has been assigned to the QMdiSubWindow.
    /// Otherwise the shaders will not be bound.
    currentROpenGLWidget->resize(oglSize);
    qDebug() << "RMainWindow::loadSubWindow():: oglSize=" << oglSize;

    currentSubWindow->resize(scrollArea->size());
    currentSubWindow->show();

}



void RMainWindow::updateDoubleSpinBox(int value)
{
    if (this->sender() == ui->sliderHigh)
    {
        float valueF = convertSliderToScale(value);
        ui->doubleSpinBoxHigh->setValue(valueF);
        qDebug("updateDoubleSpinBox:: (high) valueF = %f", valueF);
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
    if (currentROpenGLWidget == NULL)
        return;

    if (this->sender() == ui->sliderHigh)
    {
        currentROpenGLWidget->setNewMax(convertSliderToScale(value));
    }
    else if (this->sender() == ui->sliderLow)
    {
        currentROpenGLWidget->setNewMin(convertSliderToScale(value));
    }
    else
    {
        return;
    }

    if (currentROpenGLWidget->getNewMin() == currentROpenGLWidget->getNewMax())
    {
        return;
    }

    updateCurrentROpenGLWidget();
}



void RMainWindow::setupSliders(ROpenGLWidget* rOpenGLWidget)
{
    // The slider range of integer values must be consistent with the maximum data range.
    // For USET for example, it is 4096.
    // We also need to define the number of decimals in the spinBox High and Low.
    int matType = rOpenGLWidget->getRMatImageList().at(0)->matImage.type();
    int decimals = 0;
    ui->sliderHigh->blockSignals(true);
    ui->sliderLow->blockSignals(true);
    if (rOpenGLWidget->getRMatImageList().at(0)->getInstrument() == instruments::USET)
    {
        qDebug("RMainWindow::setupSliders():: setting sliders for USET data.");        
        ui->sliderHigh->setRange(1, 4096);
        ui->sliderLow->setRange(1, 4096);
    }
    else if (matType == CV_16U || matType == CV_16UC3 || matType == CV_32F || matType == CV_32FC3)
    {
        ui->sliderHigh->setRange(1, 65536);
        ui->sliderLow->setRange(1, 65536);
    }

    else if (matType == CV_8U || matType == CV_8UC3)
    {
        ui->sliderHigh->setRange(1, 256);
        ui->sliderLow->setRange(1, 256);
    }
    else
    {
        qDebug("RMainWindow::setupSliders():: ERROR. Unknown image type. Could not setup sliders");
        exit(1);
    }


    ui->sliderHigh->blockSignals(false);
    ui->sliderLow->blockSignals(false);

    sliderRange = (float) ui->sliderHigh->maximum() - ui->sliderHigh->minimum() + 1;

    if (rOpenGLWidget->getRMatImageList().at(0)->getBscale() == 1
            && rOpenGLWidget->getRMatImageList().at(0)->getInstrument() != instruments::MAG)
    {
        sliderScale = 1.0;
        sliderToScaleMinimum = 0;
    }
    else
    {
         /// Here the image is assumed to be seen as scientific data
         /// for which scaling needs to be tightly set around the min and max
         /// so we can scan through with maximum dynamic range.
        float dataMax = (float) rOpenGLWidget->getRMatImageList().at(0)->getDataMax();
        float dataMin = (float) rOpenGLWidget->getRMatImageList().at(0)->getDataMin();
        float dataRange = dataMax - dataMin;
        sliderScale = (float) dataRange / sliderRange;
        sliderToScaleMinimum = dataMin;
        qDebug() << "RMainWindow::setupSliders() dataMax =" << dataMax;
        qDebug() << "RMainWindow::setupSliders() dataMin =" << dataMin;
        qDebug() << "RMainWindow::setupSliders() sliderScale =" << sliderScale;
        // update the number of decimals in the spinBox High and Low
        decimals = 2;
    }

    ui->doubleSpinBoxHigh->blockSignals(true);
    ui->doubleSpinBoxLow->blockSignals(true);

    ui->doubleSpinBoxHigh->setMaximum(convertSliderToScale(ui->sliderHigh->maximum()));
    ui->doubleSpinBoxHigh->setMinimum(convertSliderToScale(ui->sliderHigh->minimum()));
    ui->doubleSpinBoxHigh->setDecimals(decimals);
    ui->doubleSpinBoxLow->setMaximum(convertSliderToScale(ui->sliderLow->maximum()));
    ui->doubleSpinBoxLow->setMinimum(convertSliderToScale(ui->sliderLow->minimum()));
    ui->doubleSpinBoxLow->setDecimals(decimals);

    ui->doubleSpinBoxHigh->blockSignals(false);
    ui->doubleSpinBoxLow->blockSignals(false);
}

void RMainWindow::setupSubImage()
{
    subQImage = new QImage(ui->subFrame->getNaxis1(), ui->subFrame->getNaxis2(), QImage::Format_ARGB32);
    subQImage->fill(QColor(125, 125, 125, 125));

    ui->subFrame->setImage(subQImage);
    ui->subFrame->setDrawCross(true);
    ui->subFrame->update();

    //subPaintWidget = new PaintWidget(subQImage, subMagnification);
    //subPaintWidget->setDrawCross(true);


}


void RMainWindow::cannyEdgeDetection()
{
    if (ui->treeWidget->rMatLightList.isEmpty() && ui->treeWidget->getLightUrls().empty())
    {
        ui->statusBar->showMessage(QString("No lights for Canny edge detection"), 3000);
        return;
    }

    processing->setShowContours(ui->contoursCheckBox->isChecked());
    processing->setShowLimb(ui->limbCheckBox->isChecked());

    if (!ui->treeWidget->rMatLightList.isEmpty())
    {
        processing->cannyEdgeDetection(ui->cannySlider->value());
        /// Results need to be displayed as ROpenGLWidget as we will, de facto,
        /// deal with time series
        createNewImage(processing->getContoursRMatList());
        cannyScrollArea = currentScrollArea;
        cannySubWindow = ui->mdiArea->currentSubWindow();
        //cannySubWindow->move(currentPos);
        cannyContoursROpenGLWidget = currentROpenGLWidget;
    }
    else if (!ui->treeWidget->getLightUrls().empty())
    {
        processing->cannyEdgeDetectionOffScreen(ui->cannySlider->value());
    }



}


void RMainWindow::updateCannyDetection()
{
    if (processing->getContoursRMatList().isEmpty())
    {
        return;
    }

    processing->setShowContours(ui->contoursCheckBox->isChecked());
    processing->setShowLimb(ui->limbCheckBox->isChecked());

    processing->cannyEdgeDetection(ui->cannySlider->value());

    quint32 imageCoordX = currentROpenGLWidget->getImageCoordX();
    quint32 imageCoordY = currentROpenGLWidget->getImageCoordY();

    delete cannyContoursROpenGLWidget;
    //ROpenGLWidget *ptr = cannyContoursROpenGLWidget;
    cannyContoursROpenGLWidget = new ROpenGLWidget(processing->getContoursRMatList(), this);
    //delete ptr;


    currentROpenGLWidget = cannyContoursROpenGLWidget;

    cannyScrollArea->setWidget(cannyContoursROpenGLWidget);
    resizeScrollArea(cannyContoursROpenGLWidget, currentScrollArea);
    cannyContoursROpenGLWidget->resize(oglSize);

    cannyContoursROpenGLWidget->update();

    displayPlotWidget(cannyContoursROpenGLWidget);
    setupSliders(cannyContoursROpenGLWidget);
    autoScale(cannyContoursROpenGLWidget);

    connect(cannyContoursROpenGLWidget, SIGNAL(gotSelected(ROpenGLWidget*)), this, SLOT(changeROpenGLWidget(ROpenGLWidget*)));
    connect(cannyContoursROpenGLWidget, SIGNAL(sendSubQImage(QImage*,float,int,int)), this, SLOT(updateSubFrame(QImage*,float,int,int)));


    currentScrollArea->update();
    cannySubWindow->update();

    currentROpenGLWidget->setImageCoordX(imageCoordX);
    currentROpenGLWidget->setImageCoordY(imageCoordY);
    currentROpenGLWidget->updateSubQImage();

//    if (limbFittingROpenGLWidget != NULL)
//    {
//        limbFittingROpenGLWidget->initialize();
//        limbFittingROpenGLWidget->update();
//    }

//    /// Update the texture in the existing resultROpenGLWidget
//    resultROpenGLWidget->setRMatImageList(processing->getEllipseRMat());
//    resultROpenGLWidget->initialize();
//    resultROpenGLWidget->update();

}

void RMainWindow::cannyRegisterSeries()
{
    processing->setShowContours(ui->contoursCheckBox->isChecked());
    processing->setShowLimb(ui->limbCheckBox->isChecked());
    processing->cannyRegisterSeries();
    createNewImage(processing->getResultList());
}


void RMainWindow::displayPlotWidget(ROpenGLWidget* rOpenGLWidget)
{
    ui->histoDock->setWidget(rOpenGLWidget->fetchCurrentCustomPlot());
    rOpenGLWidget->updateCustomPlotLineItems();
}

void RMainWindow::addHeaderWidget()
{
    int frameIndex = ui->sliderFrame->value();

    if (currentROpenGLWidget == NULL)
    {
        qDebug("RMainWindow::addHeaderWidget():: currentROpenGLWidget == NULL");
        return;
    }
    else if (currentROpenGLWidget->getRListImageManager()->getUrlList().empty())
    {
      qDebug("RMainWindow::addHeaderWidget():: currentROpenGLWidget->getRListImageManager() empty");
      return;
    }
    else if (currentROpenGLWidget->getRMatImageList().at(frameIndex)->isBayer())
    {
        qDebug("RMainWindow::addHeaderWidget():: current image is Bayer type");
        return;
    }

    qDebug("Showing headerWidget");

    ui->mdiArea->addSubWindow(currentROpenGLWidget->getTableRSubWindow());
    currentROpenGLWidget->getTableRSubWindow()->show();
    currentROpenGLWidget->getTableRSubWindow()->resize(currentROpenGLWidget->getRListImageManager()->getTableWidgetList().at(0)->size());

    //Connect the slider value changed with the header dock change
    //QObject::connect(ui->sliderFrame, SIGNAL(valueChanged(int)), this, SLOT(updateTableWidget(int)));

    connect(ui->sliderFrame, SIGNAL(valueChanged(int)), currentROpenGLWidget, SLOT(setupTableWidget(int)));
    connect(currentROpenGLWidget->getTableRSubWindow(), SIGNAL(gotHidden()), this, SLOT(uncheckActionHeaderState()));
}

void RMainWindow::updateSubFrame(QImage *image, float intensity, int x, int y)
{
    char format = 'f';
    int precision = 0;

    if (std::abs(currentROpenGLWidget->getRMatImageList().at(0)->getDataMax()) <= 255)
    {
        format = 'g';
        precision = 3;
    }

    QString intensityQStr = QString::number(intensity, format, precision);

    ui->subFrame->setImage(image);
    ui->subFrame->setCursorText(intensityQStr);
    ui->subFrame->setImageCoordX(x);
    ui->subFrame->setImageCoordY(y);
    ui->subFrame->update();
}

QString RMainWindow::checkExistingDir()
{
    QString openDir = QString();
    if (!ui->calibrateDirLineEdit->text().isEmpty())
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
    if (currentROpenGLWidget == NULL)
        return;

    if (currentROpenGLWidget->getRMatImageList().at(0)->getInstrument() == instruments::MAG)
    {
        sliderValueHigh = convertScaleToSlider(100.0f);
        sliderValueLow = convertScaleToSlider(-100.0f);
    }
    else if (currentROpenGLWidget->getRMatImageList().at(0)->matImage.type() == CV_8U ||
             currentROpenGLWidget->getRMatImageList().at(0)->matImage.type() == CV_8UC3)
    {
        sliderValueHigh = convertScaleToSlider(255);
        sliderValueLow = convertScaleToSlider(0);
    }
    else
    {
        /// Make the slider use the histogram autoscaling values.
        sliderValueHigh = convertScaleToSlider(currentROpenGLWidget->getRMatImageList().at(0)->getIntensityHigh());
        sliderValueLow = convertScaleToSlider(currentROpenGLWidget->getRMatImageList().at(0)->getIntensityLow());
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

void RMainWindow::autoScale(ROpenGLWidget *rOpenGLWidget)
{
    if (rOpenGLWidget == NULL)
        return;

    if (rOpenGLWidget->getRMatImageList().at(0)->getInstrument() == instruments::MAG)
    {
        rOpenGLWidget->setNewMax(100.0f);
        rOpenGLWidget->setNewMin(-100.0f);

        sliderValueHigh = convertScaleToSlider(rOpenGLWidget->getNewMax());
        sliderValueLow = convertScaleToSlider(rOpenGLWidget->getNewMin());
    }
    else if (rOpenGLWidget->getRMatImageList().at(0)->matImage.type() == CV_8U ||
             rOpenGLWidget->getRMatImageList().at(0)->matImage.type() == CV_8UC3)
    {
        rOpenGLWidget->setNewMax(255);
        rOpenGLWidget->setNewMin(0);

        sliderValueHigh = convertScaleToSlider(rOpenGLWidget->getNewMax());
        sliderValueLow = convertScaleToSlider(rOpenGLWidget->getNewMin());

    }
    else
    {
        /// Make the slider use the histogram autoscaling values.
        rOpenGLWidget->setNewMax(rOpenGLWidget->getRMatImageList().at(0)->getIntensityHigh());
        rOpenGLWidget->setNewMin(rOpenGLWidget->getRMatImageList().at(0)->getIntensityLow());

        sliderValueHigh = convertScaleToSlider(rOpenGLWidget->getNewMax());
        sliderValueLow = convertScaleToSlider(rOpenGLWidget->getNewMin());

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

void RMainWindow::minMaxScale()
{
    if (currentROpenGLWidget == NULL)
        return;

    gamma = 1.0f;
    sliderValueGamma = convertGammaToSlider(gamma);

    sliderValueHigh = convertScaleToSlider(currentROpenGLWidget->getRMatImageList().at(0)->getDataMax());
    sliderValueLow = convertScaleToSlider(currentROpenGLWidget->getRMatImageList().at(0)->getDataMin());

    ui->sliderHigh->setValue(sliderValueHigh-1);
    ui->sliderHigh->setValue(sliderValueHigh);
    ui->sliderLow->setValue(sliderValueLow +1);
    ui->sliderLow->setValue(sliderValueLow);
    ui->sliderGamma->setValue(sliderValueGamma);

    qDebug("RMainWindow::autoScale() sliderValueHigh = %i , sliderValueLow = %i", sliderValueHigh, sliderValueLow);

}

void RMainWindow::updateCurrentROpenGLWidget()
{
    float dataRange = (float) (currentROpenGLWidget->getNewMax() - currentROpenGLWidget->getNewMin()) ;
    alpha = 1.0f / dataRange;
    beta = (float) (- currentROpenGLWidget->getNewMin() / dataRange);

    currentROpenGLWidget->setAlpha(alpha);
    currentROpenGLWidget->setBeta(beta);
    currentROpenGLWidget->updateSubQImage();
    currentROpenGLWidget->update();
    currentROpenGLWidget->updateCustomPlotLineItems();
}

void RMainWindow::gammaScaleImageSlot(int value)
{
    if (currentROpenGLWidget == NULL)
        return;

    if (currentROpenGLWidget->getNewMin() == currentROpenGLWidget->getNewMax())
        return;

    gamma = convertSliderToGamma(value);
    qDebug("RMainWindow::gammaScaleImageSlot() value =%i, gamma = %f", value, gamma);
    currentROpenGLWidget->setGamma(gamma);
    currentROpenGLWidget->update();
}

void RMainWindow::updateSliderValueSlot()
{
    if (currentROpenGLWidget == NULL)
    {
        return;
    }

    if (this->sender() == ui->doubleSpinBoxHigh)
    {
        double valueD = ui->doubleSpinBoxHigh->value();
        sliderValueHigh = convertScaleToSlider((float) valueD);
        ui->sliderHigh->setValue(sliderValueHigh);
    }
    else if (this->sender() == ui->doubleSpinBoxLow)
    {
        double valueD = ui->doubleSpinBoxLow->value();
        sliderValueLow = convertScaleToSlider((float) valueD);
        ui->sliderLow->setValue(sliderValueLow);
    }
    else if (this->sender() == ui->doubleSpinBoxGamma)
    {
        double valueD = ui->doubleSpinBoxGamma->value();
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
    currentROpenGLWidget->updateSubQImage();
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

    if (currentROpenGLWidget->getRMatImageList().empty())
    {
        qDebug("currentROpenGLWidget->getRMatImageList() is empty.");
        return (float) scaledValue;
    }

    scaledValue = sliderToScaleMinimum + sliderScale * ((float) (value-1)) ;
    qDebug("RMainWindow::convertSliderToScale:: sliderToScaleMinimum = %f ; sliderScale = %f ; value = %i ; scaledValue = %f", sliderToScaleMinimum, sliderScale, value, scaledValue);
    return scaledValue;
}


int RMainWindow::convertScaleToSlider(float valueF)
{

    int sliderValue = 1;

    if (currentROpenGLWidget->getRMatImageList().empty())
    {
        qDebug("currentROpenGLWidget->getRMatImageList() is empty.");
        return sliderValue;
    }

    // If the slider Value is greater than the range of the slider, it is ignored instead of clipped.
    // So we need to clip to the max value of the slider which is the slider range.
    sliderValue = (int) std::min(std::round((valueF - sliderToScaleMinimum)/sliderScale) + 1, sliderRange);

    return sliderValue;
}

float RMainWindow::convertSliderToGamma(int value)
{
    float valueF = (float) (gammaMin +  gammaScale * (value -1));
    return valueF;
}


void RMainWindow::changeROpenGLWidget(ROpenGLWidget *rOpenGLWidget)
{

    currentROpenGLWidget = rOpenGLWidget;

    // Restore slider values
    qDebug("RMainWindow::changeROpenGLWidget:: rOpenGLWidget->getNewMin() = %f", rOpenGLWidget->getNewMin());
    setupSliders(rOpenGLWidget);
    // Scaling
    qDebug("RMainWindow::changeROpenGLWidget:: rOpenGLWidget->getNewMin() = %f", rOpenGLWidget->getNewMin());
    ui->sliderHigh->setValue(convertScaleToSlider(currentROpenGLWidget->getNewMax()));
    ui->sliderLow->setValue(convertScaleToSlider(currentROpenGLWidget->getNewMin()));

    qDebug("RMainWindow::changeROpenGLWidget:: currentROpenGLWidget->getNewMin() = %f", currentROpenGLWidget->getNewMin());
    qDebug("RMainWindow::changeROpenGLWidget:: rOpenGLWidget->getNewMin() = %f", rOpenGLWidget->getNewMin());

    ui->sliderGamma->setValue(convertGammaToSlider(currentROpenGLWidget->getGamma()));
    // White balance
    ui->redSlider->setValue((int) std::round(currentROpenGLWidget->getWbRed() / sliderScaleWB));
    ui->greenSlider->setValue((int) std::round(currentROpenGLWidget->getWbGreen() / sliderScaleWB));
    ui->blueSlider->setValue((int) std::round(currentROpenGLWidget->getWbBlue() / sliderScaleWB));
    // Time series
    ui->sliderFrame->setRange(0, currentROpenGLWidget->getRMatImageList().size()-1);
    ui->sliderFrame->setValue(currentROpenGLWidget->getFrameIndex());
    ui->imageLabel->setText(QString::number(currentROpenGLWidget->getFrameIndex()+1) + QString("/") + QString::number(currentROpenGLWidget->getRMatImageList().size()));


    displayPlotWidget(rOpenGLWidget);

    // Refresh state of QAction buttons in toolbar
    if (currentROpenGLWidget->getTableRSubWindow() != NULL && currentROpenGLWidget->getTableRSubWindow()->isVisible())
    {
        ui->actionHeader->setChecked(true);
    }
    else
    {
        ui->actionHeader->setChecked(false);
    }

}

void RMainWindow::updateFrameInSeries(int frameIndex)
{
    /// Updates the display of the currentROpenGLWidget: change current image in the series, subImage, and histogram
    // The sliderValue minimum is 1. The frameIndex minimum is 0;
    ui->imageLabel->setText(QString::number(frameIndex+1) + QString("/") + QString::number(currentROpenGLWidget->getRMatImageList().size()));
    currentROpenGLWidget->changeFrame(frameIndex);
    ui->mdiArea->currentSubWindow()->setWindowTitle(currentROpenGLWidget->getRMatImageList().at(frameIndex)->getImageTitle());
    displayPlotWidget(currentROpenGLWidget);
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

void RMainWindow::exportMastersToFits()
{
    processing->exportMastersToFits();
}

void RMainWindow::exportFrames()
{
    QDir exportDir(ui->exportDirEdit->text());
    QString format("jpg");

    for (int i = 0 ; i < currentROpenGLWidget->getRMatImageList().size() ; i++)
    {
        QString fileName(QString("image_") + QString::number(i) + QString(".") + format);
        QFileInfo fileInfo(exportDir.filePath(fileName));
        QString filePath = processing->setupFileName(fileInfo, format);
        qDebug() << "RMainWindow::exportFrames():: filePath =" << filePath;
        QImage newImage = currentROpenGLWidget->grabFramebuffer();
        newImage.save(filePath);
        ui->sliderFrame->setValue(ui->sliderFrame->value() + 1);

    }

}

void RMainWindow::convertTo8Bit()
{
    if (currentROpenGLWidget == NULL)
    {
        return;
    }

    int idx = ui->sliderFrame->value();
    cv::Mat mat8bit = currentROpenGLWidget->getRMatImageList().at(idx)->matImage.clone();
    mat8bit.convertTo(mat8bit, CV_32F);
    mat8bit = mat8bit / currentROpenGLWidget->getRMatImageList().at(idx)->getMean();
    cv::normalize(mat8bit, mat8bit, 0, 255, cv::NORM_MINMAX);
    mat8bit.convertTo(mat8bit, CV_8U);

    RMat *rMat8bit = new RMat(mat8bit, false);
    rMat8bit->setImageTitle(QString("8-bit"));
    createNewImage(rMat8bit);
}

void RMainWindow::convertToNeg()
{
    int idx = ui->sliderFrame->value();
    cv::Mat matNeg = currentROpenGLWidget->getRMatImageList().at(idx)->matImage.clone();
    matNeg = currentROpenGLWidget->getRMatImageList().at(idx)->getDataMax() - matNeg;

    RMat *negRMat = new RMat(matNeg, false);
    negRMat->setImageTitle(QString("Negative"));
    negRMat->setInstrument(currentROpenGLWidget->getRMatImageList().at(idx)->getInstrument());
    createNewImage(negRMat);
}


void RMainWindow::calibrateOffScreenSlot()
{
    processing->calibrateOffScreen();
}

void RMainWindow::registerSeries()
{
    processing->registerSeries();
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
    int fps = ui->spinBoxFps->value();
    ulong delay = (ulong) std::round( 1000.0f / (float) fps );
    while (stopButtonStatus == false)
    {
        qApp->processEvents();
        stepForward();
        QThread::msleep(delay);
        fps = ui->spinBoxFps->value();
        delay = (ulong) std::round( 1000.0f / (float) fps );
    }

}

void RMainWindow::increaseFps()
{
    int fps = ui->spinBoxFps->value();

    if (fps < 5)
    {
        fps = 5;
    }
    else
    {
        fps += 5;
    }

    ui->spinBoxFps->setValue(fps);
}

void RMainWindow::decreaseFps()
{
    int fps = std::max(1, ui->spinBoxFps->value() -5);
    ui->spinBoxFps->setValue(fps);
}

void RMainWindow::stopButtonPressed()
{
    stopButtonStatus = true;
    qDebug("Stop Button pressed");
}


void RMainWindow::radioRMatSlot()
{
    if (currentROpenGLWidget->getRMatImageList().isEmpty())
    {
        qDebug("nothing to radio out");
        return;
    }

    if (this->sender() == ui->lightRButton)
    {
        emit radioLightImages(currentROpenGLWidget->getRMatImageList());
    }

    if (this->sender() == ui->biasRButton)
    {
        double min = 0;
        double max = 0;
        cv::minMaxLoc(currentROpenGLWidget->getRMatImageList().at(0)->matImage, &min, &max);
        qDebug("RMainWindow::radioRMatSlot():: min =%f , max =%f", min, max );
        emit radioBiasImages(currentROpenGLWidget->getRMatImageList());
    }

    if (this->sender() == ui->darkRButton)
    {
        emit radioDarkImages(currentROpenGLWidget->getRMatImageList());
    }

    if (this->sender() == ui->flatRButton)
    {
        emit radioFlatImages(currentROpenGLWidget->getRMatImageList());
    }
}



void RMainWindow::tileView()
{
    ui->mdiArea->setViewMode(QMdiArea::SubWindowView);
    ui->mdiArea->tileSubWindows();
}


void RMainWindow::on_actionHeader_triggered()
{

}

void RMainWindow::on_actionHeader_toggled(bool arg1)
{
    if (arg1)
    {
        currentROpenGLWidget->setupTableWidget(ui->sliderFrame->value());
        addHeaderWidget();
    }
    else
    {
        if (currentROpenGLWidget->getTableRSubWindow() != NULL && currentROpenGLWidget->getTableRSubWindow()->isVisible())
        {
            currentROpenGLWidget->getTableRSubWindow()->hide();
        }
    }
}

void RMainWindow::uncheckActionHeaderState()
{
    ui->actionHeader->setChecked(false);
}

void RMainWindow::resizeScrollArea(ROpenGLWidget *rOpenGLWidget, QScrollArea *scrollArea)
{

    /// Need to leave just enough space for the scroll bars. So we have to remove the scrollBarHeight from the
    /// the size of ROpenGLWidget to make sure the latter fits tightly in the QScrollArea.
    int naxis1 = rOpenGLWidget->getRMatImageList().at(0)->matImage.cols;
    int naxis2 = rOpenGLWidget->getRMatImageList().at(0)->matImage.rows;

    float widthFac = (float) ((this->centralWidget()->width() - scrollArea->horizontalScrollBar()->height()) / (float) naxis1);
    float heightFac = (float) ((this->centralWidget()->height() - scrollArea->horizontalScrollBar()->height()) / (float) naxis2);
    float resizeFac = 1;

    if (widthFac > 1 && heightFac > 1)
    {
        //then default size = image size
        oglSize = QSize(naxis1, naxis2);
    }
    else // set other default size, smaller than original image size
    {
        resizeFac = std::min(widthFac, heightFac);
        oglSize = QSize((int) ((float) naxis1 * resizeFac), (int) ((float) naxis2 * resizeFac));
    }

    rOpenGLWidget->setResizeFac(resizeFac);
    scrollArea->resize(oglSize + QSize(5, 25));

}

void RMainWindow::on_actionTileView_triggered()
{
    ui->mdiArea->tileSubWindows();
}
