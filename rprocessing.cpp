#include "rprocessing.h"

#include <QFileDialog>

// cifitsio
#include <fitsio.h>

//opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/video.hpp>
#include <opencv2/core/ocl.hpp>

#include "imagemanager.h"
#include "parallelcalibration.h"


RProcessing::RProcessing(QObject *parent): QObject(parent), masterBias(NULL), masterDark(NULL), masterFlat(NULL)
{
    if(!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    listImageManager = new RListImageManager();
}

RProcessing::~RProcessing()
{
    if (masterBias != NULL)
    {
        delete masterBias;
    }

    if (!masterDark != NULL)
    {
        delete masterDark;
    }

    if (!masterFlat != NULL)
    {
        delete masterFlat;
    }
}

void RProcessing::loadRMatLightList(QList<QUrl> urls)
{
    listImageManager->loadData(urls);
    rMatLightList = listImageManager->getRMatImageList();
}

void RProcessing::loadRMatBiasList(QList<QUrl> urls)
{
    listImageManager->loadData(urls);
    rMatBiasList = listImageManager->getRMatImageList();
}

void RProcessing::loadRMatDarkList(QList<QUrl> urls)
{
    listImageManager->loadData(urls);
    rMatDarkList = listImageManager->getRMatImageList();
}

void RProcessing::loadRMatFlatList(QList<QUrl> urls)
{
    listImageManager->loadData(urls);
    rMatFlatList = listImageManager->getRMatImageList();
}



void RProcessing::exportMastersToFits()
{
    if (masterBias != NULL && !masterBias->empty())
    {
        QFileInfo fileInfo(treeWidget->getBiasDir().filePath(QString("masterBias.fits")));
        masterBiasPath = setupFileName(fileInfo, QString("fits"));

        exportToFits(masterBias, masterBiasPath);
    }

    if (masterDark != NULL && !masterDark->empty())
    {
        QFileInfo fileInfo(treeWidget->getDarkDir().filePath(QString("masterDark.fits")));
        masterDarkPath = setupFileName(fileInfo, QString("fits"));

        exportToFits(masterDark, masterDarkPath);

    }

    if (masterFlat != NULL && !masterFlat->empty())
    {

        QFileInfo fileInfo(treeWidget->getFlatDir().filePath(QString("masterFlat.fits")));
        masterFlatPath = setupFileName(fileInfo, QString("fits"));

        exportToFits(masterFlat, masterFlatPath);
    }

    tempMessageSignal(QString("Exported master calibration frames"), 0);
}

QString RProcessing::setupFileName(QFileInfo fileInfo, QString format)
{
    QString filePath = fileInfo.filePath();

    uint fileNumber = 1;
    QFileInfo fileInfoTest(filePath);
    while (fileInfoTest.exists())
    {
        QString baseName = fileInfo.baseName() + QString("_") + QString::number(fileNumber) +QString(".");
        filePath = fileInfo.absoluteDir().filePath(baseName + format);
        fileInfoTest = QFileInfo(filePath);
        fileNumber++;
        qDebug() << "RProcessing::setupFileName():: filePath =" << filePath;
    }

    return filePath;
}


void RProcessing::exportToFits(RMat *rMatImage, QString QStrFilename)
{

    // Write fits files
    // To do: need to add FITS keyword about bayer type.
    std::string strFilename(QStrFilename.toStdString());
    fitsfile *fptr; /* pointer to the FITS file; defined in fitsio.h */
    long fpixel = 1, naxis = 2, nPixels;
    int status = 0; /* initialize status before calling fitsio routines */
    long naxes[2];
    naxes[0] = (long) rMatImage->getMatImage().cols;
    naxes[1] = (long) rMatImage->getMatImage().rows;
    nPixels = (long) (rMatImage->getMatImage().cols * rMatImage->getMatImage().rows);
    int bayer = (int) rMatImage->isBayer();
    char keyname[] = "BAYER";

    qDebug() << "RProcessing::exportToFits() bayer =" << bayer;


    // Create new file
    fits_create_file(&fptr, strFilename.c_str(), &status);

    // If the images is still bayer (CFA), should convert back to original DSLR precision of 16 bits unsigned to save memory.
    if (rMatImage->isBayer())
    {
        cv::Mat tempImage16;
        rMatImage->getMatImage().convertTo(tempImage16, CV_16U);
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)tempImage16.data, &status);
    }
    else if (rMatImage->getMatImage().type() == CV_16U)
    {
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)rMatImage->getMatImage().data, &status);
    }
    else
    {
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TFLOAT, fpixel, nPixels, (float*)rMatImage->getMatImage().data, &status);
    }

    // Write header BAYER keyword
    fits_write_key(fptr, TLOGICAL, keyname, &bayer, NULL, &status);
    // Close the file
    fits_close_file(fptr, &status);
}

void RProcessing::loadMasterDark()
{
    // We can load a master Dark file only if the url exists in the treeWidget

    if (masterDarkPath.isEmpty())
    {
        qDebug("You need at lest one Dark image in the calibration tree");
        tempMessageSignal(QString("You need at lest one Dark image in the calibration tree"));
        return;
    }

    ImageManager *imageManager = new ImageManager(masterDarkPath);
    if (imageManager->getRMatImage()->empty())
    {
        qDebug("ProcessingWidget::setupMasterDark():: master Dark is empty");
        tempMessageSignal(QString("master Dark is empty"));
        return;
    }
    masterDark = imageManager->getRMatImage();

}

void RProcessing::loadMasterFlat()
{
    // We can load a master Dark file only if the url exists in the treeWidget

    if (masterFlatPath.isEmpty())
    {
        qDebug("You need at lest one Dark image in the calibration tree");
        tempMessageSignal(QString("You need at lest one Dark image in the calibration tree"));
        return;
    }

    ImageManager *imageManager = new ImageManager(masterFlatPath);
    if (imageManager->getRMatImage()->empty())
    {
        qDebug("ProcessingWidget::setupMasterDark():: master Dark is empty");
        tempMessageSignal(QString("master Dark is empty"));
        return;
    }
    masterFlat = imageManager->getRMatImage();
    masterFlat = normalizeByMean(masterFlat);

}

void RProcessing::showMinMax(const cv::Mat &matImage)
{
    double min, max;
    cv::minMaxLoc(matImage, &min, &max);
    qDebug("RProcessing::showMinMax:: min =%f , max =%f", min, max);
}


void RProcessing::createMasters()
{
    /// The Urls of the Bias, Dark, Flat images in the treeWidget are only assigned for off-screen calibration,
    /// with drag'n'drop in the treeWidget.
    /// They remain empty if the images are dropped in the central widget.

    biasSuccess = makeMasterBias();
    darkSuccess = makeMasterDark();
    flatSuccess = makeMasterFlat();

    if (biasSuccess)
    {
        emit resultSignal(masterBias);
    }

    if (darkSuccess)
    {
        emit resultSignal(masterDark);
    }

    if (flatSuccess)
    {
        emit resultSignal(masterFlat);
    }


    if (!biasSuccess & !darkSuccess & !flatSuccess)
    {
        qDebug("ProcessingWidget:: No data processed");
        tempMessageSignal(QString("No data processed"));
        return;
    }
    tempMessageSignal(QString("Calibration masters ready. You may export. "), 0);
}

bool RProcessing::makeMasterBias()
{

    if (!treeWidget->getBiasUrls().empty())
    {
        loadRMatBiasList(treeWidget->getBiasUrls());
        qDebug() << "RProcessing::makeMasterBias() rMatBiasList.at(0)->isBayer() =" << rMatBiasList.at(0)->isBayer();
        masterBias = average(rMatBiasList);
        masterBias->setImageTitle(QString("master Bias"));
        return true;
    }

    else if (!treeWidget->rMatBiasList.empty() && !treeWidget->rMatBiasList.empty())
    {
        rMatBiasList = treeWidget->rMatBiasList;
        masterBias = average(rMatBiasList);
        masterBias->setImageTitle(QString("master Bias"));
        return true;
    }
    else
    {
        qDebug("RProcessing::rMatBiasList is empty");
        return false;
    }

}

bool RProcessing::makeMasterDark()
{

    if (!treeWidget->getDarkUrls().empty())
    {
        loadRMatDarkList(treeWidget->getDarkUrls());
        masterDark = average(rMatDarkList);
        masterDark->setImageTitle(QString("master Dark"));
        return true;
    }
    else if (!treeWidget->rMatDarkList.empty() && !treeWidget->rMatDarkList.empty())
    {
        rMatDarkList = treeWidget->rMatDarkList;
        masterDark = average(rMatDarkList);
        masterDark->setImageTitle(QString("master Dark"));
        return true;
    }
    else
    {
        qDebug("RProcessing::rMatDarkList is empty");
        return false;
    }
}

bool RProcessing::makeMasterFlat()
{
    if (!treeWidget->getFlatUrls().empty())
    {
        loadRMatFlatList(treeWidget->getFlatUrls());
        masterFlat = average(rMatFlatList);
        masterFlat->setImageTitle(QString("master Flat"));
        return true;
    }
    else if (!treeWidget->rMatFlatList.empty() && !treeWidget->rMatFlatList.empty())
    {
        rMatFlatList = treeWidget->rMatFlatList;
        masterFlat = average(rMatFlatList);
        masterFlat->setImageTitle(QString("master Flat"));
        return true;
    }
    else
    {
        qDebug("RProcessing::rMatFlatList is empty");
        return false;
    }
}


RMat* RProcessing::average(QList<RMat*> rMatList)
{
    // Averages a series of cv::Mat images using arithmetic mean.

    if (rMatList.size() == 1)
    {
        return rMatList.at(0);
    }

    int naxis1 = rMatList.at(0)->getMatImage().cols;
    int naxis2 = rMatList.at(0)->getMatImage().rows;

    // Need to create a Mat image that will host the result of the average.
    // It needs to be the same type, and has the number of channels as the Mat images of the series.
    cv::Mat avgImg = cv::Mat::zeros(naxis2, naxis1, rMatList.at(0)->getMatImage().type());
    avgImg.convertTo(avgImg, CV_32F);

     // The accumulate function in the for-loop below can only work with up to 3 color channels

    for(int i = 0; i < rMatList.size(); i++)
    {
        cv::Mat tempImage;
//        cv::Mat imageClone = rMatList->operator [](i).getMatImage().clone();

//        imageClone.convertTo(tempImage, CV_32F);

        rMatList.at(i)->getMatImage().convertTo(tempImage, CV_32F);
        //avgImg = avgImg + tempImage;
        cv::accumulate(tempImage, avgImg);
    }

    avgImg = avgImg / (float) rMatList.size();

    qDebug() << "RProcessing::average() rMatList.at(0)->isBayer() =" << rMatList.at(0)->isBayer();
    RMat *rMatAvg = new RMat(avgImg, rMatList.at(0)->isBayer(), rMatList.at(0)->getInstrument());

    return rMatAvg;
}


void RProcessing::calibrateOffScreen()
{
    if (treeWidget->getLightUrls().empty())
    {
        qDebug("ProcessingWidget::calibrateOffScreen():: No lights");
        tempMessageSignal(QString("No Light image"));
        return;
    }

    if (exportCalibrateDir.isEmpty())
    {
        qDebug("ProcessingWidget::calibrateOffScreen():: No export directory for calibrated files");
        tempMessageSignal(QString("No export directory for calibrated files."));
        return;
    }

    loadMasterDark();
    loadMasterFlat();

    for (int i = 0 ; i < treeWidget->getLightUrls().size() ; i++)
    {
        cv::Mat tempMat;
        tempMat.create(masterDark->getMatImage().rows, masterDark->getMatImage().cols, masterDark->getMatImage().type());
        resultList << new RMat(tempMat, masterDark->isBayer());
    }

    cv::parallel_for_(cv::Range(0, treeWidget->getLightUrls().size()), ParallelCalibration(exportCalibrateDir, treeWidget->getLightUrls(), masterDark, masterFlat, resultList));

    qDebug("ProcessingWidget::calibrateOffScreen():: Done.");
    tempMessageSignal(QString("Calibration completed."));

    emit listResultSignal(resultList);

}


void RProcessing::setExportMastersDir(QString dir)
{
    this->exportMastersDir = dir;
    qDebug() << "Export masters to: " << exportMastersDir;
}

void RProcessing::setExportCalibrateDir(QString dir)
{
    this->exportCalibrateDir = dir;
    qDebug() << "Export calibrated data to: " << exportCalibrateDir;
}





// Private member functions
RMat* RProcessing::normalizeByMean(RMat *rMatImage)
{
    {
    //    cv::Mat tempMat;
    //    rMatImage->getMatImage().copyTo(tempMat);
        RMat *rMatImageNorm = new RMat();
        rMatImageNorm->setBayer(rMatImage->isBayer());
        rMatImageNorm->setBscale(rMatImage->getBscale());
        rMatImageNorm->setBzero(rMatImage->getBzero());
        // Force a copy. We don't want to overwrite rMatImage->getMatImage() with the normalized one.
        rMatImage->getMatImage().copyTo(rMatImageNorm->getMatImage());
        // Normalize image by mean value
        cv::Scalar meanValue = cv::mean(rMatImage->getMatImage());
        float meanValueF = (float) meanValue.val[0];
        qDebug() << "ProcessingWidget:: rMatImageNorm.getMatImage().channels()=" << rMatImageNorm->getMatImage().channels();
        qDebug() << "ProcessingWidget:: rMatImageNorm.getMatImage().type()=" << rMatImageNorm->getMatImage().type();
        qDebug() << "ProcessingWidget:: meanValueF=" << meanValueF;
        //cv::normalize(rMatImage->getMatImage(), rMatImageNorm.getMatImage(), 0.5, 1.5, cv::NORM_MINMAX, CV_32F);
        rMatImageNorm->getMatImage() = rMatImageNorm->getMatImage() / meanValueF; // Strange...
        return rMatImageNorm;
    }
}

// setters
void RProcessing::setTreeWidget(RTreeWidget *treeWidget)
{
    this->treeWidget = treeWidget;
}

// getters

QString RProcessing::getExportMastersDir()
{
    return this->exportMastersDir = exportMastersDir;
}

QString RProcessing::getExportCalibrateDir()
{
    return this->exportCalibrateDir = exportCalibrateDir;
}

RMat* RProcessing::getMasterBias()
{
    return masterBias;
}

RMat* RProcessing::getMasterDark()
{
    return masterDark;
}

RMat* RProcessing::getMasterFlat()
{
    return masterFlat;
}
