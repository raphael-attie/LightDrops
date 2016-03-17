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


RProcessing::RProcessing(QObject *parent): QObject(parent)
{


}


void RProcessing::loadRMatLightList(QList<QUrl> urls)
{
    RListImageManager listImageManager(urls);
    rMatLightList = listImageManager.rMatImageList;
}

void RProcessing::loadRMatBiasList(QList<QUrl> urls)
{
    RListImageManager listImageManager(urls);
    rMatBiasList = listImageManager.rMatImageList;
}

void RProcessing::loadRMatDarkList(QList<QUrl> urls)
{
    RListImageManager listImageManager(urls);
    rMatDarkList = listImageManager.rMatImageList;
}

void RProcessing::loadRMatFlatList(QList<QUrl> urls)
{
    RListImageManager listImageManager(urls);
    rMatFlatList = listImageManager.rMatImageList;
}

void RProcessing::makeMasters()
{

    biasSuccess = makeMasterBias();
    darkSuccess = makeMasterDark();
    flatSuccess = makeMasterFlat();

}

bool RProcessing::makeMasterBias()
{
    if (!rMatBiasList.empty())
    {
        masterBias = average(rMatBiasList);
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
    if (!rMatDarkList.empty())
    {
        masterDark = average(rMatDarkList);
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
    if (!rMatFlatList.empty())
    {
        masterFlat = average(rMatFlatList);
        return true;
    }
    else
    {
        qDebug("RProcessing::rMatFlatList is empty");
        return false;
    }
}

RMat RProcessing::average(QList<RMat> rMatList)
{
    // Averages a series of cv::Mat images using arithmetic mean.
    RMat rImage(rMatList.at(0));

    if (rMatList.size() == 1)
    {
        return rImage;
    }

    int naxis1 = rImage.matImage.cols;
    int naxis2 = rImage.matImage.rows;

    // Need to create a Mat image that will host the result of the average.
    // It needs to be the same type, and has the number of channels as the Mat images of the series.
    cv::Mat avgImg = cv::Mat::zeros(naxis2, naxis1, rImage.matImage.type());

    // The accumulate function in the for-loop below can only work with up to 3 color channels
    cv::Mat tempImage;
    for(int i = 0; i < rMatList.size(); i++)
    {
        tempImage = rMatList.at(i).matImage;
        cv::accumulate(tempImage, avgImg);
    }

    avgImg = avgImg / (float) rMatList.size();
    rImage.matImage = avgImg;

    return rImage;

}


void RProcessing::exportMastersToFits()
{
    QDir mastersDir(exportMastersDir);
    if (!masterBias.empty())
    {
        masterBiasPath = mastersDir.filePath("masterBias.fits");
        QFile QFilePath(masterBiasPath);
        if (QFilePath.exists())
        {
            masterBiasPath = mastersDir.filePath("newMasterBias.fits");
        }
        exportToFits(masterBias, masterBiasPath);

    }

    if (!masterDark.empty())
    {
        masterDarkPath = mastersDir.filePath("masterDark.fits");
        QFile QFilePath(masterDarkPath);
        if (QFilePath.exists())
        {
            masterDarkPath = mastersDir.filePath("newMasterDark.fits");
        }
        exportToFits(masterDark, masterDarkPath);

    }

    if (!masterFlat.empty())
    {
        masterFlatPath = mastersDir.filePath("masterFlat.fits");
        QFile QFilePath(masterFlatPath);
        if (QFilePath.exists())
        {
            masterFlatPath = mastersDir.filePath("newMasterFlat.fits");
        }
        exportToFits(masterFlat, masterFlatPath);
    }

}


void RProcessing::exportToFits(RMat rImage, QString QStrFilename)
{

    // Write fits files
    // To do: need to add FITS keyword about bayer type.
    std::string strFilename(QStrFilename.toStdString());
    fitsfile *fptr; /* pointer to the FITS file; defined in fitsio.h */
    long fpixel = 1, naxis = 2, nPixels;
    int status = 0; /* initialize status before calling fitsio routines */
    long naxes[2];
    naxes[0] = (long) rImage.matImage.cols;
    naxes[1] = (long) rImage.matImage.rows;
    nPixels = (long) (rImage.matImage.cols * rImage.matImage.rows);
    bool bayer = rImage.isBayer();
    char keyname[] = "BAYER";

    // Create new file
    fits_create_file(&fptr, strFilename.c_str(), &status);

    // If the images is still bayer (CFA), should convert back to original DSLR precision of 16 bits unsigned to save memory.
    if (rImage.isBayer())
    {
        cv::Mat tempImage16;
        rImage.matImage.convertTo(tempImage16, CV_16U);
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)tempImage16.data, &status);
    }
    else
    {
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TFLOAT, fpixel, nPixels, (float*)rImage.matImage.data, &status);
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
    if (imageManager->rMatImage.empty())
    {
        qDebug("ProcessingWidget::setupMasterDark():: master Dark is empty");
        tempMessageSignal(QString("master Dark is empty"));
        return;
    }
    masterDark = RMat(imageManager->rMatImage);

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
    if (imageManager->rMatImage.empty())
    {
        qDebug("ProcessingWidget::setupMasterDark():: master Dark is empty");
        tempMessageSignal(QString("master Dark is empty"));
        return;
    }
    masterFlat = RMat(imageManager->rMatImage);
    masterFlat = normalizeByMean(masterFlat);

}


//---------------- public slots -----------------------//


void RProcessing::createMasters()
{
    if (exportMastersDir.isEmpty())
    {
        qDebug("ProcessingWidget::Master directory is empty");
        tempMessageSignal(QString("Master directory is empty"));
        return;
    }

    if (!treeWidget->getBiasUrls().empty())
        loadRMatBiasList(treeWidget->getBiasUrls());
    if (!treeWidget->getDarkUrls().empty())
        loadRMatDarkList(treeWidget->getDarkUrls());
    if (!treeWidget->getFlatUrls().empty())
        loadRMatFlatList(treeWidget->getFlatUrls());

    makeMasters();
    exportMastersToFits();

    if (biasSuccess)
        emit resultSignal(masterBias, QString("master Bias"));
    if (darkSuccess)
        emit resultSignal(masterDark, QString("master Dark"));
    if (flatSuccess)
        emit resultSignal(masterFlat, QString("master Flat"));
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

    QList<RMat> resultList;
    for (int i = 0 ; i < treeWidget->getLightUrls().size() ; i++)
    {
        cv::Mat tempMat;
        tempMat.create(masterDark.matImage.rows, masterDark.matImage.cols, masterDark.matImage.type());
        resultList << RMat(tempMat, masterDark.isBayer());
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
RMat RProcessing::normalizeByMean(RMat rImage)
{
    {
    //    cv::Mat tempMat;
    //    rImage.matImage.copyTo(tempMat);
        RMat rImageNorm(rImage);
        // Force a copy. We don't want to overwrite rImage.matImage with the normalized one.
        rImage.matImage.copyTo(rImageNorm.matImage);
        // Normalize image by mean value
        cv::Scalar meanValue = cv::mean(rImage.matImage);
        float meanValueF = (float) meanValue.val[0];
        qDebug() << "ProcessingWidget:: rImageNorm.matImage.channels()=" << rImageNorm.matImage.channels();
        qDebug() << "ProcessingWidget:: rImageNorm.matImage.type()=" << rImageNorm.matImage.type();
        qDebug() << "ProcessingWidget:: meanValueF=" << meanValueF;
        //cv::normalize(rImage.matImage, rImageNorm.matImage, 0.5, 1.5, cv::NORM_MINMAX, CV_32F);
        rImageNorm.matImage = rImageNorm.matImage / meanValueF;
        return rImageNorm;
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
