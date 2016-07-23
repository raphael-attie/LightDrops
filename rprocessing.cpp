#include "rprocessing.h"

#include <QFileDialog>

// cifitsio
#include <cfitsio/fitsio.h>

//opencv
#include <opencv2/world.hpp>
// Arrayfire
#include <arrayfire.h>

#include "imagemanager.h"
#include "parallelcalibration.h"
#include "typedefs.h"


RProcessing::RProcessing(QObject *parent): QObject(parent),
    masterBias(NULL), masterDark(NULL), masterFlat(NULL), masterFlatN(NULL), useROI(false), useXCorr(false),
    radius(0), radius1(0), radius2(0), radius3(0), masterWithMean(true), masterWithSigmaClip(false)
{
    listImageManager = new RListImageManager();
}

RProcessing::~RProcessing()
{
    if (masterBias != NULL)
    {
        delete masterBias;
    }

    /// masterDark sometimes will share the same pointer
    /// as masterBias so we need to check that we don't
    /// delete another time what is already deleted.
    if (masterDark != NULL && masterDark != masterBias)
    {
        delete masterDark;
    }

    if (masterFlat != NULL)
    {
        delete masterFlat;
    }

    if (masterFlatN != NULL)
    {
        delete masterFlatN;
    }

//    if (!resultList.empty())
//    {
//        qDeleteAll(resultList);
//        resultList.clear();
//    }

    if (!imageManagerList.isEmpty())
    {
        qDeleteAll(imageManagerList);
        resultList.clear();
    }

/// Do not to delete the rMatLightList if it comes from the treeWidget->rMatLightList
//    if (!rMatLightList.isEmpty())
//    {
//        qDeleteAll(rMatLightList);
//        rMatLightList.clear();
//    }

//    delete listImageManager;
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
    if (masterBias != NULL && !masterBias->matImage.empty())
    {
        QFileInfo fileInfo(treeWidget->getBiasDir().filePath(QString("masterBias.fits")));
        masterBiasPath = setupFileName(fileInfo, QString("fits"));

        exportToFits(masterBias, masterBiasPath);
    }

    if (masterDark != NULL && !masterDark->matImage.empty())
    {
        QFileInfo fileInfo(treeWidget->getDarkDir().filePath(QString("masterDark.fits")));
        masterDarkPath = setupFileName(fileInfo, QString("fits"));

        exportToFits(masterDark, masterDarkPath);

    }

    if (masterFlat != NULL && !masterFlat->matImage.empty())
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
    naxes[0] = (long) rMatImage->matImage.cols;
    naxes[1] = (long) rMatImage->matImage.rows;
    nPixels = (long) (rMatImage->matImage.cols * rMatImage->matImage.rows);
    int bayer = (int) rMatImage->isBayer();
    char keyname[] = "BAYER";

    // Create new file
    fits_create_file(&fptr, strFilename.c_str(), &status);

    // If the images is still bayer (CFA), should convert back to original DSLR precision of 16 bits unsigned to save memory.
    if (rMatImage->isBayer())
    {
        cv::Mat tempImage16;
        rMatImage->matImage.convertTo(tempImage16, CV_16U);
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)tempImage16.data, &status);
    }
    else if (rMatImage->matImage.type() == CV_16U)
    {
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)rMatImage->matImage.data, &status);
    }
    else
    {
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TFLOAT, fpixel, nPixels, (float*)rMatImage->matImage.data, &status);
    }

    // Write header BAYER keyword
    fits_write_key(fptr, TLOGICAL, keyname, &bayer, NULL, &status);

    if (rMatImage->getInstrument() == instruments::USET)
    {
        char keyNameTelescop[] = "TELESCOP";
        char keyValueTelescop[] = "USET";
        fits_write_key(fptr, TSTRING, keyNameTelescop, keyValueTelescop, NULL, &status);
    }
    else
    {
        char keyNameTelescop[] = "TELESCOP";
        char keyValueTelescop[] = "UNDEFINED";
        fits_write_key(fptr, TSTRING, keyNameTelescop, keyValueTelescop, NULL, &status);
    }


    // Close the file
    fits_close_file(fptr, &status);
}

void RProcessing::loadMasterDark()
{
    /// Used off-screen because no image is loaded beforehand in the ROpenGLWidget.
    /// So the url must exist in the treeWidget
    if (treeWidget->getDarkUrls().size() == 1)
    {
        masterDarkPath = treeWidget->getDarkUrls().at(0).toLocalFile();
    }
    else
    {

        qDebug("You need at lest one Dark image in the calibration tree");
        tempMessageSignal(QString("You need at lest one Dark image in the calibration tree"));
        return;
    }

    /// Let the ImageManager on the stack, (so we don't have to call delete).
    ///  and make a deep copy of the data in RMat before leaving this function.
    ImageManager imageManager(masterDarkPath);

    masterDark = new RMat(*imageManager.getRMatImage());

}

void RProcessing::loadMasterFlat()
{
    /// Used off-screen because no image is loaded beforehand in the ROpenGLWidget.
    /// So the url must exist in the treeWidget
    if (treeWidget->getFlatUrls().size() == 1)
    {
        masterFlatPath = treeWidget->getFlatUrls().at(0).toLocalFile();
    }
    else
    {
        qDebug("You need at lest one Flat image in the calibration tree");
        tempMessageSignal(QString("You need at lest one Flat image in the calibration tree"));
        return;
    }

    /// Let the ImageManager on the stack, (so we don't have to call delete).
    ///  and make a deep copy of the data in RMat before leaving this function.
    ImageManager imageManager(masterFlatPath);

    masterFlat = new RMat(*imageManager.getRMatImage());

    normalizeFlat();
    qDebug("RProcessing::loadMasterFlat():: min / max of normalized flat:");
    showMinMax(masterFlatN->matImage);
    //emit resultSignal(masterFlatN);
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
    {   /// Images not yet in memory,
        /// load the bias files from their urls found in the treeWidget
        loadRMatBiasList(treeWidget->getBiasUrls());

        if (masterWithMean)
        {
            masterBias = average(rMatBiasList);
            masterBias->setImageTitle( QString("master Bias (arithmetic mean)") );
        }
        else if (masterWithSigmaClip)
        {
            masterBias = sigmaClipAverage(rMatBiasList);
            masterBias->setImageTitle( QString("master Bias (sigma-clipped average)") );
        }
        return !masterBias->matImage.empty();
    }

    else if (treeWidget->getBiasUrls().empty() && !treeWidget->rMatBiasList.empty())
    {   /// Images already in memory,
        /// get the pointer to the biases (list of Mat) from the treeWidget
        rMatBiasList = treeWidget->rMatBiasList;

        if (masterWithMean)
        {
            masterBias = average(rMatBiasList);
            masterBias->setImageTitle( QString("master Bias (arithmetic mean)") );
        }
        else if (masterWithSigmaClip)
        {
            masterBias = sigmaClipAverage(rMatBiasList);
            masterBias->setImageTitle( QString("master Bias (sigma-clipped average)") );
        }

        return !masterBias->matImage.empty();
    }
    else
    {
        emit messageSignal(QString("Bias images not found"));
        return false;
    }

}

bool RProcessing::makeMasterDark()
{

    if (!treeWidget->getDarkUrls().empty())
    {   /// Images not yet in memory,
        /// load the dark files from their urls found in the treeWidget
        loadRMatDarkList(treeWidget->getDarkUrls());
        if (masterWithMean)
        {
            masterDark = average(rMatDarkList);
            masterDark->setImageTitle( QString("master Dark (arithmetic mean)") );
        }
        else if (masterWithSigmaClip)
        {
            masterDark = sigmaClipAverage(rMatDarkList);
            masterDark->setImageTitle( QString("master Dark (sigma-clipped average)") );
        }
        masterDark->setImageTitle(QString("master Dark"));
        return !masterDark->matImage.empty();
    }
    else if (treeWidget->getDarkUrls().empty() && !treeWidget->rMatDarkList.empty())
    {   /// Images already in memory,
        /// get the pointer to the darks (list of Mat) from the treeWidget
        rMatDarkList = treeWidget->rMatDarkList;

        if (masterWithMean)
        {
            masterDark = average(rMatDarkList);
            masterDark->setImageTitle( QString("master Dark (arithmetic mean)") );
        }
        else if (masterWithSigmaClip)
        {
            masterDark = sigmaClipAverage(rMatDarkList);
            masterDark->setImageTitle( QString("master Dark (sigma-clipped average)") );
        }
        return !masterDark->matImage.empty();
    }
    else
    {
        emit messageSignal(QString("Dark images not found"));
        return false;
    }
}

bool RProcessing::makeMasterFlat()
{
    if (!treeWidget->getFlatUrls().empty())
    {   /// Images not yet in memory,
        /// load the flat files from their urls found in the treeWidget
        loadRMatFlatList(treeWidget->getFlatUrls());

        if (masterWithMean)
        {
            masterFlat = average(rMatFlatList);
            masterFlat->setImageTitle( QString("master Flat (arithmetic mean)") );
        }
        else if (masterWithSigmaClip)
        {
            masterFlat = sigmaClipAverage(rMatFlatList);
            masterFlat->setImageTitle( QString("master Flat (sigma-clipped average)") );
        }

        if (masterBias != NULL)
        {
            cv::Mat tempMatFlat;
            cv::Mat tempMatBias;

            masterBias->matImage.convertTo(tempMatFlat, CV_32F);
            masterFlat->matImage.convertTo(tempMatBias, CV_32F);

            cv::subtract(tempMatFlat, tempMatBias, masterFlat->matImage);
        }
    }
    else if (treeWidget->getFlatUrls().empty() && !treeWidget->rMatFlatList.empty())
    {   /// Images already in memory,
        /// get the pointer to the flats (list of Mat) from the treeWidget
        rMatFlatList = treeWidget->rMatFlatList;

        if (masterWithMean)
        {
            masterFlat = average(rMatFlatList);
            masterFlat->setImageTitle( QString("master Flat (arithmetic mean)") );
        }
        else if (masterWithSigmaClip)
        {
            masterFlat = sigmaClipAverage(rMatFlatList);
            masterFlat->setImageTitle( QString("master Flat (sigma-clipped average)") );
        }
    }
    else
    {
        emit messageSignal(QString("Flat images not found"));
        return false;
    }

    if (masterFlat->matImage.empty())
    {
        return false;
    }

    /// At this stage, the masterFlat is available.
    /// Need to subtract the Bias if there's one.
    if (masterBias == NULL)
    {
        return true;
    }

    cv::Mat tempMatFlat;
    cv::Mat tempMatBias;
    masterBias->matImage.convertTo(tempMatFlat, CV_32F);
    masterFlat->matImage.convertTo(tempMatBias, CV_32F);
    cv::subtract(tempMatFlat, tempMatBias, masterFlat->matImage);
    //masterFlat->matImage.convertTo(masterFlat->matImage, rMatFlatList.at(0)->matImage.type());
    masterFlat->setInstrument(rMatFlatList.at(0)->getInstrument());
    /// Update masterFlat statistics

    masterFlat->calcStats();
    /// Convert masterFlat to original type
    return true;
}


RMat* RProcessing::average(QList<RMat*> rMatList)
{
    /// Averages a series of cv::Mat images using arithmetic mean.
    int naxis1 = rMatList.at(0)->matImage.cols;
    int naxis2 = rMatList.at(0)->matImage.rows;

    /// Need to create a Mat image that will host the result of the average.
    /// It needs to be the same type, and has the number of channels as the Mat images of the series.
    cv::Mat avgImg = cv::Mat::zeros(naxis2, naxis1, CV_32F);

    for(int i = 0; i < rMatList.size(); i++)
    {
        cv::Mat tempImage;
        rMatList.at(i)->matImage.convertTo(tempImage, CV_32F);
        /// Sum the images with cv::accumulate
        /// This function can only work with up to 3 color channels
        cv::accumulate(tempImage, avgImg);
    }

    avgImg = avgImg / (float) rMatList.size();
    RMat *rMatAvg = new RMat(avgImg, rMatList.at(0)->isBayer(), rMatList.at(0)->getInstrument());
    return rMatAvg;
}

RMat *RProcessing::sigmaClipAverage(QList<RMat*> rMatImageList)
{
    if (rMatImageList.size() == 1)
    {
        qDebug() << "Only 1 image. Returning input Mat Image.";
        return rMatImageList.at(0);
    }
    /// Get sigma coefficient from sigmaEdit
    float sigmaFactor = 1.0f;
    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();
    emit messageSignal(QString("Stacking %1 frames with sigma-clipping").arg(nFrames));

    /// Initialize the GPU array series
    af::array arfSeries(naxis2, naxis1, nFrames);
    /// Store the matImageList (converted to float) in the GPU array
    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        //cv::cvtColor(rMatImageList.at(k).matImage, tempMat, CV_RGB2GRAY);
        af::array tempArf(naxis2, naxis1, (float*) tempMat.data);
        //arfSeries(af::span, af::span, k) = transpose(tempArf);
        arfSeries(af::span, af::span, k) = tempArf;
    }

    af::timer start2 = af::timer::start();

    af::array meanArf = af::median(arfSeries, 2);
    qDebug("RProcessing::sigmaClipAverage::  elapsed seconds: %f us", af::timer::stop(start2));
    af::array meanArfTiled = af::tile(meanArf, 1, 1, nFrames);
    //af::array stdevArf = af::moddims(af::stdev(arfSeries, 2), naxis2, naxis1);
    af::array stdevArf = af::stdev(arfSeries, 2);
    af::array stdevArfTiled = af::tile(stdevArf, 1, 1, nFrames);
    af::array arfMaskReject = af::abs(arfSeries - meanArfTiled) > stdevArfTiled;
    arfSeries(arfMaskReject) = meanArfTiled(arfMaskReject);
    meanArf = af::mean(arfSeries, 2);
    qDebug("RProcessing::sigmaClipAverage::  elapsed seconds: %f us", af::timer::stop(start2));
    /// Copy an array from the device to the host:
    float *hostArf = meanArf.host<float>();
    qDebug("RProcessing::sigmaClipAverage::  elapsed seconds: %f us", af::timer::stop(start2));
    /// Prepare output Mat image
    cv::Mat matImage(naxis2, naxis1, CV_32F, hostArf);
    RMat *rMatAvg = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());

    return rMatAvg;
}

double RProcessing::pixelDistance(double u, double v)
{
    return cv::sqrt(u*u + v*v);
}

cv::Mat RProcessing::make2DGaussian(int matSize, double sigma)
{
    /// This assumes usage on square 2D Fourier transform
    ///
//    [x,y]=meshgrid((-nc/2+0.5):(nc/2-0.5),(-nr/2+0.5):(nr/2-0.5));
//    sigma=1/sigma;

//    r2=(2*x/nc).^2+(2*y/nr).^2;   %  Mapping the pixels in k-space

//    g=exp(-r2/(2*sigma^2)); % Gaussian Window

//    g=g/max(g(:));  % Normalized

    cv::Mat gauss2d = cv::Mat::zeros(matSize, matSize, CV_32F);
    double origin = matSize/2.0;

    double pxDist;
    for (int i = 0; i < matSize; i++)
    {
        for (int j = 0; j < matSize; j++)
        {
            pxDist = pixelDistance(i - origin, j - origin);
            gauss2d.at<float>(i, j) = exp(-pxDist / (2.0 * sigma * sigma));
        }
    }

    double min = 0;
    double max = 0;
    cv::minMaxLoc(gauss2d, &min, &max);
    gauss2d = gauss2d / max;

//    cv::Mat kernelX = cv::getGaussianKernel(matSize, sigma, CV_32F);
//    cv::Mat kernelY = cv::getGaussianKernel(matSize, sigma, CV_32F);
//    cv::Mat kernelXY = kernelX * kernelY.t();

    return gauss2d;
}

cv::Mat RProcessing::fftshift(cv::Mat matFourier)
{
    cv::Mat matShifted = matFourier.clone();

    // rearrange the quadrants of Fourier image  so that the origin is at the image center
    int cx = matShifted.cols/2;
    int cy = matShifted.rows/2;

    cv::Mat q0(matShifted, cv::Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
    cv::Mat q1(matShifted, cv::Rect(cx, 0, cx, cy));  // Top-Right
    cv::Mat q2(matShifted, cv::Rect(0, cy, cx, cy));  // Bottom-Left
    cv::Mat q3(matShifted, cv::Rect(cx, cy, cx, cy)); // Bottom-Right

    cv::Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);

    q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
    q2.copyTo(q1);
    tmp.copyTo(q2);

    return matShifted;
}

cv::Mat RProcessing::makePowerSpectrumFFT(cv::Mat matImage)
{
    cv::Mat planes[] = {cv::Mat_<float>(matImage), cv::Mat::zeros(matImage.size(), CV_32F)};
    cv::Mat complexI;
    cv::merge(planes, 2, complexI);         // Add to the expanded another plane with zeros

    cv::dft(complexI, complexI);
    // compute the magnitude and switch to logarithmic scale
    // => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
    cv::split(complexI, planes);  // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
    cv::magnitude(planes[0], planes[1], planes[0]);// planes[0] = magnitude
    cv::Mat powerSpectrum = planes[0];

    return powerSpectrum;
}

cv::Mat RProcessing::makeImageHPF(cv::Mat matImage, double sigma)
{
    /// Apply high-pass filter on given image.
    cv::Mat matImage32 = matImage.clone();
    matImage32.convertTo(matImage32, CV_32F);
    cv::Mat matImageFFT;
    cv::dft(matImage32, matImageFFT, cv::DFT_COMPLEX_OUTPUT);
    cv::Mat matImageFFTshifted = fftshift(matImageFFT);

    cv::Mat gaussianMat = make2DGaussian(matImage.rows, sigma);
    cv::Mat invGaussian = 1.0 - gaussianMat;

    vector<cv::Mat> fftMasked_ReIm;
    cv::split(matImageFFTshifted, fftMasked_ReIm);
    fftMasked_ReIm[0] = fftMasked_ReIm[0].mul(invGaussian);
    fftMasked_ReIm[1] = fftMasked_ReIm[1].mul(invGaussian);
    cv::Mat fftMasked;
    cv::merge(fftMasked_ReIm, fftMasked);

    fftMasked = fftshift(fftMasked);

    cv::Mat matImageHPF;
    cv::idft(fftMasked, matImageHPF, cv::DFT_SCALE| cv::DFT_REAL_OUTPUT );
    double min, max;
    cv::minMaxLoc(matImageHPF, &min, &max);
    matImageHPF = matImageHPF - min;

    fftMasked_ReIm.clear();

    return matImageHPF;
}


void RProcessing::calibrateOffScreen()
{
    if (treeWidget->getLightUrls().empty())
    {
        qDebug("ProcessingWidget::calibrateOffScreen():: No lights");
        tempMessageSignal(QString("No Light image"));
        return;
    }

//    if (exportCalibrateDir.isEmpty())
//    {
//        qDebug("ProcessingWidget::calibrateOffScreen()::export directory empty");
//        tempMessageSignal(QString("No export directory"));
//        return;
//    }

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    if (masterDark == NULL)
    {
        loadMasterDark();
    }

    // It's possible that after the above, masterDark is still NULL.
    if (masterDark == NULL)
    {
        if (masterBias !=NULL)
        {
            // We may need a proper copy constructor for RMat.
            masterDark = new RMat(*masterBias);
        }
    }

    if (masterFlat == NULL)
    {
        loadMasterFlat();
    }

//    for (int i = 0 ; i < treeWidget->getLightUrls().size() ; i++)
//    {
//        cv::Mat tempMat;
//        tempMat.create(masterDark->matImage.rows, masterDark->matImage.cols, masterDark->matImage.type());
//        resultList << new RMat(tempMat, masterDark->isBayer());
//    }

    //cv::parallel_for_(cv::Range(0, treeWidget->getLightUrls().size()), ParallelCalibration(treeWidget->getLightUrls(), masterDark, masterFlat, resultList));
    calibrate();
    qDebug("ProcessingWidget::calibrateOffScreen():: Done.");
    tempMessageSignal(QString("Calibration completed."), 0);

    rMatLightList = resultList; // Add to tree widget?
    resultList.clear();

    emit listResultSignal(rMatLightList);


}

void RProcessing::calibrate()
{
    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    for(int i = 0; i < treeWidget->getLightUrls().size(); i++)
    {
        QString filePathQStr = treeWidget->getLightUrls().at(i).toLocalFile();
        ImageManager lightManager(filePathQStr);
        cv::Mat matResult;
        cv::Mat lightMat = lightManager.getRMatImage()->matImage;

        if (lightMat.type() != CV_32F)
        {
            lightMat.convertTo(lightMat, CV_32F);
        }
        if (masterDark->matImage.type() != CV_32F)
        {
            masterDark->matImage.convertTo(masterDark->matImage, CV_32F);
        }

        cv::subtract(lightMat, masterDark->matImage, matResult);

        if (!masterFlatN->matImage.empty())
        {
            cv::divide(matResult, masterFlatN->matImage, matResult);
        }

        cv::threshold(matResult, matResult, 0, 0, cv::THRESH_TOZERO);

        /// Copy the result into the list.
        /// This may be an overkill. We could output directly into the elements of the list
        resultList << new RMat(matResult, masterDark->isBayer());
        resultList.at(i)->calcMinMax();
        resultList.at(i)->calcStats();
        resultList.at(i)->setBscale(masterFlatN->getDataMax()/masterFlat->getDataMax());
        resultList.at(i)->setImageTitle(QString("Limb-registered image # %i").arg(i+1));
    }
    /// We can assign the result as the rMatLightList since the processing is off-screen,
    /// and thus we do not overwrite anything.

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


void RProcessing::registerSeries()
{
    /// Here we process the rMatLightList. It is assigned in two ways:
    /// 1) By Drag and Drop in the QMdiArea
    /// 2) After a calibration like in calibrate()
    /// We use resultList as the (temporary?) output list.
    if (!treeWidget->rMatLightList.isEmpty())
    {
        rMatLightList = treeWidget->rMatLightList;
    }
    else
    {
        emit tempMessageSignal(QString("No lights to register"));
        return;
    }

    int nFrames = rMatLightList.size();

    // Define the motion model
    const int warp_mode_1 = cv::MOTION_TRANSLATION;
    const int warp_mode_2 = cv::MOTION_EUCLIDEAN;

    // Specify the number of iterations.
    int number_of_iterations_1 = 50; // stellar
    int number_of_iterations_2 = 50; // stellar
//    int number_of_iterations_1 = 100; // solar
//    int number_of_iterations_2 = 200; // solar

    // Specify the threshold of the increment
    // in the correlation coefficient between two iterations
    double termination_eps_1 = 1e-1; // stellar
    double termination_eps_2 = 1e-2; // stellar
//    double termination_eps_1 = 1e-2; // solar
//    double termination_eps_2 = 1e-3; // solar

    // Define termination criteria
    //cv::TermCriteria criteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, number_of_iterations, termination_eps);
    //cv::TermCriteria criteria(cv::TermCriteria::EPS, number_of_iterations, termination_eps);
    cv::TermCriteria criteria1(cv::TermCriteria::COUNT, number_of_iterations_1, termination_eps_1);
    cv::TermCriteria criteria2(cv::TermCriteria::EPS, number_of_iterations_2, termination_eps_2);


    // Get the 1st image of the rMatLightList as the reference image (and put as 1st element of resultList)

    cv::Mat refMat = rMatLightList.at(0)->matImage;
    refMat.convertTo(refMat, CV_32F);
    // Normalize to a multiple of the exposure time * median?
    float expTime = rMatLightList.at(0)->getExpTime();
    float mean = rMatLightList.at(0)->getMean();
    float normFactor = 1.0f / (mean);

    qDebug("RProcessing::registerSeries() image 1/%i", nFrames);
    qDebug("expTime = %f", expTime);
    qDebug("mean = %f", mean);
    qDebug("normFactor = %f", normFactor);

    refMat = refMat * normFactor;
    cv::Mat refMat2;
    cv::threshold(refMat, refMat2, 0.5, 0.5, cv::THRESH_TRUNC);

    resultList << new RMat(refMat, false);
    resultList.at(0)->setBscale(normFactor);
    // Get a resampled version. 1/4 on each axis.
    cv::Mat refMat2R;
    cv::resize(refMat2, refMat2R, cv::Size(), 0.25, 0.25, CV_INTER_AREA);


    for (int i = 1 ; i < nFrames; ++i)
    {
        qDebug("Registering image #%i/%i", i, nFrames);

        // Normalize to a multiple of the exposure time
        expTime = rMatLightList.at(i)->getExpTime();
        mean = rMatLightList.at(i)->getMean();
        normFactor = 1.0f / (mean);
        qDebug("expTime = %f", expTime);
        qDebug("mean = %f", mean);
        qDebug("normFactor = %f", normFactor);

        cv::Mat registeredMat = rMatLightList.at(i)->matImage.clone();
        registeredMat.convertTo(registeredMat, CV_32F);
        registeredMat = registeredMat * normFactor;
        cv::Mat registeredMat2;
        cv::threshold(registeredMat, registeredMat2, 0.6, 0.6, cv::THRESH_TRUNC);
        // Get a resampled version. 1/4 on each axis.
        cv::Mat registeredMat2R;
        cv::resize(registeredMat2, registeredMat2R, cv::Size(), 0.25, 0.25, CV_INTER_AREA);

        cv::Mat warp_matrix_1 = cv::Mat::eye(2, 3, CV_32F);
        double eccEps = 0;
        // 1st pass of the ECC algorithm on the decimated image. The results are stored in warp_matrix.
        eccEps = cv::findTransformECC(
                    refMat2R,
                    registeredMat2R,
                    warp_matrix_1,
                    warp_mode_1,
                    criteria1
                    );
        qDebug() << "eccEps 1 =" << eccEps / 0.25;
        warp_matrix_1.at<float>(0, 2) /= 0.25;
        warp_matrix_1.at<float>(1, 2) /= 0.25;
        std::cout << "result warp_matrix 1 =" << std::endl << warp_matrix_1 << std::endl << std::endl;

        // 2nd pass on the full resolution images.
        eccEps = cv::findTransformECC(
                    refMat2,
                    registeredMat2,
                    warp_matrix_1,
                    warp_mode_2,
                    criteria2
                    );

        qDebug() << "eccEps 2 =" << eccEps;
        std::cout << "result warp_matrix 2 =" << std::endl << warp_matrix_1 << std::endl << std::endl;

        cv::warpAffine(registeredMat, registeredMat, warp_matrix_1, registeredMat.size(), cv::INTER_NEAREST + CV_WARP_INVERSE_MAP);
        resultList << new RMat(registeredMat, false); // This RMat is necessarily non-bayer.
        resultList.at(i)->setBscale(normFactor);

    }

}

void RProcessing::registerSeriesOnLimbFit()
{
    if (limbFitWarpMat.empty())
    {
        emit messageSignal(QString("No results from limb fitting."));
        return;
    }
    // This overload is used [optionally] on top of the limb-fitting algorithm to improve the co-alignment.
    int nFrames = limbFitResultList1.size();

    if (!limbFitResultList2.isEmpty())
    {
        limbFitResultList2.clear();
    }

    // Define the motion model
    const int warp_mode_1 = cv::MOTION_TRANSLATION;
    const int warp_mode_2 = cv::MOTION_EUCLIDEAN;

    // Specify the number of iterations.
//    int number_of_iterations_1 = 50; // stellar
//    int number_of_iterations_2 = 50; // stellar
    int number_of_iterations_1 = 100; // solar
    int number_of_iterations_2 = 200; // solar

    // Specify the threshold of the increment
    // in the correlation coefficient between two iterations
//    double termination_eps_1 = 1e-1; // stellar
//    double termination_eps_2 = 1e-2; // stellar
    double termination_eps_1 = 1e-2; // solar
    double termination_eps_2 = 1e-3; // solar

    // Define termination criteria
    //cv::TermCriteria criteria(cv::TermCriteria::COUNT + cv::TermCriteria::EPS, number_of_iterations, termination_eps);
    //cv::TermCriteria criteria(cv::TermCriteria::EPS, number_of_iterations, termination_eps);
    cv::TermCriteria criteria1(cv::TermCriteria::COUNT, number_of_iterations_1, termination_eps_1);
    cv::TermCriteria criteria2(cv::TermCriteria::EPS, number_of_iterations_2, termination_eps_2);


    // Get the 1st image of the rMatLightList as the reference image (and put as 1st element of resultList)

    cv::Mat refMat;
    limbFitResultList1.at(0)->matImage.convertTo(refMat, CV_16U);


    RMat *refRMat = new RMat(refMat, false, limbFitResultList1.at(0)->getInstrument());
    limbFitResultList2 << refRMat;
    limbFitResultList2.at(0)->setImageTitle(QString("X-corr registered image # 1"));

    cv::Mat refMat0 = normalizeByThresh(limbFitResultList1.at(0)->matImage, limbFitResultList1.at(0)->getIntensityLow(), limbFitResultList1.at(0)->getIntensityHigh(), limbFitResultList1.at(0)->getNormalizeRange());
//    emit resultSignal(refMat0, false, rMatList.at(0)->getInstrument(), QString("BLAH 0"));

    // Set low and high threshold values to properly saturate the disk.
    float lowThresh = 100.0f;
    float highThresh = 2000.0f;
    // Get and normalize/saturate the new reference image (at i = 0).
    // Need to take into account the new range within which the image has been normalized, and clip it.
    cv::Mat refMat1 = normalizeClipByThresh(refMat0, lowThresh, highThresh, limbFitResultList1.at(0)->getDataRange());
//    emit resultSignal(refMat1, false, instruments::generic, QString("BLAH 1"));

    // Get ROI if any
    cv::Mat refMat2;
    if (useROI)
    {
        refMat2 = refMat1(cvRectROI);
    }
    else
    {
        refMat2 = refMat1;
    }
    refMat2.convertTo(refMat2, CV_32F);

    // Get a resampled version of the reference image. 1/4 on each axis.
    cv::Mat refMat2R;
    cv::resize(refMat2, refMat2R, cv::Size(), 0.25, 0.25, CV_INTER_AREA);


    for (int i = 1 ; i < nFrames; ++i)
    {

        qDebug("Registering image #%i/%i", i, nFrames);

        cv::Mat tempMat = normalizeByThresh(limbFitResultList1.at(i)->matImage, limbFitResultList1.at(i)->getIntensityLow(), limbFitResultList1.at(i)->getIntensityHigh(), limbFitResultList1.at(i)->getNormalizeRange());

        // Get the image to co-align with respect to the reference image
        cv::Mat registeredMat = normalizeClipByThresh(tempMat, lowThresh, highThresh, limbFitResultList1.at(i)->getDataRange());

        cv::Mat registeredMat2;
        if (useROI)
        {
            registeredMat2 = registeredMat(cvRectROI);
        }
        else
        {
            registeredMat2 = registeredMat;
        }


        registeredMat2.convertTo(registeredMat2, CV_32F);

        //resultSignal(registeredMat2, false, instruments::generic);

        // Get a resampled version. 1/4 on each axis.
        cv::Mat registeredMat2R;
        cv::resize(registeredMat2, registeredMat2R, cv::Size(), 0.25, 0.25, CV_INTER_AREA);

        cv::Mat warp_matrix_1 = cv::Mat::eye(2, 3, CV_32F);
        double eccEps = 0;
        // 1st pass of the ECC algorithm on the decimated image. The results are stored in warp_matrix.

        eccEps = cv::findTransformECC(
                    refMat2R,
                    registeredMat2R,
                    warp_matrix_1,
                    warp_mode_1,
                    criteria1
                    );

        qDebug() << "eccEps 1 =" << eccEps / 0.25;
        warp_matrix_1.at<float>(0, 2) /= 0.25;
        warp_matrix_1.at<float>(1, 2) /= 0.25;
        std::cout << "result warp_matrix 1 =" << std::endl << warp_matrix_1 << std::endl << std::endl;

        // 2nd pass on the full resolution images.
        eccEps = cv::findTransformECC(
                    refMat2,
                    registeredMat2,
                    warp_matrix_1,
                    warp_mode_1,
                    criteria2
                    );

        qDebug() << "eccEps 2 =" << eccEps;
        std::cout << "result warp_matrix 2 =" << std::endl << warp_matrix_1 << std::endl << std::endl;

        //cv::warpAffine(rMatList.at(i)->matImage, registeredMat, warp_matrix_1, registeredMat.size(), cv::INTER_CUBIC + CV_WARP_INVERSE_MAP);
        //cv::Mat warp_matrix_total = limbFitWarpMat + warp_matrix_1 ;
        cv::Mat warp_matrix_total = cv::Mat::eye( 2, 3, CV_32FC1 );
        warp_matrix_total.at<float>(0, 2) = limbFitWarpMat.at<float>(0, 2) + warp_matrix_1.at<float>(0, 2);
        warp_matrix_total.at<float>(1, 2) = limbFitWarpMat.at<float>(1, 2) + warp_matrix_1.at<float>(1, 2);
        std::cout << "warp_matrix_total = " << std::endl << " " << warp_matrix_total << std::endl << std::endl;

        cv::warpAffine(rMatLightList.at(i)->matImage, registeredMat, warp_matrix_total, registeredMat.size(),cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);

        registeredMat.convertTo(registeredMat, CV_16U);
        limbFitResultList2 << new RMat(registeredMat, false, rMatLightList.at(i)->getInstrument()); // This RMat is necessarily non-bayer.
        limbFitResultList2.at(i)->setImageTitle(QString("X-corr registered image # %1").arg(i));
    }

}

void RProcessing::registerSeriesByPhaseCorrelation()
{
    if (!treeWidget->rMatLightList.isEmpty())
    {
        rMatLightList = treeWidget->rMatLightList;
    }
    else
    {
        emit tempMessageSignal(QString("No lights to register"));
        return;
    }

    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    // Normalize the data to have equivalent statistics / histograms
    cv::Mat refMat = normalizeByThresh(rMatLightList.at(0)->matImage, rMatLightList.at(0)->getIntensityLow(), rMatLightList.at(0)->getIntensityHigh(), rMatLightList.at(0)->getNormalizeRange());
    refMat.convertTo(refMat, CV_32F);

    double sigma = 80.0;
    cv::Mat refMatHPF = makeImageHPF(refMat, sigma);

    RMat *refRMat = new RMat(*rMatLightList.at(0));
    refRMat->setImageTitle(QString("Phase-registered image # ") + QString::number(0));
    refRMat->setDate_time(rMatLightList.at(0)->getDate_time());
    resultList << refRMat;

    for (int i=1; i < rMatLightList.size(); i++)
    {
        cv::Mat matImage = normalizeByThresh(rMatLightList.at(i)->matImage, rMatLightList.at(i)->getIntensityLow(), rMatLightList.at(i)->getIntensityHigh(), rMatLightList.at(i)->getNormalizeRange());
        matImage.convertTo(matImage, CV_32F);
        cv::Mat matImageHPF = makeImageHPF(matImage, sigma);

        cv::Point2d shift = cv::phaseCorrelate(refMatHPF, matImageHPF);
        std::cout << "Shifts = " << shift << std::endl;

        cv::Mat shiftedMatImage;
        cv::Mat warpMat = cv::Mat::eye(2, 3, CV_32F);
        warpMat.at<float>(0, 2) = shift.x;
        warpMat.at<float>(1, 2) = shift.y;

        cv::warpAffine(rMatLightList.at(i)->matImage, shiftedMatImage, warpMat, matImage.size(), cv::INTER_LANCZOS4 + cv::WARP_INVERSE_MAP);

        shiftedMatImage.convertTo(shiftedMatImage, rMatLightList.at(i)->matImage.type());

        RMat *resultMat = new RMat(shiftedMatImage, false, rMatLightList.at(i)->getInstrument());
        resultMat->setImageTitle(QString("Phase-registered image # ") + QString::number(i));
        resultMat->setDate_time(rMatLightList.at(i)->getDate_time());
        resultList << resultMat;
    }


}

void RProcessing::cannyEdgeDetectionOffScreen(int thresh)
{
    if (treeWidget->getLightUrls().empty())
    {
        qDebug("ProcessingWidget::cannyEdgeDetectionOffScreen():: No lights");
        tempMessageSignal(QString("No Light image(s)"));
        return;
    }

    if (!contoursRMatList.isEmpty())
    {
        qDeleteAll(contoursRMatList);
        contoursRMatList.clear();
    }

    if (!centers.isEmpty())
    {
        centers.clear();
    }

    centers.reserve(treeWidget->getLightUrls().size());
    radius = 0;

    for(int i = 0; i < treeWidget->getLightUrls().size(); i++)
    {
        qDebug("RProcessing:: cannyEdgeDetection() on image #%i", i);
        QString filePathQStr = treeWidget->getLightUrls().at(i).toLocalFile();
        ImageManager *newImageManager = new ImageManager(filePathQStr);
        imageManagerList << newImageManager;
        rMatLightList << newImageManager->getRMatImage();
        setupCannyDetection(i);
        cannyDetect(thresh);


        limbFit(i);

        /// Get results showing contours of all the edges
        contoursRMat = new RMat(contoursMat.clone(), false);
        contoursRMat->setImageTitle(QString("All canny edges contours"));
        contoursRMatList << contoursRMat;
    }

    tempMessageSignal(QString("Canny detection done."), 0);
}



bool RProcessing::cannyEdgeDetection(int thresh)
{
    // Check if data exist
    if (rMatLightList.isEmpty())
    {
        emit tempMessageSignal(QString("No lights for Canny edge detection"));
        return false;
    }


    if (!contoursRMatList.isEmpty())
    {
        qDeleteAll(contoursRMatList);
        contoursRMatList.clear();
    }

    if (!centers.isEmpty())
    {
        centers.clear();
    }

    if (!circleOutList.isEmpty())
    {
        circleOutList.clear();
    }

    radius = 0;
    instruments instrument = rMatLightList.at(0)->getInstrument();

    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        qDebug("RProcessing:: cannyEdgeDetection() on image # %i", i+1);
        setupCannyDetection(i);
        cannyDetect(thresh);
        if (instrument == instruments::USET)
        {
            fixUset(contoursMat);
        }

        bool success = limbFit(i);
        qDebug("RProcessing:: success on image # %i", (int) success);
        if (!success)
        {
            qDebug("RProcessing:: returning.");
            return false;
        }

        /// Get results showing contours of all the edges
        contoursRMat = new RMat(contoursMat.clone(), false);
        contoursRMat->setImageTitle(QString("Canny edges: Image # %1").arg(i+1));
        contoursRMatList << contoursRMat;

    }

    radius = radius / (float) rMatLightList.size();
    radius1 = radius1 / (float) rMatLightList.size();
    radius2 = radius2 / (float) rMatLightList.size();
    qDebug("RProcessing:: average radius (fitEllipse) = %f", radius);
    qDebug("RProcessing:: average radius (HyperEllipse) = %f", radius1);
    qDebug("RProcessing:: average radius (L-M) = %f", radius2);
    // 917.05 px from cv::fitEllipse
    // 916.86 px from Hyper
    // 916.89 px from L-M

    return true;
}




void RProcessing::setupCannyDetection(int i)
{

    cv::Mat normalizedMat;

    if (useHPF)
    {
        cv::Mat matImageHPF = makeImageHPF(rMatLightList.at(i)->matImage, hpfSigma);
        RMat *rMatHPF = new RMat(matImageHPF, false, instruments::generic);
        normalizedMat = normalizeByThresh(matImageHPF, rMatHPF->getIntensityLow(), rMatHPF->getIntensityHigh(), rMatHPF->getNormalizeRange());
        normalizedMat.convertTo(sampleMatN, CV_8U, 256.0f / rMatHPF->getNormalizeRange());
        delete rMatHPF;
    }
    else
    {
        // Normalize the data to have equivalent statistics / histograms
        normalizedMat = normalizeByThresh(rMatLightList.at(i)->matImage, rMatLightList.at(i)->getIntensityLow(), rMatLightList.at(i)->getIntensityHigh(), rMatLightList.at(i)->getNormalizeRange());
        /// convert to 32-bit
        normalizedMat.convertTo(sampleMatN, CV_32F, 256.0f / rMatLightList.at(i)->getNormalizeRange());

       /// Setting for contrast stretching to boost contrast between limb and off-limb
       /// We also saturate/clip above newMin and newMax to homogenize the disk and remove unwanted features
       float newMin = 10;
       float newMax = 150;

       float newDataRange = newMax - newMin;
       float alpha = 256.0f / newDataRange;
       float beta = -newMin * 256.0f /newDataRange;

       // Convert to 8 bit with contrast stretching to boost contrast between limb and off-limb
       sampleMatN.convertTo(sampleMatN, CV_8U, alpha, beta);

       cv::threshold(sampleMatN, sampleMatN, newMax, newMax, cv::THRESH_TRUNC);
   //    sampleMatN = newMax - sampleMatN;
   //    cv::threshold(sampleMatN, sampleMatN, newMax - newMin, newMax - newMin, cv::THRESH_TRUNC);
   //    sampleMatN = newMax - sampleMatN;

    }

    //emit resultSignal(sampleMatN, false, instruments::generic, QString("Normalized for Canny"));
}

void RProcessing::cannyDetect(int thresh)
{
    /// Detect edges using cannycompareContourAreas
    double thresh1 = ((double) thresh) / 2.0;
    double thresh2 = (double) thresh;

    // Canny edge detection works best when blurring a bit.
    //cv::blur(sampleMatN, sampleMatN, cv::Size(3,3));
    //cv::blur(sampleMatN, sampleMatN, cv::Size(9,9));
    //cv::blur(sampleMatN, sampleMatN, cv::Size(13,13));
    //cv::blur(sampleMatN, sampleMatN, cv::Size(15,15));
    //cv::blur(sampleMatN, sampleMatN, cv::Size(17,17));

    cv::blur(sampleMatN, sampleMatN, cv::Size(blurSigma,blurSigma));
    cv::Canny(sampleMatN, contoursMat, thresh1, thresh2, 3);

/// void adaptiveThreshold(InputArray src, OutputArray dst, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C)
    //cv::adaptiveThreshold(sampleMatN, contoursMat, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 5, -5);

}

bool RProcessing::limbFit(int i)
{
    // Find contours
    vector< vector <cv::Point> > contours;
    vector< cv::Vec4i > hierarchy;
    cv::findContours(contoursMat.clone(), contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE, cv::Point(0, 0));

    // For comparing with Werner limb fit
    cv::Mat matImage = rMatLightList.at(i)->matImage.clone(); //normalizeByThresh(rMatLightList.at(i)->matImage, rMatLightList.at(i)->getIntensityLow(), rMatLightList.at(i)->getIntensityHigh(), rMatLightList.at(i)->getNormalizeRange());
    matImage.convertTo(contoursMat, CV_8U, 256.0f /  rMatLightList.at(0)->getNormalizeRange());

    cv::cvtColor(contoursMat, contoursMat, CV_GRAY2RGB);

    if (contours.size() == 0)
    {
        qDebug("No contours found at image %i", i+1);
        tempMessageSignal(QString("No contours found at image %1").arg(i+1));
        return false;
    }

    // sort contours
    std::sort(contours.begin(), contours.end(), compareContourAreas);
    // grab contours
    //vector< cv::Point > biggestContour = contours[contours.size()-1];
    //std::vector<cv::Point> smallestContour = contours[0];

    // gather points of all contours in one big vector

//    vector< vector<cv::Point> > selectedContours;
//    size_t nSelectedContours = std::min(contours.size(), (size_t) 200);
//    for (size_t ii = 1 ; ii <= nSelectedContours ; ++ii)
//    {
//        selectedContours.push_back(contours[contours.size() - ii]);
//    }

//    for (size_t ii = 1 ; ii < contours.size() ; ++ii)
//    {
//        selectedContours.push_back(contours[contours.size() - ii]);
//    }

    vector< vector<cv::Point> > selectedContours = contours;

    selectedContoursList << selectedContours;

    size_t nContourPoints = 0;

       for (int ii = 0; ii < selectedContours.size(); ++ii)
       {
            nContourPoints += selectedContours[ii].size();
       }

    vector< cv::Point > contours1D;
    //contours1D.reserve(nContourPoints);

    for (int ii = 0; ii < selectedContours.size(); ++ii)
      {
        const vector< cv::Point > & v = selectedContours[ii];
        contours1D.insert( contours1D.end() , v.begin() , v.end() );
        qDebug()<< ii << ": contours1D[%d] = "<< contours1D[ii].x << contours1D[ii].y;
      }

    if (showContours)
    {

                cv::Scalar color = cv::Scalar( 0, 255, 0);
//                cv::drawContours( contoursMat, selectedContours, 0, color, 2, 8);

                for( int ii = 0; ii < selectedContours.size(); ++ii )
                {
                    //cv::drawContours( contoursMat, selectedContours, ii, color, 2, 8, hierarchy, 0, cv::Point() );
                    for (int jj = 0; jj < selectedContours[ii].size(); jj++)
                    {
                        cv::Point limbPoint = selectedContours[ii][jj];
                        cv::line(contoursMat, limbPoint, limbPoint, color, 2, 8, 0);
                    }

                }

//        cv::RNG rng(12345);

//        for( int i = 0; i< contours.size(); i++ )
//        {
//            cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) ); // random rgb value
//            //cv::Scalar color = cv::Scalar( 0, 255, 0);
//            cv::drawContours( contoursMat, contours, i, color, 2, 8, hierarchy, 0, cv::Point() );
//        }
    }


    cv::RotatedRect ellRect = cv::fitEllipse(contours1D);

    ellRectList << ellRect;

    radius += (ellRect.size.height + ellRect.size.width) / 4.0f;
//    qDebug("RADIUS 1 = %f", ellRect.size.width/2.0f);
//    qDebug("RADIUS 2 = %f", ellRect.size.height/2.0f);
    qDebug("RADIUS R = %f", (ellRect.size.height + ellRect.size.width) / 4.0f);


/// ------------------------ FIT CIRCLE ---------------------------------

    reals *X = new reals[nContourPoints];
    reals *Y = new reals[nContourPoints];

    for (size_t ii = 0; ii < contours1D.size(); ++ii)
    {
        X[ii] = contours1D[ii].x;
        Y[ii] = contours1D[ii].y;
    }

    Data contourData((int) contours1D.size(), X, Y);


    Circle circleOut1 = CircleFitByHyper(contourData);
    radius1 += circleOut1.r;

//    reals circleX = ellRect.center.x;
//    reals circleY = ellRect.center.y;
//    reals circleR = 913.0f;
//    Circle circleInit(circleX, circleY, circleR);

    Circle circleInit = circleOut1;
//    circleInit.r = 913.0f;

    reals lambdaIni = 0.001;

//    qDebug("");
//    qDebug(" ----------- Circle fitting initial parameters:  ----------");
//    qDebug("");
//    qDebug("nContourPoints = %i", nContourPoints);
//    qDebug("Center at (%f ; %f)", circleInit.a, circleInit.b);
//    qDebug("");

    CircleFitByLevenbergMarquardtFull(contourData, circleInit, lambdaIni, circleOut);
    radius2 += circleOut.r;
    centers.append(cv::Point2f((float) circleOut.a, (float) circleOut.b));
    circleOutList << circleOut;

//    qDebug("LMA  output:");
//    qDebug("status = %i", status);
//    qDebug("Center at [%f ; %f]", circleOut2.a, circleOut2.b);
//    qDebug("Radius R = %f", circleOut2.r);
//    qDebug("");

//    CircleFitByLevenbergMarquardtReduced(contourData, circleInit, lambdaIni, circleOut);
//    radius2 += circleOut.r;
//    centers.append(cv::Point2f((float) circleOut.a, (float) circleOut.b));
//    circleOutList << circleOut;
/// ---------------------------------------------------------------------------- ///

    if (showLimb)
    {
//        cv::RotatedRect cvRect;
//        cvRect.center = centers.at(i);
//        cvRect.size = cv::Size2f(radius2, radius2);
//        cvRect.angle = 0.0f;

        cv::Scalar red = cv::Scalar(255, 0, 0);
//        cv::Scalar green = cv::Scalar(0, 255, 0);

        //cv::ellipse(contoursMat, ellRect, red, 2, 8);

        // draw the fitted circle
        cv::Point2f circleCenter(circleOut.a, circleOut.b);
        cv::circle(contoursMat, circleCenter, circleOut.r, red, 2, 8);
    }

    return true;
}


void RProcessing::cannyRegisterSeries()
{
    if (centers.isEmpty())
    {
        tempMessageSignal(QString("Run limb fitting first."));
        return;
    }

    if (!limbFitResultList1.isEmpty())
    {
        qDeleteAll(limbFitResultList1);
        limbFitResultList1.clear();
    }


    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        qDebug("Canny-registering image # %i ", i+1);
        /// Register series
        limbFitWarpMat = cv::Mat::eye( 2, 3, CV_32FC1 );
        cv::Point2f origin(rMatLightList.at(i)->matImage.cols / 2.0f, rMatLightList.at(i)->matImage.rows / 2.0f);

        cv::Point2f delta = centers.at(i) - origin;
        limbFitWarpMat.at<float>(0, 2) = delta.x;
        limbFitWarpMat.at<float>(1, 2) = delta.y;
        std::cout << "limbFitWarpMat = " << std::endl << " " << limbFitWarpMat << std::endl << std::endl;

        cv::Mat registeredMat = rMatLightList.at(i)->matImage.clone();
        registeredMat.convertTo(registeredMat, CV_32F);

        cv::warpAffine(registeredMat, registeredMat, limbFitWarpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);

        registeredMat.convertTo(registeredMat, CV_16U);

        RMat *resultMat = new RMat(registeredMat, false, rMatLightList.at(i)->getInstrument());
        resultMat->setImageTitle(QString("Registered image # ") + QString::number(i));
        resultMat->setDate_time(rMatLightList.at(i)->getDate_time());
        limbFitResultList1 << resultMat;
    }

    limbFitPlot = new QCustomPlot();
    limbFitPlot->addGraph();

    // Prepare the plot data
    QVector<double> frameNumbers(rMatLightList.size());
    QVector<double> radius(rMatLightList.size());
    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        frameNumbers << i;
        radius << circleOutList.at(i).r;
    }

    limbFitPlot->graph(0)->setData(frameNumbers, radius);
    limbFitPlot->rescaleAxes();
    limbFitPlot->xAxis->setRange(0, rMatLightList.size());
}

bool RProcessing::wernerLimbFit()
{
    // Check if data exist
    if (rMatLightList.isEmpty())
    {
        emit tempMessageSignal(QString("No images"));
        return false;
    }

    if (!contoursRMatList.isEmpty())
    {
        qDeleteAll(contoursRMatList);
        contoursRMatList.clear();
    }

    if (!circleOutList.isEmpty())
    {
        circleOutList.clear();
    }

//    int i = 0;
    for (int i = 0; i < rMatLightList.size(); i++)
    {

        cv::Mat matImage0 = rMatLightList.at(i)->matImage.clone(); //normalizeByThresh(rMatLightList.at(i)->matImage, rMatLightList.at(i)->getIntensityLow(), rMatLightList.at(i)->getIntensityHigh(), rMatLightList.at(i)->getNormalizeRange());
        cv::Mat matImage;
        //matImage.create(matImage0.rows, matImage0.cols, CV_32S);
        cv::Mat matImageHPF;
        matImageHPF.create(matImage0.rows, matImage0.cols, CV_32F);
        if (useHPF)
        {
            matImageHPF = makeImageHPF(matImage0, hpfSigma);
            qDebug("RPRocessing::wernerLimbFit  hpfSigma = %f", hpfSigma);
            matImageHPF.convertTo(matImage, CV_32F);
        }
        else
        {
            matImage0.convertTo(matImage, CV_16U);
        }

        int numDots = 64;

        Data wernerPoints = Data(numDots*4);
        raphFindLimb2( matImage, &wernerPoints, numDots);

        circleOut = CircleFitByTaubin(wernerPoints);

        //    Circle circleOut1 = CircleFitByHyper(wernerPoints);
        //    Circle circleInit = circleOut1;
        //    reals lambdaIni = 0.001;
        //    CircleFitByLevenbergMarquardtFull(wernerPoints, circleInit, lambdaIni, circleOut);

        /// 2nd pass of fitting:
        /// Here, the circle might still be off because of outliers (clouds, ...)
        /// Try sigma-clipping on the set of detected points.

        /// 2) Identify the outlyers
        /// 3) Reject them to define a cleaner set of points.
        /// 4) Fit this cleaner set of points (2nd pass)

        cv::Point2f circleCenter(circleOut.a, circleOut.b);
        std::vector<float> distances(wernerPoints.n);
        for (int j = 0; j < wernerPoints.n; j++)
        {
            cv::Point2f point( (float) wernerPoints.X[j], (float) wernerPoints.Y[j]);
            float distance = cv::norm(point-circleCenter);
            distances[j] = distance;

        }
        /// documentation: void meanStdDev(InputArray src, OutputArray mean, OutputArray stddev, InputArray mask=noArray())
        cv::Scalar mean, stddev;
        cv::meanStdDev(distances, mean, stddev);

        std::vector<cv::Point> newPoints;
        for (int j = 0; j < wernerPoints.n; j++)
        {
            if ( abs(distances[j] - mean[0] ) <  stddev[0] )
            {
                cv::Point2f point( (float) wernerPoints.X[j], (float) wernerPoints.Y[j]);
                newPoints.push_back(point);

            }
        }

        Data cleanDataPoints(newPoints.size());
        Circle circleOut2;
        if (!newPoints.empty())
        {
            for (int j = 0; j < newPoints.size() ; j++)
            {
                cleanDataPoints.X[j] = newPoints.at(j).x;
                cleanDataPoints.Y[j] = newPoints.at(j).y;
            }

            /// 2nd pass at fitting
            circleOut2 = CircleFitByTaubin(cleanDataPoints);
        }
        /// End of 2nd pass
        if (!newPoints.empty())
        {   /// Store the 2nd pass
            circleOutList << circleOut2;
        }
        else
        {   /// if there were only 1 pass, store it.
            circleOutList << circleOut;
        }


        /// Display results
        matImage0.convertTo(matImage, CV_16U);
        //emit resultSignal(matImage, false, instruments::USET, QString("HPF"));

        contoursMat = cv::Mat::zeros(matImage.rows, matImage.cols, CV_8U);
        matImage0.convertTo(contoursMat, CV_8U, 256.0f /  rMatLightList.at(i)->getNormalizeRange());

        cv::cvtColor(contoursMat, contoursMat, CV_GRAY2RGB);

        /// Display points in contoursMat
        for (int j = 0; j < wernerPoints.n; j++)
        {
            /// 1st pass in green
            cv::Point point( wernerPoints.X[j], wernerPoints.Y[j]);
            cv::line(contoursMat, point, point, cv::Scalar(0, 255, 0), 8, 8, 0);
        }
        if (!newPoints.empty())
        {
            for (int j = 0; j < cleanDataPoints.n; j++)
            {
                /// 2nd pass in red
                cv::Point point( cleanDataPoints.X[j], cleanDataPoints.Y[j]);
                cv::line(contoursMat, point, point, cv::Scalar(255, 0, 0), 6, 6, 0);
            }
        }


        if (showLimb)
        {

            // draw the 1st fitted circle in red
            cv::Scalar green = cv::Scalar(0, 255, 0);
            cv::Scalar red = cv::Scalar(255, 0, 0);
            cv::circle(contoursMat, circleCenter, circleOut.r, green, 2, 8);

            if (!newPoints.empty())
            {
                // draw the 2nd fitted circle in green
                cv::Point2f circleCenter2(circleOut2.a, circleOut2.b);
                cv::circle(contoursMat, circleCenter2, circleOut2.r, red, 2, 8);
            }

        }

        /// Pack the results showing contours of all the edges
        contoursRMat = new RMat(contoursMat.clone(), false);
        contoursRMat->setImageTitle(QString("werner Limb Detection: Image # %1").arg(i+1));
        contoursRMatList << contoursRMat;

    }

        return true;
}

void RProcessing::raphFindLimb(cv::Mat matImage, Data *dat, int numDots)
{
    cv::Mat blurMat, matSlice1, matSlice2, matSlice3, matSlice4;
    cv::blur(matImage, blurMat, cv::Size(7, 7));
    cv::Mat gradXLeft, gradXRight, gradYBottom, gradYTop;
    cv::Point minLoc, maxLoc;
    double minVal, maxVal;

    int naxis1 = matImage.cols;
    int naxis2 = matImage.rows;

    /// Kernels for image gradients.
    /// Over x-axis, forward (from left edge) and backward (from right-edge) direction
    cv::Mat kernelXLeft = (cv::Mat_<float>(1,3)<<-0.5, 0, 0.5);
    cv::Mat kernelXRight = (cv::Mat_<float>(1,3)<<0.5, 0, -0.5);
    /// Over y-axis, forward and backward direction
    cv::Mat kernelYBottom = (cv::Mat_<float>(3,1)<<-0.5, 0, 0.5);
    cv::Mat kernelYTop = (cv::Mat_<float>(3,1)<<0.5, 0, -0.5);

    cv::filter2D(blurMat, gradXLeft, -1, kernelXLeft, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
    cv::filter2D(blurMat, gradXRight, -1, kernelXRight, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
    cv::filter2D(blurMat, gradYBottom, -1, kernelYBottom, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
    cv::filter2D(blurMat, gradYTop, -1, kernelYTop, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);

    gradXLeft(cv::Range::all(), cv::Range(0,15)) = 0;
    gradYBottom(cv::Range(0,15), cv::Range::all()) = 0;

    for (int ii = 0; ii < numDots; ii++)
    {
        int X = naxis1/4 + ii*naxis1/(2*numDots);
        int Y = naxis2/4 + ii*naxis2/(2*numDots);
        /// Left-hand slices (no copy)
        matSlice1 = gradXLeft(cv::Range(Y, Y+1), cv::Range(0, naxis1/4));
        cv::minMaxLoc(matSlice1, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii] = maxLoc.x;
        dat->Y[ii] = Y;
        /// Right-hand slices (no copy)
        matSlice2 = gradXRight(cv::Range(Y, Y+1), cv::Range(3*naxis1/4, naxis1));
        cv::minMaxLoc(matSlice2, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + numDots] = maxLoc.x + 3*naxis1/4;
        dat->Y[ii + numDots] = Y;
        ///Bottom slices (no copy)
        matSlice3 = gradYBottom(cv::Range(0, naxis2/4), cv::Range(X, X+1));
        cv::minMaxLoc(matSlice3, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + 2*numDots] = X;
        dat->Y[ii + 2*numDots] = maxLoc.y;
        ///Top slices (no copy)
        matSlice4 = gradYTop(cv::Range(3*naxis2/4, naxis2), cv::Range(X, X+1));
        cv::minMaxLoc(matSlice4, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + 3*numDots] = X;
        dat->Y[ii + 3*numDots] = maxLoc.y + 3*naxis2/4;

    }
}

void RProcessing::raphFindLimb2(cv::Mat matImage, Data *dat, int numDots)
{
    cv::Mat matSlice1, matSlice2;
    cv::Mat gradXLeft, gradXRight, gradYBottom, gradYTop;
    cv::Point minLoc, maxLoc;
    double minVal, maxVal;
    int naxis1 = matImage.cols;
    int naxis2 = matImage.rows;
    /// Kernels for image gradients.
    /// Over x-axis, forward (from left edge) and backward (from right-edge) direction
    cv::Mat kernelXLeft = (cv::Mat_<float>(1,3)<<-0.5, 0, 0.5);
    cv::Mat kernelXRight = (cv::Mat_<float>(1,3)<<0.5, 0, -0.5);
    /// Over y-axis, forward and backward direction
    cv::Mat kernelYBottom = (cv::Mat_<float>(3,1)<<-0.5, 0, 0.5);
    cv::Mat kernelYTop = (cv::Mat_<float>(3,1)<<0.5, 0, -0.5);


    for (int ii = 0; ii < numDots; ii++)
    {
        int X = naxis1/4 + ii*naxis1/(2*numDots);
        int Y = naxis2/4 + ii*naxis2/(2*numDots);
        /// Left-hand slices (no copy)
        //matSlice1 = matImage(cv::Range(Y, Y+1), cv::Range(0, naxis1/4));
        matSlice1 = matImage(cv::Range(Y, Y+1), cv::Range(0, naxis1/4));
        cv::blur(matSlice1, matSlice2, cv::Size(7, 7));
        cv::filter2D(matSlice2, gradXLeft, -1, kernelXLeft, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradXLeft, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii] = maxLoc.x;
        dat->Y[ii] = Y;
        /// Right-hand slices (no copy)
        matSlice1 = matImage(cv::Range(Y, Y+1), cv::Range(3*naxis1/4, naxis1));
        cv::blur(matSlice1, matSlice2, cv::Size(7, 7));
        cv::filter2D(matSlice2, gradXRight, -1, kernelXRight, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradXRight, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + numDots] = maxLoc.x + 3*naxis1/4;
        dat->Y[ii + numDots] = Y;
        ///Bottom slices (no copy)
        matSlice1 = matImage(cv::Range(0, naxis2/4), cv::Range(X, X+1));
        cv::blur(matSlice1, matSlice2, cv::Size(7, 7));
        cv::filter2D(matSlice2, gradYBottom, -1, kernelYBottom, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradYBottom, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + 2*numDots] = X;
        dat->Y[ii + 2*numDots] = maxLoc.y;
        ///Top slices (no copy)
        matSlice1 = matImage(cv::Range(3*naxis2/4, naxis2), cv::Range(X, X+1));
        cv::blur(matSlice1, matSlice2, cv::Size(7, 7));
        cv::filter2D(matSlice2, gradYTop, -1, kernelYTop, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradYTop, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + 3*numDots] = X;
        dat->Y[ii + 3*numDots] = maxLoc.y + 3*naxis2/4;

    }
}


QList<RMat *> RProcessing::normalizeSeriesByStats(QList<RMat*> rMatImageList)
{
    QList<RMat*> normalizedRMatImageList;

    for (int i =0 ; i < rMatImageList.size() ; ++i)
    {
        normalizedRMatImageList << normalizeByStats(rMatImageList.at(i));
        normalizedRMatImageList.at(i)->setImageTitle(normalizedRMatImageList.at(i)->getImageTitle() + QString("# %1").arg(i));
        normalizedRMatImageList.at(i)->setDate_time(rMatImageList.at(i)->getDate_time());
    }

    return normalizedRMatImageList;
}

void RProcessing::blurRMat(RMat *rMat)
{
    cv::Mat blurMat;
    cv::blur(rMat->matImage, blurMat, cv::Size(5, 5));

    blurMat.convertTo(blurMat, CV_8U, 256.0f/rMat->getDataRange());
    QImage blurImage(blurMat.data, blurMat.cols, blurMat.rows, QImage::Format_Grayscale8);
    emit resultQImageSignal(blurImage);
}

RMat* RProcessing::normalizeByStats(RMat *rMat)
{
    /// This is like normalizeByThresh() using newMin and newMax as intensityLow and
    /// intensityHigh from RMat::calcStats()

    cv::Mat matImage = normalizeByThresh(rMat->matImage, rMat->getIntensityLow(), rMat->getIntensityHigh(), rMat->getNormalizeRange());

    RMat* normalizeRMat = new RMat(matImage, rMat->isBayer(), rMat->getInstrument());
    normalizeRMat->setImageTitle(QString("Normalized image "));
    showMinMax(matImage);

    return normalizeRMat;
}

void RProcessing::normalizeByStatsInPlace(RMat *rMat)
{
    /// Overload of normalizeByStats(RMat *rMat) for the series "in place", no copy.

    int type = rMat->matImage.type();
    float newMin = rMat->getIntensityLow();
    float newMax = rMat->getIntensityHigh();
    float newDataRange = newMax - newMin;
    float alpha = rMat->getDataRange() / newDataRange;
    float beta = -newMin * rMat->getDataRange() /newDataRange;

    rMat->matImage.convertTo(rMat->matImage, CV_32F);
    rMat->matImage.convertTo(rMat->matImage, type, alpha, beta);
}



cv::Mat RProcessing::normalizeByThresh(cv::Mat matImage, float oldMin, float oldMax, float newRange)
{
    /// This function do contrast stretching and clips the intensity between newMin and newMax.
    /// Contrast Stretching formula from : http://homepages.inf.ed.ac.uk/rbf/HIPR2/stretch.htm
    /// This assumes the lower Range in newRange equals 0.

    float oldRange = oldMax - oldMin;
    float alpha = newRange / oldRange;
    float beta = -oldMin * newRange /oldRange;

    cv::Mat normalizedMatImage;
    matImage.convertTo(normalizedMatImage, CV_32F);

    if (matImage.type() == CV_32F || matImage.type() == CV_32FC3)
    {
        normalizedMatImage.convertTo(normalizedMatImage, CV_16U, alpha, beta);
    }
    else
    {
        normalizedMatImage.convertTo(normalizedMatImage, matImage.type(), alpha, beta);
    }

    return normalizedMatImage;
}

cv::Mat RProcessing::normalizeClipByThresh(cv::Mat matImage, float newMin, float newMax, float dataRange)
{
    /// This function do contrast stretching and clips the intensity between newMin and newMax.
    float newDataRange = newMax - newMin;
    float alpha = dataRange / newDataRange;
    float beta = -newMin * dataRange /newDataRange;

    cv::Mat normalizedMatImage;
    matImage.convertTo(normalizedMatImage, CV_32F);

    // Now we need to clip the image between the max and min of the extrema of the instrument data type range.
    cv::threshold(normalizedMatImage, normalizedMatImage, newMax, newMax, cv::THRESH_TRUNC);
//    normalizedMatImage = newMax - normalizedMatImage;
//    cv::threshold(normalizedMatImage, normalizedMatImage, newMax - newMin, newMax - newMin, cv::THRESH_TRUNC);
//    normalizedMatImage = newMax - normalizedMatImage;

    normalizedMatImage.convertTo(normalizedMatImage, matImage.type(), alpha, beta);
    return normalizedMatImage;
}

void RProcessing::fixUset(cv::Mat matImage)
{
    for (int x = 0; x < 2048; x++)
    {
        matImage.at<uchar>(0, x) = 0;
        matImage.at<uchar>(1, x) = 0;
        matImage.at<uchar>(2, x) = 0;
    }
}


// Private member functions
void RProcessing::normalizeFlat()
{
    //    cv::Mat tempMat;
    //    rMatImage->matImage.copyTo(tempMat);
    masterFlatN = new RMat(*masterFlat);
    // Normalize image by mean value
    cv::Scalar meanValue = cv::mean(masterFlat->matImage);
    float meanValueF = (float) meanValue.val[0];
    qDebug() << "ProcessingWidget:: masterFlatN->matImage.channels()=" << masterFlatN->matImage.channels();
    qDebug() << "ProcessingWidget:: masterFlatN.matImage.type()=" << masterFlatN->matImage.type();
    qDebug() << "ProcessingWidget:: meanValueF=" << meanValueF;
    //cv::normalize(rMatImage->matImage, rMatImageNorm.matImage, 0.5, 1.5, cv::NORM_MINMAX, CV_32F);
    masterFlatN->matImage.convertTo(masterFlatN->matImage, CV_32F);
    masterFlatN->matImage = masterFlatN->matImage / meanValueF; // Strange...
    masterFlatN->calcMinMax();
    masterFlatN->setImageTitle(QString("Master Flat normalized"));
    float bScale = (float) masterFlat->getDataMax() / meanValueF;
    masterFlatN->setBscale(bScale);
    masterFlatN->calcMinMax();
    masterFlatN->calcStats();
}



// setters
void RProcessing::setTreeWidget(RTreeWidget *treeWidget)
{
    this->treeWidget = treeWidget;
}

void RProcessing::setCurrentROpenGLWidget(ROpenGLWidget *rOpenGLWidget)
{
    currentROpenGLWidget = rOpenGLWidget;
}

void RProcessing::setShowContours(bool status)
{
    showContours = status;
}

void RProcessing::setShowLimb(bool status)
{
    showLimb = status;
}

void RProcessing::setUseXCorr(bool useXCorr)
{
    this->useXCorr = useXCorr;
}

void RProcessing::setCvRectROI(cv::Rect cvRect)
{
    this->cvRectROI = cvRect;
}

void RProcessing::setUseROI(bool status)
{
    this->useROI = status;
}

void RProcessing::setBlurSigma(double sigma)
{
    this->blurSigma = sigma;
}

void RProcessing::setUseHPF(bool status)
{
    this->useHPF = status;
}

void RProcessing::setHPFSigma(double sigma)
{
    this->hpfSigma = sigma;
}

void RProcessing::setupMasterWithSigmaClip(bool enabled)
{
    this->masterWithSigmaClip = enabled;
    this->masterWithMean = !masterWithSigmaClip;
    qDebug() << "RProcessing:: masterWithSigmaClip = " << masterWithSigmaClip;
    qDebug() << "RProcessing:: masterWithMean = " << masterWithMean;
}

void RProcessing::setupMasterWithMean(bool enabled)
{
    this->masterWithMean = enabled;
    this->masterWithSigmaClip = !masterWithMean;
    qDebug() << "RProcessing:: masterWithSigmaClip = " << masterWithSigmaClip;
    qDebug() << "RProcessing:: masterWithMean = " << masterWithMean;
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

RMat *RProcessing::getEllipseRMat()
{
    return ellipseRMat;
}

QImage *RProcessing::getCannyQImage()
{
    return cannyQImage;
}

QList<RMat *> RProcessing::getContoursRMatList()
{
    return contoursRMatList;
}

QList<RMat *> RProcessing::getResultList()
{
    return resultList;
}

QList<RMat *> RProcessing::getResultList2()
{
    return resultList2;
}

QList<RMat *> RProcessing::getLimbFitResultList1()
{
    return limbFitResultList1;
}

QList<RMat *> RProcessing::getLimbFitResultList2()
{
    return limbFitResultList2;
}

QVector<Circle> RProcessing::getCircleOutList()
{
    return circleOutList;
}

RMat *RProcessing::getCannyRMat()
{
    return cannyRMat;
}

RMat *RProcessing::getContoursRMat()
{
    return contoursRMat;
}

bool RProcessing::compareContourAreas(std::vector<cv::Point> contour1, std::vector<cv::Point> contour2)
{
    double i = fabs( cv::contourArea(cv::Mat(contour1)) );
    double j = fabs( cv::contourArea(cv::Mat(contour2)) );
    return ( i < j );
}
