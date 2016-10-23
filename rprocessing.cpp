#include "rprocessing.h"

#include <QFileDialog>

// cifitsio
#include <cfitsio/fitsio.h>

//opencv
#include <opencv2/world.hpp>

#include "imagemanager.h"
#include "parallelcalibration.h"
#include "typedefs.h"

using namespace af;

RProcessing::RProcessing(QObject *parent): QObject(parent),
    masterBias(NULL), masterDark(NULL), masterFlat(NULL), masterFlatN(NULL), stackedRMat(NULL), useROI(false), useXCorr(false),
    radius(0), radius1(0), radius2(0), radius3(0), meanRadius(0), masterWithMean(true), masterWithSigmaClip(false), stackWithMean(true),
    stackWithSigmaClip(false)
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
        char keyNameRadius[] = "SOLAR_R";
        float keyValueRadius = meanRadius;
        fits_write_key(fptr, TFLOAT, keyNameRadius, &keyValueRadius, NULL, &status);

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

cv::Mat RProcessing::histogram(cv::Mat matVector, int &nBins, float &width)
{
    cv::Mat matHistogram;
    double min, max;
    cv::minMaxLoc(matVector, &min, &max);
    float range[2];
    range[0] = (float) min;
    range[1] = (float) max;
    nBins = roundf((max - min)/ width);
    const float* histRange = { range };
    bool uniform = true;
    bool accumulate = false;

    cv::calcHist( &matVector, 1, 0, cv::Mat(), matHistogram, 1, &nBins, &histRange, uniform, accumulate);

    return matHistogram;
}

float RProcessing::calcMedian(std::vector<float> data, float width)
{
    cv::Mat matVector(data, false);
    double min, max;
    cv::minMaxLoc(matVector, &min, &max);

    int nBins;
    cv::Mat matHistogram = histogram(matVector, nBins, width);

    float cdf = 0;
    int totalCounts = data.size();
    float medianVal;

    for (int i = 1; i < nBins ; i++)
    {
        cdf += matHistogram.at<float>(i);

        if (cdf / totalCounts >= 0.5)
        {
            medianVal = i;
            break;
        }
    }
    /// We have to recover the actual data value
    /// that falls within that median bin value
    medianVal = min + medianVal*width;

    return medianVal;
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

void RProcessing::stack(QList<RMat *> rMatImageList)
{
    if (rMatImageList.isEmpty())
    {
        tempMessageSignal(QString("No image to stack"));
        return;
    }

    if (stackWithMean)
    {
        stackedRMat = average(rMatImageList);
    }
    else if (stackWithSigmaClip)
    {
        stackedRMat = sigmaClipAverage(rMatImageList);
    }

    if (stackedRMat == NULL)
    {
        tempMessageSignal(QString("Stacking failed"));
        return;
    }

    emit resultSignal(stackedRMat);

    return;
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
    avgImg.convertTo(avgImg, rMatList.at(0)->matImage.type());
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

    try {
        af::setBackend(AF_BACKEND_OPENCL);
//        af::array tempArray(10, 10, 2);
//        tempArray(af::span, af::span, 0) = 1;
//        tempArray(af::span, af::span, 1) = 0;

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
//    qDebug("RProcessing::sigmaClipAverage::  elapsed seconds: %f us", af::timer::stop(start2));
//    af::array meanArfTiled = af::tile(meanArf, 1, 1, nFrames);
//    //af::array stdevArf = af::moddims(af::stdev(arfSeries, 2), naxis2, naxis1);
//    af::array stdevArf = af::stdev(arfSeries, 2);
//    af::array stdevArfTiled = af::tile(stdevArf, 1, 1, nFrames);
//    af::array arfMaskReject = af::abs(arfSeries - meanArfTiled) > stdevArfTiled;
//    arfSeries(arfMaskReject) = meanArfTiled(arfMaskReject);
//    meanArf = af::mean(arfSeries, 2);
//    qDebug("RProcessing::sigmaClipAverage::  elapsed seconds: %f s", af::timer::stop(start2));
    /// Copy an array from the device to the host:
    float *hostArf = meanArf.host<float>();

    qDebug("RProcessing::sigmaClipAverage::  elapsed seconds: %f s", af::timer::stop(start2));
    /// Prepare output Mat image
    cv::Mat matImage(naxis2, naxis1, CV_32F, hostArf);
    matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
    RMat *rMatAvg = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());

    rMatAvg->setSOLAR_R(rMatImageList.at(0)->getSOLAR_R());


    return rMatAvg;

    } catch (af::exception& e) {
        printf("caught exception when trying CPU backend\n");
        fprintf(stderr, "%s\n", e.what());
    }
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

QList<RMat *> RProcessing::sharpenSeries(QList<RMat*> rMatImageList, float weight1, float weight2)
{
    QList<RMat*> rMatSharpList;

    for (int i = 0 ; i < rMatImageList.size() ; i++)
    {
        RMat *rMatSharp = sharpenCurrentImage(rMatImageList.at(i), weight1, weight2);
        rMatSharpList << normalizeByStats(rMatSharp);
        rMatSharpList.at(i)->setImageTitle(QString("Sharpened image # %1").arg(i));
    }

    return rMatSharpList;
}

RMat* RProcessing::sharpenCurrentImage(RMat *rMatImage, float weight1, float weight2)
{
    cv::Mat blurredImage;
    cv::Mat matImageSharpened;

    cv::GaussianBlur(rMatImage->matImage, blurredImage, cv::Size(0, 0), 3);

    cv::addWeighted(rMatImage->matImage, weight1, blurredImage, -weight2, 0, matImageSharpened);
//        cv::Mat matDiff = rMatImageList.at(i)->matImage - blurredImage;
//        cv::addWeighted(rMatImageList.at(i)->matImage, 1.0, matDiff, weight2, 0, matImageSharpened);
    RMat *rMatSharp = new RMat(matImageSharpened, false, rMatImage->getInstrument());
    RMat *rMatSharpN = normalizeByStats(rMatSharp);
    rMatSharpN->setSOLAR_R(rMatImage->getSOLAR_R());

    return rMatSharpN;
}

void RProcessing::blockProcessingLocal(QList<RMat*> rMatImageList)
{

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    setBackend(AF_BACKEND_OPENCL);

    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();

    /// Initialize the GPU array series
    int binning = 2;
    int nBAxis1 = naxis1/binning;
    int nBAxis2 = naxis2/binning;
    af::array arfSeries(naxis2, naxis1, nFrames);
    //af::array arfSeriesBinned(nBAxis2, nBAxis1, nFrames);
    af::array gradientMagSeries(nBAxis2, nBAxis1, nFrames);


    /// gaussiankernel (const int rows, const int cols, const double sig_r=0, const double sig_c=0)
//    af::array kernel = gaussianKernel(3, 3, 1, 1);
//    af::array kernel = af::constant(0, 3, 3);
//    kernel(1, 1) = 1;
//    kernel(1, 2) = 1;
//    kernel(2, 1) = 1;
//    kernel(2, 2) = 1;

    af::array kernel = af::constant(0, 7, 7);
    kernel(seq(3,6), seq(3,6)) = 1;
    af_print(kernel);

    /// Store the matImageList (converted to float) in the GPU array, and get a binned version
    ///

    af::timer afTimer1 = af::timer::start();

    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempArf(naxis2, naxis1, (float*) tempMat.data);
        arfSeries(span, span, k) = tempArf;
        af::array  arfTemp = convolve2(tempArf, kernel);
        af::array binnedFrame = arfTemp(seq(0, af::end, binning), seq(0, af::end, binning));
        af::array dx, dy;
        af::grad(dx, dy, binnedFrame);
        /// Sum of absolute of the X- and Y- gradient (look for a norm-L2 function in ArrayFire?)
        af::array gradientMagnitude = af::abs(dx) + af::abs(dy);
        gradientMagSeries(span, span, k) = gradientMagnitude;
        /// Sobel gradient magnitude
        ///gradientMagSeries(span, span, k) = af::sobel(binnedFrame);
    }
    qDebug("RProcessing::blockProcessingLocal:: afTimer1 elapsed seconds: %f s", af::timer::stop(afTimer1));

    /// Sobel gradient. Comment this if using the classical gradient above.
//    af::timer afTimer2 = af::timer::start();
//    af::array gradientMagSeries = af::sobel(arfSeriesBinned);
//    qDebug("RProcessing::blockProcessing:: afTimer2 elapsed seconds: %f s", af::timer::stop(afTimer2));


    /// Calculate number of blocks in the binned image. There must be enough block to cover the whole image, considering an overlap of 50%;
    int blkSize = 32;
    int blkSizeL = blkSize * 3;
    int binnedBlkSize = blkSize/binning;
    int nBBlksX = nBAxis1 / binnedBlkSize * 2; /// *2 from an overalp of 50% in the X-direction;
    int nBBlksY = nBAxis2 / binnedBlkSize;
    //int nBPixels = nBAxis1 * nBAxis2;
    int nBBlks = nBBlksX * nBBlksY;



//    af::timer afTimer3 = af::timer::start();

//    af::array canvas = constant(0, naxis2, naxis1);

//    /// Loop over all the block series
//    for (int k = 0; k < nBBlks; k++)
//    {
//        int blkOffsetX = (k % nBBlksY) * binnedBlkSize;
//        int blkOffsetY = (k / nBBlksX) * binnedBlkSize;

//        af::array binnedBlkSeries = moddims(gradientMagSeries(seq(blkOffsetY, blkOffsetY + binnedBlkSize -1), seq(blkOffsetX, blkOffsetX + binnedBlkSize -1), seq(0, nFrames-1)), dim4(binnedBlkSize*binnedBlkSize, nFrames));
//        af::array gradSum = flat(af::sum( binnedBlkSeries, 0 ));
//        /// Sort the sum of the gradient-norm
//        af::array sortedArray;
//        af::array sortIndices;
//        af::sort(sortedArray, sortIndices, gradSum, 0, false);

//        blkOffsetX *= binning;
//        blkOffsetY *= binning;

//        af::array blkSeries = arfSeries(seq(blkOffsetY, blkOffsetY + blkSize -1), seq(blkOffsetX, blkOffsetX + blkSize -1), seq(0, nFrames-1));
//        af::array sortedBlk = blkSeries(span, span, sortIndices);
//        af::array bestBlock = sortedBlk(span, span, 0);
//        canvas(seq(blkOffsetY, blkOffsetY + blkSize -1), seq(blkOffsetX, blkOffsetX + blkSize -1)) = bestBlock;
//    }

//    qDebug("RProcessing::blockProcessing:: afTimer3 elapsed seconds: %f s", af::timer::stop(afTimer3));


    af::timer afTimer4 = af::timer::start();

    af::array canvas = constant(0, naxis1, naxis2, f32);
    af::array canvasWeights = constant(0, naxis1, naxis2, f32);

    /// Loop over all the block series
    ///
    int offsetX = 500;
    int offsetY = 500;

    int x1 = offsetX;
    int x2 = x1;
    int y1 = offsetY;
    int y2 = y1;
    int xL1, yL1, xL2, yL2;
    int dx = 0;
    int dy = 0;
    int x1_Bottom = offsetX;
    int y1_Bottom = offsetY;
    int dx_Bottom = 0;
    int dy_Bottom = 0;

    QPoint blkPos1(x1, y1);
//    qDebug("[x1 , y1] = [%d, %d]", x1, y1);
    qDebug("y1 = %d", y1);
    af::array searchBlk1, searchBlk2, bestBlk1, bestBlk2;
    /// Create a swappable template for dealing with the blocks at the edge of the image who need a bottom template
    af::array searchBlk1Bottom;
    sortBestBlocks(bestBlk1, searchBlk1, arfSeries, gradientMagSeries, blkSize, blkSizeL, blkPos1, nFrames, binning);
    searchBlk1Bottom = searchBlk1;

    /// Put the best block on the new canvas. Careful with the dimensions, they are not consistent with cv::Mat
    /// We have to assimilate the rows as the x-direction, so we have 2D af::array as array[x,y]
    canvas(seq(blkPos1.x(), blkPos1.x() + blkSize -1), seq(blkPos1.x(), blkPos1.x() + blkSize -1)) = bestBlk1;
    //canvasWeights(seq(blkPos1.x(), blkPos1.x() + blkSize -1), seq(blkPos1.x(), blkPos1.x() + blkSize -1)) += 1;


    int k = 0;
    //while (k < 2)
    while ((x2 < 1500) && (y2 < 600))
    {
        /// If the end of the 2nd block is off the edge, we go back to the beginning of the x-axis but step-up in the y-axis;
        if ((x2 + blkSize > 1500) )
        {

            x1 = offsetX;
            y1 += blkSize/2;

            if (y1 > 600)
            {
                break;
            }
            blkPos1 = QPoint(x1, y1);

            //qDebug("Changing to row at y1 = %d", y1);
            //qDebug("y1 = %d", y1);

            sortBestBlocks(bestBlk1, searchBlk1, arfSeries, gradientMagSeries, blkSize, blkSizeL, blkPos1, nFrames, binning);

            /// Template matching
            af::array SAD = matchTemplate(searchBlk1Bottom, bestBlk1);
            //af_print(SAD);
            af::array idx, idy;
            af::array minValuesX, minValuesY;
            af::min(minValuesX, idx, SAD);
            af::min(minValuesY, idy, minValuesX, 1);
            af::array minY = flat(idy);
            af::array minX = flat(idx(minY));

            unsigned int *dx1 = minX.host<unsigned int>();
            unsigned int *dy1 = minY.host<unsigned int>();
            //qDebug("[dx1, dy1] = [%u, %u]", dx1[0], dy1[0]);

            xL1 = x1_Bottom - blkSize;
            yL1 = y1_Bottom - blkSize; /// Investigate possible correction needed for fixing the relative position of the top block with respect to the bottom block, on the 3rd row, because the bottom block starts shifts as of 2nd row...
            int x1Shifted = xL1 + (int) dx1[0] + dx_Bottom;
            int y1Shifted = yL1 + (int) dy1[0] + dy_Bottom;

            QPoint newBlkPos1(x1Shifted, y1Shifted);

            canvas(seq(newBlkPos1.x(), newBlkPos1.x() + blkSize -1), seq(newBlkPos1.y(), newBlkPos1.y() + blkSize -1)) = bestBlk1;
            //canvasWeights(seq(newBlkPos1.x(), newBlkPos1.x() + blkSize -1), seq(newBlkPos1.y(), newBlkPos1.y() + blkSize -1)) += 1;

            searchBlk1Bottom = searchBlk1;

            dx_Bottom = x1Shifted - x1;
            dy_Bottom = y1Shifted - y1;
            x1_Bottom = x1;
            y1_Bottom = y1;

            /// The new block on the left-side of the row is also the new template for the next (2nd) block in that row. So its correct position need to be adujusted by an updated dx- and dy- shift.
            dx = dx_Bottom;
            dy = dy_Bottom;
        }

        /// X-position of block 2, overlapping with block 1 by 50%. Y-position of block 2, same as block 1;
        x2 = x1 + blkSize/2;
        y2 = y1;
        QPoint blkPos2(x2, y2);

        sortBestBlocks(bestBlk2, searchBlk2, arfSeries, gradientMagSeries, blkSize, blkSizeL, blkPos2, nFrames, binning);

        /// Template matching
        af::array SAD = matchTemplate(searchBlk1, bestBlk2);
        //af_print(SAD);
        af::array idx, idy;
        af::array minValuesX, minValuesY;
        af::min(minValuesX, idx, SAD);
        af::min(minValuesY, idy, minValuesX, 1);
        af::array minY = flat(idy);
        af::array minX = flat(idx(minY));

        unsigned int *dx2 = minX.host<unsigned int>();
        unsigned int *dy2 = minY.host<unsigned int>();
        /// new x2 = xL1 + d2 + dx (see Lucky Imaging slides 3-4)
        xL1 = x1 - blkSize;
        yL1 = y1 - blkSize;
        int x2shifted = xL1 + (int) dx2[0] + dx;
        int y2shifted = yL1 + (int) dy2[0] + dy;
        QPoint shiftedBlkPos2(x2shifted, y2shifted);

//        qDebug("[x2shifted, y2shifted] = [%d, %d]", x2shifted, y2shifted);

        canvas(seq(shiftedBlkPos2.x(), shiftedBlkPos2.x() + blkSize -1), seq(shiftedBlkPos2.y(), shiftedBlkPos2.y() + blkSize -1)) = bestBlk2;
        //canvasWeights(seq(shiftedBlkPos2.x(), shiftedBlkPos2.x() + blkSize -1), seq(shiftedBlkPos2.y(), shiftedBlkPos2.y() + blkSize -1)) += 1;

        /// Advance to block2 as base of the new template block.
        searchBlk1 = searchBlk2;
        dx = (x2shifted - x2);
        dy = (y2shifted - y2);


        x1 = x2;
        blkPos1 = QPoint(x1, y1);

//        qDebug("[x2shifted, y2shifted] = [%d, %d]", x2shifted, y2shifted);
//        qDebug("[dx, dy] = [%d, %d]", dx, dy);

        delete[] dx2;
        delete[] dy2;

        //qDebug("[dx, dy] = %d, %d", dx, dy);
        k++;

//        if (kk == 1)
//        {
//            break;
//        }
    }

    //canvas /= canvasWeights;

    qDebug("RProcessing::blockProcessingLocal:: afTimer4 elapsed seconds: %f s", af::timer::stop(afTimer4));

    float *hostArf = canvas.host<float>();
    cv::Mat matImage(naxis2, naxis1, CV_32F, hostArf);
//    cv::Mat matImage(N2, N1, CV_32F, hostArf);
    matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
    RMat *blockRMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
    blockRMat->setImageTitle("with local templates");
    resultList << blockRMat;
    delete[] hostArf;

    /// Display one block series
//    for ( int k=0; k < nFrames; k++)
//    {
//        af::array resultSlice = sortedBlk(span, span, k);
//        //af::array resultSlice = arfSeriesBinned(span, span, k);
//        //cv::Mat matImage(naxis2, naxis1, CV_32F, hostArf);

//        float *hostArf = resultSlice.host<float>();
//        cv::Mat matImage(blkSize, blkSize, CV_32F, hostArf);
//        matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
//        RMat *blockMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
//        luckyList << blockMat;
//        delete[] hostArf;
//    }

    /// Test a problematic template / block pair
//    x1 = 1332;
//    y1 = 596;
//    blkPos1 = QPoint(x1, y1);

//    x2 = x1 + blkSize/2;
//    y2 = y1;
//    QPoint blkPos2(x2, y2);

//    sortBestBlocks(bestBlk1, searchBlk1, arfSeries, gradientMagSeries, blkSize, blkSizeL, blkPos1, nFrames, binning);
//    sortBestBlocks(bestBlk2, searchBlk2, arfSeries, gradientMagSeries, blkSize, blkSizeL, blkPos2, nFrames, binning);

//    /// Template matching
//    af::array SAD = matchTemplate(searchBlk1, bestBlk2);
//    //af_print(SAD);
//    af::array idx, idy;
//    af::array minValuesX, minValuesY;
//    af::min(minValuesX, idx, SAD);
//    af::min(minValuesY, idy, minValuesX, 1);
//    af::array minY = flat(idy);
//    af::array minX = flat(idx(minY));

//    unsigned int *dx2 = minX.host<unsigned int>();
//    unsigned int *dy2 = minY.host<unsigned int>();

//    int xShift = (int) dx2[0] - (blkSize + blkSize/2);
//    int yShift = (int) dy2[0] - blkSize;

//    int x2Shifted = blkSize/2 + xShift;
//    int Yoffset2 = 10;
//    int y2Shifted = Yoffset2 + yShift;


//    qDebug("[x2Shifted , y2Shifted] = [%d, %d]", x2Shifted, y2Shifted);
//    qDebug("[xShift , yShift] = [%d, %d]", xShift, yShift);

//    af::array testCanvas = constant(0, 3*blkSize, 3*blkSize);
//    testCanvas(seq(0, blkSize-1), seq(Yoffset2, Yoffset2 + blkSize-1)) = bestBlk1;
//    testCanvas(seq(x2Shifted, x2Shifted + blkSize -1), seq(y2Shifted, y2Shifted + blkSize -1)) = bestBlk2;


//    float *hostArfTest = testCanvas.host<float>();
//    cv::Mat matImageTest(3*blkSize, 3*blkSize, CV_32F, hostArfTest);
//    matImageTest.convertTo(matImageTest, rMatImageList.at(0)->matImage.type());
//    RMat *blockRMatTest = new RMat(matImageTest, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());

//    if (!resultList2.empty())
//    {
//        qDeleteAll(resultList2);
//        resultList2.clear();
//    }

//    resultList2 << blockRMatTest;
//    delete[] hostArfTest;

}



void RProcessing::sortBestBlocks(af::array & bestBlk, af::array & searchBlk, af::array & arraySeries, af::array & arrayBinnedSeries, int blkSize, int blkSizeL, QPoint blkPos, int nFrames, int binning)
{
    QPoint blkPosB = blkPos/binning;
    int binnedBlkSize = blkSize / binning;

    af::array binnedBlkSeries = moddims(arrayBinnedSeries(seq(blkPosB.x(), blkPosB.x() + binnedBlkSize -1), seq(blkPosB.y(), blkPosB.y() + binnedBlkSize -1), seq(0, nFrames-1)), dim4(binnedBlkSize*binnedBlkSize, nFrames));
    af::array gradSum = flat(af::sum( binnedBlkSeries, 0 ));
    //af_print(gradSum);

    /// Sort the sum of the gradient-norm
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, gradSum, 0, false);
    //af_print(sortedArray);

    /// Extract the larger version of the best block.
    int xL = blkPos.x() - blkSize;
    int yL = blkPos.y() - blkSize;

    QPoint blkPosL(xL, yL);
    /// Index of the best block is supposed to be at sortIndices(0)
    //af::array blkSeriesL = arraySeries(seq(blkPosL.x(), blkPosL.x() + blkSizeL -1), seq(blkPosL.y(), blkPosL.y() + blkSizeL -1), seq(0, nFrames-1));
    unsigned int *bestIndex = sortIndices.host<unsigned int>();
    searchBlk = arraySeries(seq(blkPosL.x(), blkPosL.x() + blkSizeL -1), seq(blkPosL.y(), blkPosL.y() + blkSizeL -1), bestIndex[0]);

    //af::array blkSeries = arraySeries(seq(blkPos.x(), blkPos.x() + blkSize -1), seq(blkPos.y(), blkPos.y() + blkSize -1), seq(0, nFrames-1));
    bestBlk = searchBlk(seq(blkSize, 2*blkSize -1), seq(blkSize, 2*blkSize -1));
}

void RProcessing::extractLuckySample(QList<RMat *> rMatImageList, const int &blkSize, const int &binning)
{
    if (!luckyBlkList.empty())
    {
        qDeleteAll(luckyBlkList);
        luckyBlkList.clear();
    }

    af::array arfSeries;
    af::array qualityBinnedSeries;
    blockProcessingSetup(rMatImageList, arfSeries, qualityBinnedSeries, binning);

    int binnedBlkSize = blkSize/binning;

    int nBest = 3;
    int offsetX = 500;
    int offsetY = 500;
    int x = 0;
    int y = 0;

    int k = 0;
    x = offsetX + (k*blkSize/2 % 1000);
    y = offsetY + ((k*blkSize/2) / 1000)*blkSize/2;

    int bufferSpace = 8;
    af::array stackedBlks;
    makeAlignedStack(stackedBlks, arfSeries, qualityBinnedSeries,
                     nBest, blkSize, binnedBlkSize, x, y, bufferSpace, binning);

    /// Display one block series
    for ( int k=0; k < nBest; k++)
    {
        af::array resultSlice = stackedBlks(span, span, k);
        float *hostArf = resultSlice.host<float>();
        cv::Mat matImage(blkSize, blkSize, CV_32F, hostArf);
        matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
        RMat *blockRMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
        blockRMat->setImageTitle("Blocks stack");
        luckyBlkList << blockRMat;
        delete[] hostArf;
    }
}

void RProcessing::blockProcessingGlobal(QList<RMat *> rMatImageList)
{

    if (!resultList2.empty())
    {
        qDeleteAll(resultList2);
        resultList2.clear();
    }

    setBackend(AF_BACKEND_OPENCL);

    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();

    /// Initialize the GPU array series
    int binning = 2;
    int nBAxis1 = naxis1/binning;
    int nBAxis2 = naxis2/binning;
    int blkSize = 32;
    int binnedBlkSize = blkSize/binning;
    int bufferSpace = 10;

    af::array arfSeries(naxis2, naxis1, nFrames);
    //af::array arfSeriesBinned(nBAxis2, nBAxis1, nFrames);
    af::array gradientMagSeries(nBAxis2, nBAxis1, nFrames);

    af::array kernel = af::constant(0, 7, 7);
    kernel(seq(3,6), seq(3,6)) = 1;

    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempArf(naxis2, naxis1, (float*) tempMat.data);
        arfSeries(span, span, k) = tempArf;
        af::array  arfTemp = convolve2(tempArf, kernel);
        af::array binnedFrame = arfTemp(seq(0, af::end, binning), seq(0, af::end, binning));
        af::array dx, dy;
        af::grad(dx, dy, binnedFrame);
        /// Sum of absolute of the X- and Y- gradient (look for a norm-L2 function in ArrayFire?)
        af::array gradientMagnitude = af::abs(dx) + af::abs(dy);
        gradientMagSeries(span, span, k) = gradientMagnitude;
        /// Sobel gradient magnitude
        ///gradientMagSeries(span, span, k) = af::sobel(binnedFrame);
    }


    //    af::array globalTemplate = arfSeries(span, span, 0); //af::median(arfSeries, 2);
    af::array globalTemplate = af::median(arfSeries, 2);
    qDebug("Median template ready...");

    int nBest = 3;
    //af::array canvas0 = constant(0, naxis1, naxis2, f32);
    af::array canvas = constant(0, naxis1, naxis2, f32);
    //af::array canvas2 = constant(0, naxis1, naxis2, nBest, f32);

    //af::array canvasWeights = constant(0, naxis1, naxis2, f32);

    /// Loop over all the block series
    ///
//    int offsetX = blkSize;
//    int offsetY = blkSize;
    int offsetX = 500;
    int offsetY = 500;
    int k = 0;

    af::array searchBlk;
    af::array bestBlk; //, bestBlk2;

    double timeA1 = 0;
    double timeB1 = 0;
    double timeA2 = 0;
    double timeB2 = 0;
    double timeC1 = 0;
    double timeC2 = 0;
    double timeC3 = 0;
    double loopTime = 0;

//    af::array arrOnes = constant(1, dim4(blkSize, blkSize));

    QElapsedTimer timer;
    timer.start();

    int x = 0;
    int y = 0;

    bool processBlocks = true;

    int nBlocks = 437;
       for (int k = 0 ; k <= nBlocks; k++)
       {


//    for (int y = offsetY; y< naxis2 - 1 - blkSizeSeach; y += blkSize/2)
//    {
//        for (int x = offsetX; x< naxis1 - 1 - blkSizeSeach; x += blkSize/2)
//        {
//                for (int y = offsetY; y< 600; y += blkSize/2)
//                {
//                    for (int x = offsetX; x< 1500; x += blkSize/2)
//                    {



            af::sync();
            af::timer afTimer4 = af::timer::start();

            x = offsetX + (k*blkSize/2 % 1000);
            y = offsetY + ((k*blkSize/2) / 1000)*blkSize/2;

            qDebug("Block k = %d; [x, y] = [%d, %d]", k, x, y);


            int xRef1 = x - bufferSpace/2;
            int yRef1 = y - bufferSpace/2;
            int xRef2 = xRef1 + blkSize + bufferSpace - 1;
            int yRef2 = yRef1 + blkSize + bufferSpace - 1;


            /// Extract the best block from the current block series

            //            af::sync();
            //            af::timer afTimerA2 = af::timer::start();
            extractBestBlock(bestBlk, arfSeries, gradientMagSeries,
                              blkSize, binnedBlkSize, x, y, nFrames, binning);
            //            af::sync();
            //            timeA2 += af::timer::stop(afTimerA2);

            /// Need to extract the template block from the global template.
            searchBlk  = globalTemplate(seq(xRef1, xRef2), seq(yRef1, yRef2));

//            af::sync();
//            af::timer afTimerC2 = af::timer::start();
            af::array SAD;
            matchTemplate2(SAD, searchBlk, bestBlk);
//            af::sync();
//            timeC2 += af::timer::stop(afTimerC2);
            //qDebug("[dx2, dy2] = [%d, %d]", dx2, dy2);

            int dx, dy, xShifted, yShifted;
            fetchTMatch2Shifts(SAD, dx, dy, bestBlk.dims());
            xShifted = x + (dx - bufferSpace/2);
            yShifted = y + (dy - bufferSpace/2);
            qDebug("[dx, dy] = [%d, %d]", dx, dy);

            canvas(seq(xShifted, xShifted + blkSize -1), seq(yShifted, yShifted + blkSize -1)) = bestBlk;

            af::sync();
            loopTime += af::timer::stop(afTimer4);


//        }
//    }
       }

    qDebug("blockProcessingGlobal:: QTimer : total time: %f s", (float) timer.elapsed() /1000);


    qDebug("blockProcessingGlobal:: Total blocks processed: %d", k);
    qDebug("blockProcessingGlobal:: templateMatching time per block timeA1: %f ms", timeA1 / nBlocks * 1000);
    qDebug("blockProcessingGlobal:: templateMatching time per block timeA2: %f ms", timeA2 / nBlocks * 1000);

    qDebug("blockProcessingGlobal::  time per block timeB1: %f ms", timeB1 / nBlocks * 1000);
    qDebug("blockProcessingGlobal::  time per block timeB2: %f ms", timeB2 / nBlocks * 1000);

    qDebug("blockProcessingGlobal:: timeC1: %f ms", timeC1 / nBlocks * 1000);
    qDebug("blockProcessingGlobal:: timeC2: %f ms", timeC2 / nBlocks * 1000);
    qDebug("blockProcessingGlobal:: timeC3: %f ms", timeC3 / nBlocks *1000);

    qDebug("blockProcessingGlobal:: Loop time per block: %f ms", loopTime/nBlocks * 1000);

//    float *hostArf0 = canvas0.host<float>();
//    cv::Mat matImage0(naxis2, naxis1, CV_32F, hostArf0);
//    matImage0.convertTo(matImage0, rMatImageList.at(0)->matImage.type());
//    RMat *blockRMat0 = new RMat(matImage0, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
//    blockRMat0->setImageTitle("with global templates - matchTemplate");
//    resultList2 << blockRMat0;
//    delete[] hostArf0;


    float *hostArf = canvas.host<float>();
    cv::Mat matImage(naxis2, naxis1, CV_32F, hostArf);
    matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
    RMat *blockRMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
    blockRMat->setImageTitle("with global templates (1 block)");
    resultList2 << blockRMat;
    delete[] hostArf;



//    float *hostArf2 = bestBlk.host<float>();
//    cv::Mat matImage2(blkSize, blkSize, CV_32F, hostArf2);
//    matImage2.convertTo(matImage2, rMatImageList.at(0)->matImage.type());
//    RMat *blockRMat2 = new RMat(matImage2, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
//    resultList2 << blockRMat2;
    //    delete[] hostArf2;
}

void RProcessing::blockProcessingSetup(QList<RMat *> rMatImageList, af::array &arfSeries, af::array &qualityBinnedSeries, const int &binning)
{
    setBackend(AF_BACKEND_OPENCL);

    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();

    /// Initialize the GPU array series
    int nBAxis1 = naxis1/binning;
    int nBAxis2 = naxis2/binning;

    arfSeries = constant(0, naxis1, naxis2, nFrames);
    qualityBinnedSeries = constant(0, nBAxis1, nBAxis2, nFrames);

    af::array kernel = af::constant(0, 7, 7);
    kernel(seq(3,6), seq(3,6)) = 1;

    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempArf(naxis2, naxis1, (float*) tempMat.data);
        arfSeries(span, span, k) = tempArf;
        af::array  arfTemp = convolve2(tempArf, kernel);
        af::array binnedFrame = arfTemp(seq(0, af::end, binning), seq(0, af::end, binning));
        af::array dx, dy;
        af::grad(dx, dy, binnedFrame);
        /// Sum of absolute of the X- and Y- gradient (look for a norm-L2 function in ArrayFire?)
        af::array gradientMagnitude = af::abs(dx) + af::abs(dy);
        qualityBinnedSeries(span, span, k) = gradientMagnitude;
        /// Sobel gradient magnitude
        ///gradientMagSeries(span, span, k) = af::sobel(binnedFrame);
    }

    //    af::array globalTemplate = arfSeries(span, span, 0); //af::median(arfSeries, 2);


}

void RProcessing::blockProcessingGlobalStack1(QList<RMat *> rMatImageList, const int &blkSize, const int &binning)
{


    if (!luckyBlkList.empty())
    {
        qDeleteAll(luckyBlkList);
        luckyBlkList.clear();
    }

    af::array arfSeries;
    af::array qualityBinnedSeries;
    blockProcessingSetup(rMatImageList, arfSeries, qualityBinnedSeries, binning);
    af::array globalRefAr = af::median(arfSeries, 2);
    qDebug("Median template ready...");
    af::array canvas = constant(0, globalRefAr.dims(0), globalRefAr.dims(1), f32);
    af::array canvas2 = constant(0, globalRefAr.dims(0), globalRefAr.dims(1), f32);
    af::array canvasWeights = constant(0, globalRefAr.dims(0), globalRefAr.dims(1), f32);
    af::array checkBlks;

    float radius = rMatImageList.at(0)->getSOLAR_R();
    float radiusMargin = 200;
    int binnedBlkSize = blkSize/binning;
    int bufferSpace = 8;

    int nBest = 3;
    int offsetX = blkSize;//500;
    int offsetY = blkSize;//500;
    int x = 0;
    int y = 0;

    int nBlocks = 8000;//437;
    int naxis1 = arfSeries.dims(0);
    int naxis2 = arfSeries.dims(1);
    int naxis1End = naxis1 - 1 - 3*blkSize;
    int naxis2End = naxis2 - 1 - 3*blkSize;

    //af::array cropRange = tile(af::range(blkSize), 1, nBest) + bufferSpace/2;

    int k = 0;
    af::sync();
    af::timer afTimer = af::timer::start();

    //for (int k = 0 ; k <= nBlocks; k++)
    //while ( (y < 500) || (x < naxis1End) )
    while ( (y < naxis2End) || (x < naxis1End) )
    {
        x = offsetX + (k*blkSize/2 % (naxis1End));
        y = offsetY + ((k*blkSize/2) / (naxis1End))*blkSize/2;

        float distToCenter = sqrt( pow(x - naxis1/2, 2) + pow(y - naxis2/2, 2) );

        if ( distToCenter + blkSize > radius + radiusMargin)
        {
            canvas(seq(x, x + blkSize -1), seq(y, y + blkSize -1)) = globalRefAr(seq(x, x + blkSize -1), seq(y, y + blkSize -1));
            canvasWeights(seq(x, x + blkSize -1), seq(y, y + blkSize -1)) = 1;
            k++;
            continue;
        }

        //qDebug("[x, y] = [%d, %d]", x, y);

        af::array stackedBlks;

        makeAlignedStack(stackedBlks, arfSeries, qualityBinnedSeries, nBest,
                         blkSize, binnedBlkSize, x, y, bufferSpace, binning);

        af::array avgStackedBlk = median(stackedBlks, 2);
        if (k == 1)
        {
            checkBlks = stackedBlks;
        }

        int xRef1 = x - bufferSpace/2;
        int xRef2 = xRef1 + blkSize + bufferSpace - 1;
        int yRef1 = y - bufferSpace/2;
        int yRef2 = yRef1 + blkSize + bufferSpace - 1;
        //        qDebug("[xRef1, xRef2] = [%d, %d]", xRef1, xRef2);
        //        qDebug("[yRef1, yRef2] = [%d, %d]", yRef1, yRef2);


        af::array searchBlk  = globalRefAr(seq(xRef1, xRef2), seq(yRef1, yRef2));

        af::array SAD;
        matchTemplate2(SAD, searchBlk, avgStackedBlk);

        int dx, dy;
        fetchTMatch2Shifts(SAD, dx, dy, avgStackedBlk.dims());
        int xShifted1 = x + dx - bufferSpace/2;
        int xShifted2 = xShifted1 + blkSize -1;
        int yShifted1 = y + dy - bufferSpace/2;
        int yShifted2 = yShifted1 + blkSize -1;

        //        qDebug("[xShifted1, xShifted2] = [%d, %d]", xShifted1, xShifted2);
        //        qDebug("[yShifted1, yShifted2] = [%d, %d]", yShifted1, yShifted2);

        //canvas(seq(xShifted1, xShifted2), seq(yShifted1, yShifted2)) = avgStackedBlk;
        canvas(seq(xShifted1, xShifted2), seq(yShifted1, yShifted2)) += avgStackedBlk;
        canvasWeights(seq(xShifted1, xShifted2), seq(yShifted1, yShifted2)) += 1;

        /// Increment the block number.
        k++;
    }

    af::array mask = canvasWeights == 0;
    canvasWeights(mask) = 1;
    canvas /= canvasWeights;
    // Fill in the borders of the canvas that were not populated with the stacked blocks.
    // Left-hand side border
    canvas(seq(0, offsetX + blkSize/2), span) = globalRefAr(seq(0, offsetX + blkSize/2), span);
    // Bottom side border
    canvas(span, seq(0, offsetY + blkSize/2)) = globalRefAr(span, seq(0, offsetY + blkSize/2));
    // Right-hand side border
    canvas(seq(naxis1End, af::end), span) = globalRefAr(seq(naxis1End, af::end), span);
    // Top side border
    canvas(span, seq(naxis2End, af::end)) = globalRefAr(span, seq(naxis2End, af::end));

    af::sync();
    double totalTime = af::timer::stop(afTimer);


    int k2 = k -1;
    qDebug("ProcessingGlobalStack1:: total time for %d loops = %f s", k2, totalTime);
    qDebug("ProcessingGlobalStack1:: average time per loop = %f ms", totalTime / k2 * 1000.0);

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    /// Slice in the 1st original image
    resultList << rMatImageList.at(0);
    populateResultListWithAr(rMatImageList, canvas, QString("with averaged overlap"));

}

void RProcessing::extractBestBlock(af::array &bestBlk, af::array &arfSeries, af::array &arrayBinnedSeries,
                                     const int &blkSize, const int &binnedBlkSize, const int &x, const int &y,
                                    const int &nFrames, const int &binning)
{
    int xB = x/binning;
    int yB = y/binning;

    af::array binnedBlkSeries = moddims(arrayBinnedSeries(seq(xB, xB + binnedBlkSize -1), seq(yB, yB + binnedBlkSize -1), seq(0, nFrames-1)), dim4(binnedBlkSize*binnedBlkSize, nFrames));
    af::array gradSum = flat(af::sum( binnedBlkSeries, 0 ));
    //af_print(gradSum);

    /// Sort the sum of the gradient-norm
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, gradSum, 0, false);

    unsigned int *bestIndices = sortIndices.host<unsigned int>();
    bestBlk = arfSeries(seq(x, x + blkSize -1), seq(y, y + blkSize -1), bestIndices[0]);
}

void RProcessing::makeAlignedStack(af::array &stackedBlks, const af::array &arfSeries, const af::array &qualityBinnedSeries, const int nBest, const int &blkSize, const int &binnedBlkSize, int &x, int &y, const int bufferSpace, const int &binning)
{
    int nFrames = arfSeries.dims(2);
    int xB = x/binning;
    int yB = y/binning;

    af::array binnedBlkSeries = moddims(qualityBinnedSeries(seq(xB, xB + binnedBlkSize -1), seq(yB, yB + binnedBlkSize -1), seq(0, nFrames-1)), dim4(binnedBlkSize*binnedBlkSize, nFrames));
    af::array gradSum = flat(af::sum( binnedBlkSeries, 0 ));
    //af_print(gradSum);

    /// Sort the sum of the gradient-norm. Sorted Array is the array ordered in decreasing quality
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, gradSum, 0, false);

    //int bufferSpace = 8;
    int blkSizeL = blkSize + bufferSpace;
    int xL = x - bufferSpace/2;
    int yL = y - bufferSpace/2;


    unsigned int *bestIndices = sortIndices.host<unsigned int>();
    /// Maybe the slice function would be better;
    af::array refBlk = arfSeries(seq(xL, xL + blkSizeL -1), seq(yL, yL + blkSizeL -1), bestIndices[0]);
    /// The refBlk is the reference block. But it will be used as the template to search in the bigger block
    /// (bestSearchBlks below) of size blkSizeL that will be cropped afterwards.
    af::array selectedInds = sortIndices(seq(0, nBest-1));


    af::array bestSearchBlks = arfSeries(seq(x, x + blkSize -1), seq(y, y + blkSize -1), selectedInds);
    /// Now let's start aligning the template onto those bigger bloks.
    af::array SADs;
    matchTemplate2(SADs, refBlk, bestSearchBlks);
    af::array dxArr, dyArr;
    fetchTMatch2Shifts(SADs, dxArr, dyArr, bestSearchBlks.dims());
    /// Dimensions of dxArr and dyArr: [1, nBest]
//    qDebug("Block shift: ");
//    af_print(dxArr - bufferSpace/2);
//    af_print(dyArr - bufferSpace/2);

    /// We crop over the refBlk which itself must be cropped by its extra bufferSpace
    x += bufferSpace/2;
    y += bufferSpace/2;

    int64 *dxHost = dxArr.host<int64>();
    int64 *dyHost = dyArr.host<int64>();

    /// Populate the stack.

    stackedBlks = constant(0, blkSize, blkSize, nBest);

//    af::sync();
//    af::timer timerFor = af::timer::start();

    for (int i = 0; i < nBest; i++)
    {
//        af::array xCropi = xCropRange(span, i);
//        af::array yCropi = yCropRange(span, i);
//        stackedBlks(span, span, i) = arfSeries(xCropi, yCropi, bestIndices[i]);
        int x1 = x - dxHost[i] ;
        int x2 = x1 + blkSize - 1;
        int y1 = y - dyHost[i];
        int y2 = y1 + blkSize - 1;
        stackedBlks(span, span, i) = arfSeries(seq(x1, x2), seq(y1, y2), bestIndices[i]);
    }
    stackedBlks.eval();

    //af::sync();
    //double loopTime = af::timer::stop(timerFor);
    //qDebug("ProcessingGlobalStack1:: for loop = %f ms", loopTime * 1000.0);

    /// Here, we have a stack of nBest blocks, in descending quality.
    /// We can either average them, with arithmetic, sigma-clip, or median,
    /// or mask-average them with a quality mask.
}

void RProcessing::populateResultListWithAr(QList<RMat *> rMatImageList, af::array &canvas, QString title)
{

    float *hostArf = canvas.host<float>();
    cv::Mat matImage(canvas.dims(1), canvas.dims(0), CV_32F, hostArf);
    matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
    RMat *blockRMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
    blockRMat->setImageTitle(title);
    blockRMat->setSOLAR_R(rMatImageList.at(0)->getSOLAR_R());
    resultList << blockRMat;
    delete[] hostArf;

}

void RProcessing::matchTemplate2(af::array &res, af::array &searchImg, af::array &img)
{
    af::array onesAr = constant(1, img.dims());
    af::array A2 = convolve2(searchImg * searchImg, onesAr);

    af::array imgF = flip(flip(img, 0), 1);
    af::array convImgs = convolve2(searchImg, imgF);
    res = A2 - 2 * convImgs;
}

void RProcessing::matchTemplate3(af::array &res, af::array &A, af::array &k, af::array &arrOnes)
{

    af::array A2 = convolve2(A * A, arrOnes);

    af::array k1 = flip(flip(k, 0), 1);
    af::array AK = convolve2(A, k1);
    res = A2 - 2 * AK;
}


void RProcessing::findMinLoc(af::array &ar, int &dx, int &dy)
{
    /// Get the location of the minimum in each image or 2D array
    /// Works in 2D only
    float minVal = 0;
    unsigned minLoc = 0;
    af::min<float>(&minVal, &minLoc, ar);

    dx = minLoc % ar.dims(0);
    dy = minLoc / ar.dims(0);
}

void RProcessing::findMinLoc(af::array & a, af::array &dx, af::array &dy)
{
    /// Get the location of the minimum in each image or 2D array in a series.
    /// Works in 3D.
    af::array minVal, idx0;
    af::min(minVal, idx0, moddims(a, a.dims(0)*a.dims(1), a.dims(2)));
    dx = idx0 % a.dims(0);
    dy = idx0 / a.dims(0);
}

void RProcessing::fetchTMatch2Shifts(af::array & a, int &dx, int &dy, dim4 dims)
{
    findMinLoc(a, dx, dy);
    dx -= (dims[0]/2 -1);
    dy -= (dims[0]/2 -1);
}

void RProcessing::fetchTMatch2Shifts(af::array & a, af::array &dx, af::array &dy, dim4 dims)
{
    findMinLoc(a, dx, dy);
    dx -= (dims[0]/2 -1);
    dy -= (dims[0]/2 -1);
}

void RProcessing::fetchTMatchShiftsGfor(af::array & ar, af::array &dxAr, af::array &dyAr)
{
    float minVal = 0;
    unsigned minLoc = 0;

    gfor(seq i, ar.dims(2))
    {
        af::min<float>(&minVal, &minLoc, ar(span, span, i));
        dxAr(i) = minLoc % ar.dims(0);
        dyAr(i) = minLoc / ar.dims(0);
    }

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
{   /// X-correlate upon the results of the limb-based registration
    /// Uses output variable "limbFitResultList1" from solarLimbRegisterSeries();

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
{   /// Limb-fitting based on canny edge detection

    /// Find contours
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

bool RProcessing::wernerLimbFit(QList<RMat*> rMatImageList, bool smooth, int smoothSize)
{
    /// Check if data exist
    if (rMatImageList.isEmpty())
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
    for (int i = 0; i < rMatImageList.size(); i++)
    {

        cv::Mat matImage0 = rMatImageList.at(i)->matImage.clone(); //normalizeByThresh(rMatLightList.at(i)->matImage, rMatLightList.at(i)->getIntensityLow(), rMatLightList.at(i)->getIntensityHigh(), rMatLightList.at(i)->getNormalizeRange());
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
        raphFindLimb( matImage, &wernerPoints, numDots, smooth, smoothSize);
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

        /// 2nd pass
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
        float median = calcMedian(distances, 0.1);

        std::vector<cv::Point> newPoints;
        for (int j = 0; j < wernerPoints.n; j++)
        {
//            if ( abs(distances[j] - mean[0] ) <  stddev[0] )
            if ( abs(distances[j] - median ) <  stddev[0] )
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

        /// 3rd pass
        cv::Point2f circleCenter2(circleOut2.a, circleOut2.b);
        std::vector<float> distances2(cleanDataPoints.n);
        for (int j = 0; j < cleanDataPoints.n; j++)
        {
            cv::Point2f point( (float) cleanDataPoints.X[j], (float) cleanDataPoints.Y[j]);
            float distance = cv::norm(point-circleCenter2);
            distances2[j] = distance;
        }
        cv::meanStdDev(distances2, mean, stddev);
        median = calcMedian(distances2, 0.1);

        std::vector<cv::Point> newPoints2;
        for (int j = 0; j < cleanDataPoints.n; j++)
        {
            //if ( abs(distances2[j] - mean[0] ) <  stddev[0] )
            if ( abs(distances2[j] - median ) <  stddev[0] )
            {
                cv::Point2f point( (float) cleanDataPoints.X[j], (float) cleanDataPoints.Y[j]);
                newPoints2.push_back(point);

            }
        }
        Data cleanDataPoints2(newPoints2.size());
        Circle circleOut3;
        if (!newPoints2.empty())
        {
            for (int j = 0; j < newPoints2.size() ; j++)
            {
                cleanDataPoints2.X[j] = newPoints2.at(j).x;
                cleanDataPoints2.Y[j] = newPoints2.at(j).y;
            }

            /// 3rd pass at fitting
            circleOut3 = CircleFitByTaubin(cleanDataPoints2);
        }
        /// End of 3rd pass

        if (newPoints2.empty())
        {   /// if there were only 1 pass, store it.
            circleOutList << circleOut;
            centers.append(cv::Point2f((float) circleOut.a, (float) circleOut.b));
        }
        else
        {   /// Store the 3rd pass
            circleOutList << circleOut3;
            centers.append(cv::Point2f((float) circleOut3.a, (float) circleOut3.b));
        }


        /// Display results
        matImage0.convertTo(matImage, CV_16U);
        //emit resultSignal(matImage, false, instruments::USET, QString("HPF"));

        contoursMat = cv::Mat::zeros(matImage.rows, matImage.cols, CV_8U);
        matImage0.convertTo(contoursMat, CV_8U, 256.0f /  rMatLightList.at(i)->getNormalizeRange());

        cv::cvtColor(contoursMat, contoursMat, CV_GRAY2RGB);

        if (showContours)
        {
            /// Display points in contoursMat
            for (int j = 0; j < wernerPoints.n; j++)
            {
                /// 1st pass in green
                cv::Point point( wernerPoints.X[j], wernerPoints.Y[j]);
                cv::line(contoursMat, point, point, cv::Scalar(0, 255, 0), 8, 8, 0);
            }
            if (!newPoints2.empty())
            {
                for (int j = 0; j < cleanDataPoints2.n; j++)
                {
                    /// 2nd pass in red
                    cv::Point point( cleanDataPoints2.X[j], cleanDataPoints2.Y[j]);
                    cv::line(contoursMat, point, point, cv::Scalar(255, 0, 0), 6, 6, 0);
                }
            }
        }

        if (showLimb)
        {


            cv::Scalar green = cv::Scalar(0, 255, 0);
            cv::Scalar red = cv::Scalar(255, 0, 0);
            cv::Scalar orange = cv::Scalar(255, 128, 0);
            // draw the 1st fitted circle in green
            cv::circle(contoursMat, circleCenter, circleOut.r, green, 2, 8);

            if (!newPoints2.empty())
            {
                // draw the 2nd fitted circle in orange
                cv::Point2f circleCenter3(circleOut3.a, circleOut3.b);
                cv::circle(contoursMat, circleCenter3, circleOut3.r, orange, 2, 8);
            }

        }

        /// Pack the results showing contours of all the edges
        contoursRMat = new RMat(contoursMat.clone(), false);
        contoursRMat->setImageTitle(QString("werner Limb Detection: Image # %1").arg(i+1));
        contoursRMatList << contoursRMat;

    }

    limbFitPlot = new QCustomPlot();
    limbFitPlot->addGraph();

    // Prepare the plot data
    QVector<double> frameNumbers;
    QVector<double> radius;
    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        frameNumbers << i;
        radius << circleOutList.at(i).r;
        qDebug("radius = %f", circleOutList.at(i).r);
    }
    std::vector<double> radii = radius.toStdVector();
    cv::Scalar tempRadius = cv::mean(radii);
    meanRadius = (float) tempRadius[0];
    qDebug() << "vector radius = " << radius;
    qDebug("meanRadius() = %f", meanRadius);

    limbFitPlot->graph(0)->setData(frameNumbers, radius);
    limbFitPlot->rescaleAxes();
    limbFitPlot->xAxis->setRange(0, rMatLightList.size());


    return true;
}

bool RProcessing::solarLimbRegisterSeries(QList<RMat*> rMatImageList)
{
    /// Align the image series based on the chosen limb fitting algorithm and the values
    /// in this->centers


    if (centers.isEmpty())
    {
        tempMessageSignal(QString("Run limb fitting first."));
        return false;
    }

    if (!limbFitResultList1.isEmpty())
    {
        qDeleteAll(limbFitResultList1);
        limbFitResultList1.clear();
    }



    for (int i = 0 ; i < rMatImageList.size() ; ++i)
    {
        qDebug("Canny-registering image # %i ", i+1);
        /// Register series
        limbFitWarpMat = cv::Mat::eye( 2, 3, CV_32FC1 );
        cv::Point2f origin(rMatImageList.at(i)->matImage.cols / 2.0f, rMatImageList.at(i)->matImage.rows / 2.0f);

        cv::Point2f delta = centers.at(i) - origin;
        limbFitWarpMat.at<float>(0, 2) = delta.x;
        limbFitWarpMat.at<float>(1, 2) = delta.y;
        std::cout << "limbFitWarpMat = " << std::endl << " " << limbFitWarpMat << std::endl << std::endl;

        cv::Mat registeredMat = rMatImageList.at(i)->matImage.clone();
        registeredMat.convertTo(registeredMat, CV_32F);

        cv::warpAffine(registeredMat, registeredMat, limbFitWarpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);

        registeredMat.convertTo(registeredMat, CV_16U);

        RMat *resultMat = new RMat(registeredMat, false, rMatImageList.at(i)->getInstrument());
        resultMat->setImageTitle(QString("Registered image # ") + QString::number(i));
        resultMat->setDate_time(rMatImageList.at(i)->getDate_time());
        limbFitResultList1 << resultMat;
    }


    return true;
}



void RProcessing::raphFindLimb(cv::Mat matImage, Data *dat, int numDots, bool smooth, int smoothSize)
{
    /// Fix USET images for bad 1st column and 1st row
    /// This fix won't be needed after using Emil's updated version
    /// of Suncap, but still needed for the older files in the USET database.
    fixUset(matImage);
    cv::Mat matSlice;
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
    cv::Size kernelSize(smoothSize, smoothSize);
    if (smoothSize == 0) { smooth = false; }

    for (int ii = 0; ii < numDots; ii++)
    {
        int X = naxis1/4 + ii*naxis1/(2*numDots);
        int Y = naxis2/4 + ii*naxis2/(2*numDots);
        /// Left-hand slices (no copy)
        matSlice = matImage(cv::Range(Y, Y+1), cv::Range(0, naxis1/4));
        if (smooth){ cv::blur(matSlice, matSlice, kernelSize); }
        cv::filter2D(matSlice, gradXLeft, -1, kernelXLeft, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradXLeft, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii] = maxLoc.x;
        dat->Y[ii] = Y;
        /// Right-hand slices (no copy)
        matSlice = matImage(cv::Range(Y, Y+1), cv::Range(3*naxis1/4, naxis1));
        if (smooth){ cv::blur(matSlice, matSlice, kernelSize); }
        cv::filter2D(matSlice, gradXRight, -1, kernelXRight, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradXRight, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + numDots] = maxLoc.x + 3*naxis1/4;
        dat->Y[ii + numDots] = Y;
        ///Bottom slices (no copy)
        matSlice = matImage(cv::Range(0, naxis2/4), cv::Range(X, X+1));
        if (smooth){ cv::blur(matSlice, matSlice, kernelSize); }
        cv::filter2D(matSlice, gradYBottom, -1, kernelYBottom, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        cv::minMaxLoc(gradYBottom, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + 2*numDots] = X;
        dat->Y[ii + 2*numDots] = maxLoc.y;
        ///Top slices (no copy)
        matSlice = matImage(cv::Range(3*naxis2/4, naxis2), cv::Range(X, X+1));
        if (smooth){ cv::blur(matSlice, matSlice, kernelSize); }
        cv::filter2D(matSlice, gradYTop, -1, kernelYTop, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
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
    normalizeRMat->setDate_time(rMat->getDate_time());
    normalizeRMat->setSOLAR_R(rMat->getSOLAR_R());
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

    normalizedMatImage = normalizedMatImage * alpha + beta;
    cv::threshold(normalizedMatImage, normalizedMatImage, newRange - 1, newRange - 1, cv::THRESH_TRUNC);

    if (matImage.type() == CV_32F || matImage.type() == CV_32FC3)
    {
        normalizedMatImage.convertTo(normalizedMatImage, CV_16U);
    }
    else
    {
        normalizedMatImage.convertTo(normalizedMatImage, matImage.type());
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
    if (matImage.type() == CV_8U)
    {
        for (int x = 0; x < matImage.cols; x++)
        {   /// 1st rows
            matImage.at<uchar>(0, x) = matImage.at<uchar>(4, x);
            matImage.at<uchar>(1, x) = matImage.at<uchar>(4, x);
            matImage.at<uchar>(2, x) = matImage.at<uchar>(4, x);
            matImage.at<uchar>(3, x) = matImage.at<uchar>(4, x);
            /// 1st column
        }
        for (int y = 0; y < matImage.rows; y++)
        {   /// 1st rows
            matImage.at<uchar>(y, 0) = matImage.at<uchar>(y, 4);
            matImage.at<uchar>(y, 1) = matImage.at<uchar>(y, 4);
            matImage.at<uchar>(y, 2) = matImage.at<uchar>(y, 4);
            matImage.at<uchar>(y, 3) = matImage.at<uchar>(y, 4);
            /// 1st column
        }
    }
    else if (matImage.type() == CV_16U)
    {
        for (int x = 0; x < matImage.cols; x++)
        {   /// 1st rows
            matImage.at<ushort>(0, x) = matImage.at<ushort>(4, x);
            matImage.at<ushort>(1, x) = matImage.at<ushort>(4, x);
            matImage.at<ushort>(2, x) = matImage.at<ushort>(4, x);
            matImage.at<ushort>(3, x) = matImage.at<ushort>(4, x);
        }
        for (int y = 0; y < matImage.rows; y++)
        {   /// 1st rows
            matImage.at<ushort>(y, 0) = matImage.at<ushort>(y, 4);
            matImage.at<ushort>(y, 1) = matImage.at<ushort>(y, 4);
            matImage.at<ushort>(y, 2) = matImage.at<ushort>(y, 4);
            matImage.at<ushort>(y, 3) = matImage.at<ushort>(y, 4);
        }
    }

    else if (matImage.type() == CV_32F)
    {
        for (int x = 0; x < matImage.cols; x++)
        {   /// 1st rows
            matImage.at<float>(0, x) = matImage.at<float>(4, x);
            matImage.at<float>(1, x) = matImage.at<float>(4, x);
            matImage.at<float>(2, x) = matImage.at<float>(4, x);
            matImage.at<float>(3, x) = matImage.at<float>(4, x);
        }
        for (int y = 0; y < matImage.rows; y++)
        {   /// 1st rows
            matImage.at<float>(y, 0) = matImage.at<float>(y, 4);
            matImage.at<float>(y, 1) = matImage.at<float>(y, 4);
            matImage.at<float>(y, 2) = matImage.at<float>(y, 4);
            matImage.at<float>(y, 3) = matImage.at<float>(y, 4);
        }
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

void RProcessing::setSharpenLiveStatus(bool status)
{
    this->sharpenLiveStatus = status;
}

void RProcessing::setStackWithMean(bool status)
{
    this->stackWithMean = status;
}

void RProcessing::setStackWithSigmaClip(bool status)
{
    this->stackWithSigmaClip = status;
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

QList<RMat *> RProcessing::getLuckyBlkList()
{
    return luckyBlkList;
}

QVector<Circle> RProcessing::getCircleOutList()
{
    return circleOutList;
}

float RProcessing::getMeanRadius()
{
    return meanRadius;
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

void RProcessing::red_tab(int* red, int* green ,int* blue)
{
    int const redpr[] =
    {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,   0,   1,   2,   5,   7,  10,  11,  13,  15,  17,  18,  21,  23,  24,  27,
       28,  30,  33,  34,  36,  37,  40,  42,  43,  46,  47,  49,  50,  53,  55,  56,
       59,  60,  62,  63,  66,  68,  69,  70,  73,  75,  76,  78,  81,  82,  84,  85,
       88,  89,  91,  92,  95,  97,  98,  99, 102, 104, 105, 107, 108, 111, 113, 114,
       115, 118, 120, 121, 123, 126, 127, 128, 130, 131, 134, 136, 137, 139, 141, 143,
       144, 146, 147, 150, 152, 153, 155, 156, 159, 160, 162, 163, 166, 168, 169, 170,
       172, 175, 176, 178, 179, 181, 184, 185, 186, 188, 189, 192, 194, 195, 197, 198,
       201, 202, 204, 205, 207, 210, 211, 212, 214, 215, 218, 220, 221, 223, 224, 227,
       228, 230, 231, 233, 236, 237, 239, 240, 241, 243, 246, 247, 249, 250, 252, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
       255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
    int const greenpr[] =
    { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  3,  5,  7,  9, 13, 15, 17, 18, 20, 24,
      26, 28, 30, 32, 35, 37, 39, 41, 43, 47, 49, 51, 52, 54, 58, 60, 62, 64, 66, 69,
      71, 73, 75, 77, 81, 83, 85, 86, 88, 90, 94, 96, 98,100,102,105,107,109,111,113,
      115,119,120,122,124,126,130,132,134,136,137,139,143,145,147,149,151,154,156,158,
      160,162,164,168,170,171,173,175,177,181,183,185,187,188,192,194,196,198,200,202,
      205,207,209,211,213,215,219,221,222,224,226,228,232,234,236,238,239,241,245,247,
      249,251,253,255,255,255,255,255,255,255,255,255,255,255,255,255};
    int const bluepr[] =
    {  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
       0,  0,  0,  0,  0,  0,  0,  0,  3,  7, 11, 15, 23, 27, 31, 35, 39, 47, 51, 54,
       58, 62, 66, 74, 78, 82, 86, 90, 94,102,105,109,113,117,125,129,133,137,141,145,
       153,156,160,164,168,172,180,184,188,192,196,200,207,211,215,219,223,227,235,239,
       243,247,251,255,255,255,255,255,255,255,255,255,255,255,255,255};

    for (int i=0; i < 256; i++)
    {
        red[i] = redpr[i];
        blue[i] = bluepr[i];
        green[i] = greenpr[i];
    }

}

cv::Mat RProcessing::scalePreviewImage(float sunX,float sunY,float sunR, cv::Mat matImage, char filter)
{
    /* scale intensities for jpeg output */

    int naxis1 = matImage.cols;
    int naxis2 = matImage.rows;

    int xsq, ysq, x,y,maxOut=0, offset;
    long meanOut = 0;
    long meanIn = 0;
    double diskMean;
    float scalein, scaleout, rq;
    float r2 = 2;
    float Pi = 3.1415;
    int z;


    /// Squared oversized solar radius
    rq = powf(sunR + r2, 2);
    if (filter=='H')
    {
//#pragma omp parallel for
        for (x=0;x<naxis1;x++)
        {
            /// Squared x-distance to disk center
            xsq = std::pow(x-sunX,2);
            for (y=0;y<naxis2;y++)
            {
                /// Squared y-distance to disk center
                ysq = std::pow(y-sunY,2);
                //index = x+y*naxis1;
                if (xsq + ysq >= rq)
                {
                    /// Maximum value outside the disk
                    z = (int) matImage.at<ushort>(y, x);
                    if (z > maxOut)
                    {
                       maxOut =  z;
                    }
                    /// Sum intensity of pixels outside the disk (to get the mean)
                    meanOut +=  z;
                }
               meanIn += z;
            }
        }
        /// Approx. mean intensity of pixels outside the disk
        meanOut /= (naxis1*naxis2-rq*Pi);
        meanOut *= 1.2;
        meanIn /= sunR*sunR*Pi;
    }

    scaleout = (maxOut-meanOut);

    diskMean  = (cv::sum(matImage)[0]) / (sunR*sunR*Pi);


    /// Prepare output image
    cv::Mat outputMat(naxis2, naxis1, CV_8U);
    /// Clip image between 10 and 250?

    if (filter == 'H')
    {
       /// Scaling factor so that image within solar disk have a mean value of 170
       scalein = 170.0/diskMean;
       offset = 5;
    }
    else if (filter == 'P')
    {
        scalein = 165.0/diskMean;
        offset=0;
    }
    else if (filter == 'C')
    {
        scalein = 115.0/diskMean;
        offset=-10;
    }
    else
    {
        qDebug("Image filter type is not recognized. Returning.");
        return outputMat;
    }


    for (x=0;x<naxis1;x++)
    {
        xsq = std::pow(x-sunX,2);
        for (y=0;y<naxis2;y++)
        {
            ysq = std::pow(y-sunY,2);
            z = (int) matImage.at<ushort>(y, x);
            if (xsq + ysq <= rq )
            {
                int z0 = z;
                z = z*scalein-offset;
//                qDebug("scaling = %f", scalein);
//                qDebug("offset = %d", offset);
//                qDebug("before scaling z0 = %d", z0);
//                qDebug("after scaling z = %d,", z);
                if(z < 10)
                {
                    z = 10;
                }
                else if(z > 250)
                {
                    z = 250;
                }
            }
            else
            {
                /// scale corona for H-alpha
                if (filter=='H'){ z = 180*(z-meanOut)/scaleout; }
                else {z = 0;}
                if(z < 0){ z = 0; }
                if( z > 220){ z = 220; }
            }
            outputMat.at<uchar>(y, x) = (uchar) z;
        }
    }
    return outputMat;
}

void RProcessing::printArDims(af::array &ar)
{
    qDebug("ar.dims() = [%d, %d, %d]", ar.dims(0), ar.dims(1), ar.dims(2));
}

cv::Mat RProcessing::wSolarColorize(cv::Mat matImage, char filter)
{
    int red[256];
    int green[256];
    int blue[256];

    red_tab(red, green, blue);

    int naxis1 = matImage.cols;
    int naxis2 = matImage.rows;
    float sunR = meanRadius;

    cv::Mat mat8Bit(naxis2, naxis1, CV_8U);
    if (filter == 'H')
    {
        mat8Bit = scalePreviewImage(naxis1/2, naxis2/2, sunR, matImage, 'H');
    }
    else
    {
        mat8Bit = cv::Mat::zeros(naxis2, naxis1, CV_8UC3);
        return mat8Bit;
    }

    double dataMax;
    cv::minMaxLoc(mat8Bit, NULL, &dataMax);
    int fac = 255/ dataMax;

    cv::Mat coloredImg = cv::Mat::zeros(naxis2, naxis1, CV_8UC3);

    for (uint x = 0 ; x < naxis1 ; x++)
    {
        for (uint y = 0 ; y < naxis2 ; y++)
        {
            uchar z = mat8Bit.at<uchar>(y,x);
            coloredImg.at<cv::Vec3b>(y,x) = fac*cv::Vec3b(red[z], green[z], blue[z]);
        }
    }

    return coloredImg;
}

QList<RMat*> RProcessing::wSolarColorizeSeries(QList<RMat *> rMatImageList, char filter)
{
    QList<RMat*> rMat8BitList;

    for (int i = 0 ; i < rMatImageList.size() ; ++i)
    {
        cv::Mat matImage8Bit = wSolarColorize(rMatImageList.at(i)->matImage, filter);

        RMat *rMat8Bit = new RMat(matImage8Bit, false);
        rMat8Bit->setImageTitle(QString("8-bit image # %1").arg(i));
        rMat8Bit->setDate_time(rMatImageList.at(i)->getDate_time());
        rMat8BitList << rMat8Bit;
    }

    return rMat8BitList;
}



//void ushrpMask(int naxis1,int naxis2,int* data,int Datamin,int Datamax){
//   /* unsharp masking of image 2*image - smoothed image*/
//   int i,x,u,area,d=3,dd;
//   float* data2;
//   data2 = (float*) malloc(naxis1*naxis2*sizeof(float));
//   float* data3;
//   data3 = (float*) malloc(naxis1*naxis2*sizeof(float));
//   area = naxis1*naxis2;
//   dd = (2*d+1)*(2*d+1);
//   /* split into first half of image and second half -> only one "if"
//      smooth horizontally */
//   for (i=0;i<area/2;i++){
//      data2[i] = 0;
//      for (x=-d;x<=d;x++){
//         u = i+x;
//         if(u<0){u+=area;}
//         data2[i] += data[u];
//         }
//      data3[i] = data2[i];
//      }
//   for (i=area/2;i<area;i++){
//      data2[i] = 0;
//      for (x=-d;x<=d;x++){
//         u = i+x;
//         if(u>=area){u-=area;}
//         data2[i] += data[u];
//         }
//      data3[i] = data2[i];
//      }
//  /* split into first half of image and second half -> only one "if"
//     smooth vertically */
//   for (i=0;i<area/2;i++){
//     for (x=-d;x<=d;x++){
//         u = i+x*naxis1;
//         if(u<0){u+=area;}
//         data3[i] += data2[u];
//         }
//      }
//   for (i=area/2;i<area;i++){
//     for (x=-d;x<=d;x++){
//         u = i+x*naxis1;
//         if(u>=area){u-=area;}
//         data3[i] += data2[u];
//         }
//      }
//   #pragma omp parallel for
//   for (i=0;i<area;i++){
//      data[i] = 2.0*data[i]-data3[i]/dd;
//      if (data[i]<Datamin){data[i]=Datamin;}
//      if (data[i]>Datamax){data[i]=Datamax;}}
//   free(data2);
//   free(data3);
//}



