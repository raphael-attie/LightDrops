#include "rprocessing.h"

#include <QFileDialog>

// cifitsio
#include <fitsio.h>

//opencv
#include <opencv2/world.hpp>

#include "imagemanager.h"
#include "parallelcalibration.h"
#include "typedefs.h"

RProcessing::RProcessing(QObject *parent): QObject(parent),
    masterBias(NULL), masterDark(NULL), masterFlat(NULL), masterFlatN(NULL), stackedRMat(NULL), useUrlsFromTreeWidget(false), useXCorr(false),
    masterWithMean(true), masterWithSigmaClip(false), stackWithMean(true), stackWithSigmaClip(false), radius(0), radius1(0), radius2(0), radius3(0), meanRadius(0),
    useROI(false), maskCircleX(0), maskCircleY(0), maskCircleRadius(0), blkSize(32), binning(2)
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
    /// Load list of images as lights
    listImageManager->loadData(urls);
    rMatLightList = listImageManager->getRMatImageList();
}

void RProcessing::loadRMatBiasList(QList<QUrl> urls)
{
    /// Load list of images as Bias images
    listImageManager->loadData(urls);
    rMatBiasList = listImageManager->getRMatImageList();
}

void RProcessing::loadRMatDarkList(QList<QUrl> urls)
{
    /// Load list of images as Dark images
    listImageManager->loadData(urls);
    rMatDarkList = listImageManager->getRMatImageList();
}

void RProcessing::loadRMatFlatList(QList<QUrl> urls)
{
    /// Load list of images as Flat-field images
    listImageManager->loadData(urls);
    rMatFlatList = listImageManager->getRMatImageList();
}



void RProcessing::exportMastersToFits()
{
    if (masterBias != NULL && !masterBias->matImage.empty())
    {
        QFileInfo fileInfo(treeWidget->getBiasDir().filePath(QString("masterBias.fits")));
        QString masterBiasPath = setupFileName(fileInfo);

        exportToFits(masterBias, masterBiasPath);
        std::cout << "Bias dir: " << treeWidget->getBiasDir().rootPath().toStdString() << std::endl;
        std::cout << "masterBias exported at: " << masterBiasPath.toStdString() << std::endl;
    }

    if (masterDark != NULL && !masterDark->matImage.empty())
    {
        QFileInfo fileInfo(treeWidget->getDarkDir().filePath(QString("masterDark.fits")));
        QString masterDarkPath = setupFileName(fileInfo);

        exportToFits(masterDark, masterDarkPath);
        std::cout << "masterDark exported at: " << masterDarkPath.toStdString() << std::endl;

    }

    if (masterFlat != NULL && !masterFlat->matImage.empty())
    {

        QFileInfo fileInfo(treeWidget->getFlatDir().filePath(QString("masterFlat.fits")));
        QString masterFlatPath = setupFileName(fileInfo);

        exportToFits(masterFlat, masterFlatPath);
        std::cout << "masterFlat exported at: " << masterFlatPath.toStdString() << std::endl;
    }

    tempMessageSignal(QString("Exported master calibration frames"), 0);
}

void RProcessing::exportFramesToFits(QList<RMat *> rMatImageList, QDir exportDir, bool useBasename)
{

    for (int i = 0 ; i < rMatImageList.size() ; i++)
    {
        QString indexQStr = QString("%1").arg(i, 5, 10, QChar('0'));
        QString fileName(QString("image_") + indexQStr + QString(".fits"));
        if (useBasename)
        {
            fileName = rMatImageList.at(i)->getFileInfo().baseName() + QString("_C.fits");
        }
        QFileInfo fileInfo(exportDir.filePath(fileName));
        QString filePath = setupFileName(fileInfo);

        exportToFits(rMatImageList.at(i), filePath);
    }
}

QString RProcessing::setupFileName(QFileInfo fileInfo)
{
    QString filePath = fileInfo.filePath();

    uint fileNumber = 1;
    QFileInfo fileInfoTest(filePath);
    while (fileInfoTest.exists())
    {
        QString baseName = fileInfo.baseName() + QString("_") + QString::number(fileNumber) +QString(".");
        filePath = fileInfo.absoluteDir().filePath(baseName + fileInfoTest.suffix());
        fileInfoTest = QFileInfo(filePath);
        fileNumber++;
        qDebug() << "RProcessing::setupFileName():: filePath =" << filePath;
    }

    return filePath;
}


void RProcessing::exportToFits(RMat *rMatImage, QString QStrFilename)
{

    if (rMatImage->matImage.channels() == 3)
    {
        emit tempMessageSignal(QString("Multi-channel image. Select another image type (e.g Tiff)"));
        return;
    }
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

    /// UPDATE: I realize a bit late that I made a very stupid mistake with this, because of the flat-fielding.
    /// The flat-fielding change the image to floating points with a completely different range of values.
    /// Converting to integers without any other modification will totally mess up the image.
    /// I could stick to floating points in FITS but ultimately I need to convert to TIFF so must restore a range of integers
    /// That range depends on the instrument. 14 bits, 16 bits, etc...
    /// If the image is stacked, it gains precision. If it is e.g originally a 14 bit image, then stretching the image linearly
    /// so it fits within 16-bit/channel should be reasonnable.
    /// But then... when do I do it. What values do I choose when stretching?...
    /// I need revisit the flat-fielding operation to make things more rigourous.
    /// UPDATE2: The flat fielding makes much more sense now, as it gives me much more natural colors without
    /// touching any of the r,g,b white balance values. And the range is restored.
    if (rMatImage->isBayer())
    {
        cv::Mat tempImage16;
        rMatImage->matImage.convertTo(tempImage16, CV_16U);
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)tempImage16.data, &status);
    }
    else if (rMatImage->matImage.type() == CV_16UC1)
    {
        fits_create_img(fptr, USHORT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TUSHORT, fpixel, nPixels, (ushort*)rMatImage->matImage.data, &status);
    }
    else if (rMatImage->matImage.type() == CV_32FC1)
    {
        //  Create the primary array image (32-bit float  pixels)
        fits_create_img(fptr, FLOAT_IMG, naxis, naxes, &status);
        // Write the array data into the file.
        fits_write_img(fptr, TFLOAT, fpixel, nPixels, (float*)rMatImage->matImage.data, &status);
    }
    else
    {
        emit tempMessageSignal(QString("Image type not recognized"));
        return;
    }

    // Write BAYER keyword
    fits_write_key(fptr, TLOGICAL, keyname, &bayer, NULL, &status);
    // Write flipUD keyword
    fits_write_key(fptr, TLOGICAL, "FLIPUD", &rMatImage->flipUD, NULL, &status);

    if (rMatImage->getInstrument() == instruments::USET)
    {
        char keyNameTelescop[] = "TELESCOP";
        char keyValueTelescop[] = "USET";
        fits_write_key(fptr, TSTRING, keyNameTelescop, keyValueTelescop, NULL, &status);
        char keyNameRadius[] = "SOLAR_R";
        float keyValueRadius = rMatImage->getSOLAR_R();
        fits_write_key(fptr, TFLOAT, keyNameRadius, &keyValueRadius, NULL, &status);

    }
    else if (rMatImage->getInstrument() == instruments::DSLR)
    {
        char keyNameTelescop[] = "TELESCOP";
        char keyValueTelescop[] = "DSLR";
        fits_write_key(fptr, TSTRING, keyNameTelescop, keyValueTelescop, NULL, &status);

        char keyNameTEMP[] = "TEMP";
        float keyValueTEMP = rMatImage->getTEMP();
        fits_write_key(fptr, TFLOAT, keyNameTEMP, &keyValueTEMP, NULL, &status);
    }
    else
    {
        char keyNameTelescop[] = "TELESCOP";
        char keyValueTelescop[] = "UNDEFINED";
        fits_write_key(fptr, TSTRING, keyNameTelescop, keyValueTelescop, NULL, &status);
    }

    char keyNameXPOSURE[] = "XPOSURE";
    float keyValueXPOSURE = rMatImage->getXPOSURE();
    fits_write_key(fptr, TFLOAT, keyNameXPOSURE, &keyValueXPOSURE, NULL, &status);
    char keyNameEXPTIME[] = "EXPTIME";
    float keyValueEXPTIME = rMatImage->getExpTime();
    fits_write_key(fptr, TFLOAT, keyNameEXPTIME, &keyValueEXPTIME, NULL, &status);


    // Close the file
    fits_close_file(fptr, &status);

}

void RProcessing::batchExportToFits(QList<QUrl> urls, QString exportDir)
{
    if (urls.empty())
    {
        return;
    }

    bool useInputDirectory = false;
    exportQDir = QDir(exportDir);

    if (exportDir.isEmpty())
    {
        useInputDirectory = true;
    }



    // Do a batch export to FITS and keep the same basenames
    QString format("fits");
    for(int i = 0; i < urls.size(); i++)
    {

        QUrl url = urls.at(i);
        // Import file
        ImageManager iManager(url);

        if (useInputDirectory)
        {
            QFileInfo fileInfo(urls.at(i).toLocalFile());
            QString fileDir = fileInfo.path();
            exportQDir = QDir(fileDir);
        }

        QFileInfo fileInfo1(url.fileName());
        QString basename = fileInfo1.baseName();
        QString fileName(basename + QString(".") + format);
        QFileInfo fileInfo2(exportQDir.filePath(fileName));
        // Check if file does not really exist, rename if that is the cases
        QString filePath = setupFileName(fileInfo2);

        // export
        exportToFits(iManager.getRMatImage(), filePath);
    }
}

cv::Mat RProcessing::rescaleForExport8Bits(cv::Mat matImage, float alpha, float beta)
{
    int channels = matImage.channels();
    float alpha2 = (float) (255.0f * alpha);
    float beta2 = (float) (255.0f * beta);
    cv::Mat tempMat = matImage.clone();
    tempMat.convertTo(tempMat, CV_8UC(channels), alpha2, beta2);

    return tempMat;
}

void RProcessing::exportToTiff(RMat *rMatImage, QString QStrFilename)
{
    QString tiffFilename = QStrFilename + QString(".tiff");
    std::string strFilename(tiffFilename.toStdString());
    cv::Mat matImage16;
    rMatImage->matImage.convertTo(matImage16, CV_16U);
    if (rMatImage->matImage.channels() == 3)
    {
        cv::cvtColor(matImage16, matImage16, CV_RGB2BGR);
    }
    else
    {
        cv::cvtColor(matImage16, matImage16, CV_GRAY2BGR);
    }

    try {
            cv::imwrite(strFilename, matImage16);
        }
        catch (runtime_error& ex) {
            fprintf(stderr, "Exception converting image to Tiff: %s\n", ex.what());
            return;
    }
}

void RProcessing::exportToJpeg(RMat *rMatImage, QString QStrFilename)
{
    QString jpegFilename = QStrFilename + QString(".jpg");
    std::string strFilename(jpegFilename.toStdString());
    cv::Mat matImage16;

    if (rMatImage->matImage.channels() == 3)
    {
        rMatImage->matImage.convertTo(matImage16, CV_16U);
        cv::cvtColor(matImage16, matImage16, CV_RGB2BGR);
        std::cout << "exportToJpeg:: cv::cvtColor(matImage16, matImage16, CV_RGB2BGR)" << std::endl;
    }
    else
    {
        rMatImage->matImageRGB.convertTo(matImage16, CV_16U);
        cv::cvtColor(matImage16, matImage16, CV_RGB2BGR);
        std::cout << "exportToJpeg:: cv::cvtColor(matImage16, matImage16, CV_GRAY2BGR)" << std::endl;
    }
    float alpha = currentROpenGLWidget->getAlpha();
    float beta = currentROpenGLWidget->getBeta();
    cv::Mat exportMat = rescaleForExport8Bits(matImage16, alpha, beta);

    try {
            cv::imwrite(strFilename, exportMat);
        }
        catch (runtime_error& ex) {
            fprintf(stderr, "Exception converting image to jpeg: %s\n", ex.what());
            return;
    }

}

// Make a function to write image

void RProcessing::loadMasterBias()
{
    /// Used off-screen because no image is loaded beforehand in the ROpenGLWidget.
    /// So the url must exist in the treeWidget
    if (treeWidget->getBiasUrls().size() == 1)
    {
        masterBiasUrl = treeWidget->getBiasUrls().at(0);
        std::cout << "RProcessing::loadMasterBias() masterBiasUrl = " << masterBiasUrl.toLocalFile().toStdString() << std::endl;
    }
    else if (treeWidget->getBiasUrls().empty())
    {

        qDebug("You need at least one Bias image in the calibration tree");
        tempMessageSignal(QString("You need at least one Bias image in the calibration tree"));
        return;
    }
    else
    {
        qDebug("Master Bias unknown");
        tempMessageSignal(QString("master Bias unknown"));
        return;
    }

    /// Let the ImageManager on the stack, (so we don't have to call delete).
    ///  and make a deep copy of the data in RMat before leaving this function.
    ImageManager imageManager(masterBiasUrl);

    masterBias = new RMat(*imageManager.getRMatImage());
    if (masterBias->matImage.type() != CV_32F)
    {
        masterBias->matImage.convertTo(masterBias->matImage, CV_32F);
    }
}

void RProcessing::loadMasterDark()
{
    /// Used off-screen because no image is loaded beforehand in the ROpenGLWidget.
    /// So the url must exist in the treeWidget
    if (treeWidget->getDarkUrls().size() == 1)
    {
        masterDarkUrl = treeWidget->getDarkUrls().at(0);
    }
    else if (treeWidget->getDarkUrls().empty())
    {
        return;
    }
    else
    {
        qDebug("Master Dark unknown");
        tempMessageSignal(QString("master Dark unknown"));
    }

    /// Let the ImageManager on the stack, (so we don't have to call delete).
    ///  and make a deep copy of the data in RMat before leaving this function.
    ImageManager imageManager(masterDarkUrl);

    masterDark = new RMat(*imageManager.getRMatImage());
    if (masterDark->matImage.type() != CV_32F)
    {
        masterDark->matImage.convertTo(masterDark->matImage, CV_32F);
    }

}

void RProcessing::loadMasterFlat()
{
    /// Used off-screen because no image is loaded beforehand in the ROpenGLWidget.
    /// So the url must exist in the treeWidget
    if (treeWidget->getFlatUrls().size() == 1)
    {
        masterFlatUrl = treeWidget->getFlatUrls().at(0);
    }
    else if (treeWidget->getFlatUrls().empty())
    {
        return;
    }
    else
    {
        qDebug("Master Flat unknown");
        tempMessageSignal(QString("master Flat unknown"));
    }

    /// Let the ImageManager on the stack, (so we don't have to call delete).
    ///  and make a deep copy of the data in RMat before leaving this function.
    ImageManager imageManager(masterFlatUrl);

    masterFlat = new RMat(*imageManager.getRMatImage());

    if (masterFlat->matImage.type() != CV_32F)
    {
        masterFlat->matImage.convertTo(masterFlat->matImage, CV_32F);
    }
}

void RProcessing::showMinMax(const cv::Mat &matImage)
{
    std::cout<< "Show Min/Max" << std::endl;
    double min, max;
    cv::minMaxLoc(matImage, &min, &max);
    qDebug("[Min , Max] = [%f , %f]", min, max);
    emit tempMessageSignal(QString("[Min , Max] = [%1 , %2]").arg(min, 0, 'f', 2).arg(max, 0, 'f', 2), 0);
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
    float medianVal = 0;

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

af::array RProcessing::rMatMul(const af::array &lhs, const af::array &rhs)
{
    return lhs * rhs;
}

void RProcessing::rebin(const af::array &a, af::array &b, int binning)
{
    if (binning == 1){b = a;}
    if (binning == 2){rebin2(a, b);}
    if (binning == 4){rebin4(a, b);}
}

void RProcessing::rebin2(const af::array &a, af::array &b2)
{
    int binning = 2;
    af::array b1 = a(af::seq(0, af::end, binning), af::seq(0, af::end)) + a(af::seq(1, af::end, binning), af::seq(0, af::end));
    b2  = b1(af::seq(0, af::end), af::seq(0, af::end, binning)) + b1(af::seq(0, af::end), af::seq(1, af::end, binning));
}

void RProcessing::rebin4(const af::array &a, af::array &b4)
{
    af::array b2;
    rebin2(a, b2);
    rebin2(b2, b4);
}

void RProcessing::rebin(RMat *rMatImage, RMat *binnedRMatImage, int binning)
{
    int naxis1 = rMatImage->matImage.rows;
    int naxis2 = rMatImage->matImage.cols;
    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    rMatImage->matImage.convertTo(tempMat, CV_32F);
    af::array a(naxis2, naxis1, (float*) tempMat.data);
    af::array b;
    rebin(a, b, binning);

    float *bH = b.host<float>();
    cv::Mat matImage(naxis2/binning, naxis1/binning, CV_32F, bH);
    matImage.convertTo(matImage, rMatImage->matImage.type());
    binnedRMatImage->matImage = matImage;
    //binnedRMatImage->setInstrument(rMatImage->getInstrument());
    binnedRMatImage->setImageTitle(QString("binned"));
    binnedRMatImage->calcStats();
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

    stackedRMat->setInstrument(rMatImageList.at(0)->getInstrument());
    stackedRMat->setSOLAR_R(rMatImageList.at(0)->getSOLAR_R());
    emit resultSignal(stackedRMat);

    return;
}

RMat* RProcessing::average(QList<RMat*> rMatList)
{
    if (rMatList.size() == 1)
    {
        return rMatList.at(0);
    }
    /// Averages a series of cv::Mat images using arithmetic mean.
    int naxis1 = rMatList.at(0)->matImage.cols;
    int naxis2 = rMatList.at(0)->matImage.rows;

    /// Need to create a Mat image that will host the result of the average.
    /// It needs to be the same type, and has the number of channels as the Mat images of the series.
    cv::Mat avgImg = cv::Mat::zeros(naxis2, naxis1, rMatList.at(0)->matImage.type());
    avgImg.convertTo(avgImg, CV_32F);

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
    RMat *rMatAvg = new RMat(avgImg, rMatList.at(0)->isBayer(), rMatList.at(0)->getInstrument(), rMatList.at(0)->getXPOSURE(), rMatList.at(0)->getTEMP());
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
//    float sigmaFactor = 1.0f;
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
        rMatSharpList.at(i)->setInstrument(rMatImageList.at(i)->getInstrument());
        rMatSharpList.at(i)->setSOLAR_R(rMatImageList.at(i)->getSOLAR_R());

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

    af::setBackend(AF_BACKEND_OPENCL);

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
    kernel(af::seq(3,6), af::seq(3,6)) = 1;
    af_print(kernel);

    /// Store the matImageList (converted to float) in the GPU array, and get a binned version
    ///

    af::timer afTimer1 = af::timer::start();

    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempArf(naxis2, naxis1, (float*) tempMat.data);
        arfSeries(af::span, af::span, k) = tempArf;
        af::array  arfTemp = convolve2(tempArf, kernel);
        af::array binnedFrame = arfTemp(af::seq(0, af::end, binning), af::seq(0, af::end, binning));
        af::array dx, dy;
        af::grad(dx, dy, binnedFrame);
        /// Sum of absolute of the X- and Y- gradient (look for a norm-L2 function in ArrayFire?)
        af::array gradientMagnitude = af::abs(dx) + af::abs(dy);
        gradientMagSeries(af::span, af::span, k) = gradientMagnitude;
        /// Sobel gradient magnitude
        ///gradientMagSeries(af::span, af::span, k) = af::sobel(binnedFrame);
    }
    qDebug("RProcessing::blockProcessingLocal:: afTimer1 elapsed seconds: %f s", af::timer::stop(afTimer1));

    /// Sobel gradient. Comment this if using the classical gradient above.
//    af::timer afTimer2 = af::timer::start();
//    af::array gradientMagSeries = af::sobel(arfSeriesBinned);
//    qDebug("RProcessing::blockProcessing:: afTimer2 elapsed seconds: %f s", af::timer::stop(afTimer2));


    /// Calculate number of blocks in the binned image. There must be enough block to cover the whole image, considering an overlap of 50%;
    int blkSize = 32;
    int blkSizeL = blkSize * 3;
//    int binnedBlkSize = blkSize/binning;
//    int nBBlksX = nBAxis1 / binnedBlkSize * 2; /// *2 from an overalp of 50% in the X-direction;
//    int nBBlksY = nBAxis2 / binnedBlkSize;
    //int nBPixels = nBAxis1 * nBAxis2;
//    int nBBlks = nBBlksX * nBBlksY;



//    af::timer afTimer3 = af::timer::start();

//    af::array canvas = constant(0, naxis2, naxis1);

//    /// Loop over all the block series
//    for (int k = 0; k < nBBlks; k++)
//    {
//        int blkOffsetX = (k % nBBlksY) * binnedBlkSize;
//        int blkOffsetY = (k / nBBlksX) * binnedBlkSize;

//        af::array binnedBlkSeries = moddims(gradientMagSeries(af::seq(blkOffsetY, blkOffsetY + binnedBlkSize -1), af::seq(blkOffsetX, blkOffsetX + binnedBlkSize -1), af::seq(0, nFrames-1)), af::dim4(binnedBlkSize*binnedBlkSize, nFrames));
//        af::array gradSum = flat(af::sum( binnedBlkSeries, 0 ));
//        /// Sort the sum of the gradient-norm
//        af::array sortedArray;
//        af::array sortIndices;
//        af::sort(sortedArray, sortIndices, gradSum, 0, false);

//        blkOffsetX *= binning;
//        blkOffsetY *= binning;

//        af::array blkSeries = arfSeries(af::seq(blkOffsetY, blkOffsetY + blkSize -1), af::seq(blkOffsetX, blkOffsetX + blkSize -1), af::seq(0, nFrames-1));
//        af::array sortedBlk = blkSeries(af::span, af::span, sortIndices);
//        af::array bestBlock = sortedBlk(af::span, af::span, 0);
//        canvas(af::seq(blkOffsetY, blkOffsetY + blkSize -1), af::seq(blkOffsetX, blkOffsetX + blkSize -1)) = bestBlock;
//    }

//    qDebug("RProcessing::blockProcessing:: afTimer3 elapsed seconds: %f s", af::timer::stop(afTimer3));


    af::timer afTimer4 = af::timer::start();

    af::array canvas = af::constant(0, naxis1, naxis2, f32);
    af::array weightsCanvas = af::constant(0, naxis1, naxis2, f32);

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
    canvas(af::seq(blkPos1.x(), blkPos1.x() + blkSize -1), af::seq(blkPos1.x(), blkPos1.x() + blkSize -1)) = bestBlk1;
    //weightsCanvas(af::seq(blkPos1.x(), blkPos1.x() + blkSize -1), af::seq(blkPos1.x(), blkPos1.x() + blkSize -1)) += 1;


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

            canvas(af::seq(newBlkPos1.x(), newBlkPos1.x() + blkSize -1), af::seq(newBlkPos1.y(), newBlkPos1.y() + blkSize -1)) = bestBlk1;
            //weightsCanvas(af::seq(newBlkPos1.x(), newBlkPos1.x() + blkSize -1), af::seq(newBlkPos1.y(), newBlkPos1.y() + blkSize -1)) += 1;

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

        canvas(af::seq(shiftedBlkPos2.x(), shiftedBlkPos2.x() + blkSize -1), af::seq(shiftedBlkPos2.y(), shiftedBlkPos2.y() + blkSize -1)) = bestBlk2;
        //weightsCanvas(af::seq(shiftedBlkPos2.x(), shiftedBlkPos2.x() + blkSize -1), af::seq(shiftedBlkPos2.y(), shiftedBlkPos2.y() + blkSize -1)) += 1;

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

    //canvas /= weightsCanvas;

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
//        af::array resultSlice = sortedBlk(af::span, af::span, k);
//        //af::array resultSlice = arfSeriesBinned(af::span, af::span, k);
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

//    af::array testCanvas = af::constant(0, 3*blkSize, 3*blkSize);
//    testCanvas(af::seq(0, blkSize-1), af::seq(Yoffset2, Yoffset2 + blkSize-1)) = bestBlk1;
//    testCanvas(af::seq(x2Shifted, x2Shifted + blkSize -1), af::seq(y2Shifted, y2Shifted + blkSize -1)) = bestBlk2;


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

    af::array binnedBlkSeries = moddims(arrayBinnedSeries(af::seq(blkPosB.x(), blkPosB.x() + binnedBlkSize -1), af::seq(blkPosB.y(), blkPosB.y() + binnedBlkSize -1), af::seq(0, nFrames-1)), af::dim4(binnedBlkSize*binnedBlkSize, nFrames));
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
    //af::array blkSeriesL = arraySeries(af::seq(blkPosL.x(), blkPosL.x() + blkSizeL -1), af::seq(blkPosL.y(), blkPosL.y() + blkSizeL -1), af::seq(0, nFrames-1));
    unsigned int *bestIndex = sortIndices.host<unsigned int>();
    searchBlk = arraySeries(af::seq(blkPosL.x(), blkPosL.x() + blkSizeL -1), af::seq(blkPosL.y(), blkPosL.y() + blkSizeL -1), bestIndex[0]);

    //af::array blkSeries = arraySeries(af::seq(blkPos.x(), blkPos.x() + blkSize -1), af::seq(blkPos.y(), blkPos.y() + blkSize -1), af::seq(0, nFrames-1));
    bestBlk = searchBlk(af::seq(blkSize, 2*blkSize -1), af::seq(blkSize, 2*blkSize -1));
}

void RProcessing::extractLuckySample(QList<RMat *> rMatImageList, QList<RMat *> & blockList, int x, int y, bool isCentered)
{

    if (isCentered)
    {
        x -= blkSize/2;
        y -= blkSize/2;
    }

    af::array arDim = af::constant(blkSize, 2, nBest-1);
    af::array coordRange = af::range(af::dim4(blkSize));
    af::array xRange = coordRange + x;
    af::array yRange = coordRange + y;

    af::array arfSeries;
    af::array qualityBinnedSeries;
    blockProcessingGradient(rMatImageList, arfSeries, qualityBinnedSeries);
    af::array stackedBlks;
    makeAlignedStackGradient(stackedBlks, arfSeries, qualityBinnedSeries, arDim, xRange, yRange, x, y);

    /// Pack the blocks into a list of RMat.
    for ( int k=0; k < nBest; k++)
    {
        af::array resultSlice = stackedBlks(af::span, af::span, k);
        float *hostArf = resultSlice.host<float>();
        cv::Mat matImage(blkSize, blkSize, CV_32F, hostArf);
        matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
        RMat *blockRMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
        blockRMat->setImageTitle("Block stack method 1 (gradient)");
        blockList << blockRMat;
        delete[] hostArf;
    }
}

void RProcessing::extractLuckySample2(QList<RMat *> rMatImageList, QList<RMat *> & blockList, int x, int y, bool isCentered)
{

    if (isCentered)
    {
        x -= blkSize/2;
        y -= blkSize/2;
    }

    af::array arDim = af::constant(blkSize, 2, nBest-1);
    af::array coordRange = af::range(af::dim4(blkSize));
    af::array xRange = coordRange + x;
    af::array yRange = coordRange + y;

    af::array arfSeries;
    af::array qualitySeries;
    blockProcessingLaplace(rMatImageList, arfSeries, qualitySeries);
    af::array globalRefImage = af::median(arfSeries, 2);

    af::array stackedBlks;
    makeAlignedStackLaplace(stackedBlks, arfSeries, qualitySeries, arDim, xRange, yRange, x, y);
    /// Align the stackedBlks on the global reference image
    /// Use the 1st slice of the stackedBlks
    af::array blkSlice = stackedBlks.slice(0);
    af::array refBlk = globalRefImage(xRange, yRange);
    af::array arDim2 = af::constant(blkSize, 2);
    af::array shift = af::constant(0, 2, f32);
    phaseCorrelate(refBlk, blkSlice, arDim2, shift);
    std::printf("Shift of stack block: \n");
    af_print(shift);
    /// The shifted position here are in the reference frame of the stacked block,
    /// so we need to add the shift instead of subtracting it. We are indeed truly shifting
    /// the position of the block.

    af::array xr = xRange + tile(shift(0), blkSize);
    af::array yr = yRange + tile(shift(1), blkSize);

    /// Pack the blocks into a list of RMat.
    for ( int k=0; k < nBest; k++)
    {
        af::array resultSlice = stackedBlks(af::span, af::span, k);
        float *hostArf = resultSlice.host<float>();
        cv::Mat matImage(blkSize, blkSize, CV_32F, hostArf);
        matImage.convertTo(matImage, rMatImageList.at(0)->matImage.type());
        RMat *blockRMat = new RMat(matImage, rMatImageList.at(0)->isBayer(), rMatImageList.at(0)->getInstrument());
        blockRMat->setImageTitle("Block stack method 2 (Laplace)");
        blockList << blockRMat;
        delete[] hostArf;
    }
}

void RProcessing::blockProcessingGradient(QList<RMat *> rMatImageList, af::array &arfSeries, af::array &qualityBinnedSeries)
{

    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();

    /// Initialize the GPU array series
    int nBAxis1 = naxis1/binning;
    int nBAxis2 = naxis2/binning;
    arfSeries = af::constant(0, naxis1, naxis2, nFrames);
    qualityBinnedSeries = af::constant(0, nBAxis1, nBAxis2, nFrames);

    af::array kernel = af::constant(0, 7, 7);
    kernel(af::seq(3,6), af::seq(3,6)) = 1;

    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempArf(naxis1, naxis2, (float*) tempMat.data);
        arfSeries(af::span, af::span, k) = tempArf;
        /// Rebin the data
        af::array arfTemp = convolve2(tempArf, kernel);
        af::array binnedFrame = arfTemp(af::seq(0, af::end, binning), af::seq(0, af::end, binning));
        af::array dx, dy;
        af::grad(dx, dy, binnedFrame);
        /// Sum of absolute of the X- and Y- gradient (look for a norm-L2 function in ArrayFire?)
        af::array gradientMagnitude = af::abs(dx) + af::abs(dy);
        qualityBinnedSeries(af::span, af::span, k) = gradientMagnitude;
    }
}



void RProcessing::blockProcessingSobel(QList<RMat *> rMatImageList, af::array &arfSeries, af::array &qualitySeries)
{

    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();

    /// Initialize the GPU array series
    int nBAxis1 = naxis1/binning;
    int nBAxis2 = naxis2/binning;
    arfSeries = af::constant(0, naxis1, naxis2, nFrames);
    qualitySeries = af::constant(0, nBAxis1, nBAxis2, nFrames);

    af::array kernel = af::constant(0, 7, 7);
    kernel(af::seq(3,6), af::seq(3,6)) = 1;

    cv::Mat tempMat(naxis2, naxis1, CV_32F);
    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempArf(naxis1, naxis2, (float*) tempMat.data);
        arfSeries(af::span, af::span, k) = tempArf;
        /// Rebin the data
        af::array arfTemp = convolve2(tempArf, kernel);
        af::array binnedFrame = arfTemp(af::seq(0, af::end, binning), af::seq(0, af::end, binning));
        af::array dx, dy;
        //af::grad(dx, dy, binnedFrame);
        af::sobel(dx, dy, binnedFrame);
        /// Sum of absolute of the X- and Y- gradient (look for a norm-L2 function in ArrayFire?)
        af::array gradientMagnitude = af::abs(dx) + af::abs(dy);
        qualitySeries(af::span, af::span, k) = gradientMagnitude;
    }
}


void RProcessing::blockProcessingLaplace(QList<RMat *> rMatImageList, af::array &arfSeries, af::array &qualitySeries)
{
    /// Prepare the Laplace metric for lucky imaging
    /// function makeAlignedStackLaplace will then calculate the actual quality values
    int naxis2 = rMatImageList.at(0)->matImage.rows;
    int naxis1 = rMatImageList.at(0)->matImage.cols;
    int nFrames = rMatImageList.size();

    // Initialize the GPU array series
    arfSeries = af::constant(0, naxis1, naxis2, nFrames);
    qualitySeries = af::constant(0, naxis1/binning, naxis2/binning, nFrames);

//    float ker[] = {0, 1, 0,
//                   1,-4, 1,
//                   0, 1, 0};

    float ker[] = {1, 4, 1,
                   4,-20,4,
                   1, 4, 1};

    af::array kernel(3, 3, ker);

    cv::Mat tempMat(naxis2, naxis1, CV_32F);

    for ( int k=0; k < nFrames; k++)
    {
        rMatImageList.at(k)->matImage.convertTo(tempMat, CV_32F);
        af::array tempAr(naxis1, naxis2, (float*) tempMat.data);
        arfSeries(af::span, af::span, k) = tempAr;
        // Rebin the image. The image is in fact not rebinned if binning = 1.
        af::array binnedAr;
        rebin(tempAr, binnedAr, binning);
        // Convolve the rebinned with Laplacian kernel
        af::array arfTemp = convolve2(binnedAr, kernel);
        // Store it in the 3D array
        qualitySeries(af::span, af::span, k) = arfTemp;
    }
}

void RProcessing::blockProcessingGlobalGradients(QList<RMat *> rMatImageList)
{
    af::setBackend(AF_BACKEND_OPENCL);
    //af::setBackend(AF_BACKEND_CPU);

    af::array arfSeries;
    af::array qualityBinnedSeries;

    if (this->qualityMetric == QString("gradient"))
    {
        blockProcessingGradient(rMatImageList, arfSeries, qualityBinnedSeries);
    }
    else if (this->qualityMetric == QString("Sobel"))
    {
        blockProcessingSobel(rMatImageList, arfSeries, qualityBinnedSeries);
    }

    af::array globalRefImage = af::median(arfSeries, 2);

    int naxis1 = arfSeries.dims(0);
    int naxis2 = arfSeries.dims(1);
    af::array canvas3D = af::constant(0, naxis1, naxis2, nBest, f32);
    af::array weightCanvas = af::constant(0, naxis1, naxis2, nBest, f32);
    af::array canvas = af::constant(0, naxis1, naxis2, f32);
    // Need an array tiled with the dimensions of the block, for the shift matrix.

    af::array arDim = af::constant(blkSize, 2, nBest-1);
    af::array arDim2 = af::constant(blkSize, 2);
    af::array coordRange = af::range(af::dim4(blkSize));

    int offsetX = blkSize;
    int offsetY = blkSize;


    int naxis1End = naxis1 - 1 - 3*blkSize;
    int naxis2End = naxis2 - 1 - 3*blkSize;


    int x = 0;
    int y = 0;
    int k = 0;

    af::sync();
    af::timer afTimer = af::timer::start();

    while ( (y < naxis2End) || (x < naxis1End) )
    {
        x = offsetX + (k*blkSize/2 % naxis1End);
        y = offsetY + ((k*blkSize/2) / naxis1End)*blkSize/2;
        //qDebug("[x, y] = %d, %d", x, y);
        af::array xRange = coordRange + x;
        af::array yRange = coordRange + y;

        af::array stackedBlks;
        makeAlignedStackGradient(stackedBlks, arfSeries, qualityBinnedSeries, arDim, xRange, yRange, x, y);
        /// We should not use the median on the stackedBlks just yet...

        /// Align the stackedBlks on the global reference image
        /// Use the 1st slice of the stackedBlks
        af::array blkSlice = stackedBlks.slice(0);
        af::array refBlk = globalRefImage(xRange, yRange);
        af::array shift = af::constant(0, 2, f32);
        phaseCorrelate(refBlk, blkSlice, arDim2, shift);

        /// The shifted position here are in the reference frame of the stacked block,
        /// so we need to add the shift instead of subtracting it. We are indeed truly shifting
        /// the position of the block.

        af::array xr = xRange + tile(shift(0), blkSize);
        af::array yr = yRange + tile(shift(1), blkSize);

        canvas3D(xr, yr, af::span) += stackedBlks;
        weightCanvas(xr, yr, af::span) += 1;

        k++;

        // 20 ms
    }

    af::array mask = weightCanvas == 0;
    weightCanvas(mask) = 1;
    canvas3D /= weightCanvas;

    canvas = median(canvas3D, 2);

    int axisStart  = 1.5 * blkSize;
    canvas(af::seq(0,axisStart), af::span)   = globalRefImage(af::seq(0,axisStart), af::span);
    canvas(af::span, af::seq(0,axisStart))   = globalRefImage(af::span, af::seq(0,axisStart));
    canvas(af::seq(naxis2End, af::end), af::span)  = globalRefImage(af::seq(naxis2End, af::end), af::span);
    canvas(af::span, af::seq(naxis1End, af::end))  = globalRefImage(af::span, af::seq(naxis1End, af::end));


    af::sync();
    double totalTime = af::timer::stop(afTimer);
    int k2 = k-1;
    qDebug("ProcessingGlobalStack2:: total time for %d loops = %f s", k2, totalTime);
    qDebug("ProcessingGlobalStack2:: average time per loop = %f ms", totalTime / k2 * 1000.0);

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    /// Slice in the 1st original image
    //resultList << rMatImageList.at(0);
    populateResultListWithAr(rMatImageList, canvas, QString("with stacked blocks"));

}

void RProcessing::blockProcessingGlobalLaplace(QList<RMat *> rMatImageList)
{
    af::setBackend(AF_BACKEND_OPENCL);
    //af::setBackend(AF_BACKEND_CPU);

    af::array arfSeries;
    af::array qualityBinnedSeries;
    blockProcessingLaplace(rMatImageList, arfSeries, qualityBinnedSeries);
    af_print(qualityBinnedSeries(af::seq(0, 10), af::seq(0, 10), af::seq(0, 0)))
    af::array globalRefImage = af::median(arfSeries, 2);

    int naxis1 = arfSeries.dims(0);
    int naxis2 = arfSeries.dims(1);
    af::array canvas3D = af::constant(0, naxis1, naxis2, nBest, f32);
    af::array weightCanvas = af::constant(0, naxis1, naxis2, nBest, f32);
    af::array canvas = af::constant(0, naxis1, naxis2, f32);
    // Need an array tiled with the dimensions of the block, for the shift matrix.

    af::array arDim = af::constant(blkSize, 2, nBest-1);
    af::array arDim2 = af::constant(blkSize, 2);

    af::array coordRange = af::range(af::dim4(blkSize));

    int offsetX = blkSize;
    int offsetY = blkSize;

    int naxis1End = naxis1 - 1 - 3*blkSize;
    int naxis2End = naxis2 - 1 - 3*blkSize;

    int x = 0;
    int y = 0;
    int k = 0;

    af::sync();
    af::timer afTimer = af::timer::start();

    while ( (y < naxis2End) || (x < naxis1End) )
    {
        x = offsetX + (k*blkSize/2 % naxis1End);
        y = offsetY + ((k*blkSize/2) / naxis1End)*blkSize/2;
        //qDebug("[x, y] = %d, %d", x, y);
        af::array xRange = coordRange + x;
        af::array yRange = coordRange + y;

        af::array stackedBlks;
        makeAlignedStackLaplace(stackedBlks, arfSeries, qualityBinnedSeries, arDim, xRange, yRange, x, y);

        /// Align the stackedBlks on the global reference image
        /// Use the 1st slice of the stackedBlks
        af::array blkSlice = stackedBlks.slice(0);
        af::array refBlk = globalRefImage(xRange, yRange);
        af::array shift = af::constant(0, 2, f32);
        phaseCorrelate(refBlk, blkSlice, arDim2, shift);

        /// The shifted position here are in the reference frame of the stacked block,
        /// so we need to add the shift instead of subtracting it. We are indeed truly shifting
        /// the position of the block.

        af::array xr = xRange + tile(shift(0), blkSize);
        af::array yr = yRange + tile(shift(1), blkSize);

        canvas3D(xr, yr, af::span) += stackedBlks;
        weightCanvas(xr, yr, af::span) += 1;

        k++;
    }

    af::array mask = weightCanvas == 0;
    weightCanvas(mask) = 1;
    canvas3D /= weightCanvas;

    canvas = median(canvas3D, 2);

    int axisStart  = 1.5 * blkSize;
    canvas(af::seq(0,axisStart), af::span)   = globalRefImage(af::seq(0,axisStart), af::span);
    canvas(af::span, af::seq(0,axisStart))   = globalRefImage(af::span, af::seq(0,axisStart));
    canvas(af::seq(naxis2End, af::end), af::span)  = globalRefImage(af::seq(naxis2End, af::end), af::span);
    canvas(af::span, af::seq(naxis1End, af::end))  = globalRefImage(af::span, af::seq(naxis1End, af::end));

    af::sync();
    double totalTime = af::timer::stop(afTimer);
    int k2 = k-1;
    qDebug("ProcessingGlobalStack2:: total time for %d loops = %f s", k2, totalTime);
    qDebug("ProcessingGlobalStack2:: average time per loop = %f ms", totalTime / k2 * 1000.0);

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    /// Slice in the 1st original image
    //resultList << rMatImageList.at(0);
    populateResultListWithAr(rMatImageList, canvas, QString("with stacked blocks"));



}

void RProcessing::extractBestBlock(af::array &bestBlk, af::array &arfSeries, af::array &arrayBinnedSeries,
                                     const int &blkSize, const int &binnedBlkSize, const int &x, const int &y,
                                    const int &nFrames, const int &binning)
{
    int xB = x/binning;
    int yB = y/binning;

    af::array binnedBlkSeries = moddims(arrayBinnedSeries(af::seq(xB, xB + binnedBlkSize -1), af::seq(yB, yB + binnedBlkSize -1), af::seq(0, nFrames-1)), af::dim4(binnedBlkSize*binnedBlkSize, nFrames));
    af::array gradSum = flat(af::sum( binnedBlkSeries, 0 ));
    //af_print(gradSum);

    /// Sort the sum of the gradient-norm
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, gradSum, 0, false);

    unsigned int *bestIndices = sortIndices.host<unsigned int>();
    bestBlk = arfSeries(af::seq(x, x + blkSize -1), af::seq(y, y + blkSize -1), bestIndices[0]);
}

void RProcessing::makeAlignedStack(af::array &stackedBlks, const af::array &arfSeries, const af::array &qualityBinnedSeries, const int nBest, const int &blkSize, const int &binnedBlkSize, int &x, int &y, const int bufferSpace, const int &binning)
{
    int nFrames = arfSeries.dims(2);
    int xB = x/binning;
    int yB = y/binning;

    af::array binnedBlkSeries = moddims(qualityBinnedSeries(af::seq(xB, xB + binnedBlkSize -1), af::seq(yB, yB + binnedBlkSize -1), af::seq(0, nFrames-1)), af::dim4(binnedBlkSize*binnedBlkSize, nFrames));
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
    af::array refBlk = arfSeries(af::seq(xL, xL + blkSizeL -1), af::seq(yL, yL + blkSizeL -1), bestIndices[0]);
    /// The refBlk is the reference block. But it will be used as the template to search in the bigger block
    /// (bestSearchBlks below) of size blkSizeL that will be cropped afterwards.
    af::array selectedInds = sortIndices(af::seq(0, nBest-1));


    af::array bestBlks = arfSeries(af::seq(x, x + blkSize -1), af::seq(y, y + blkSize -1), selectedInds);
    /// Now let's start aligning the template onto those bigger bloks.
    af::array SADs;
    matchTemplate2(SADs, refBlk, bestBlks);
    af::array dxArr, dyArr;
    fetchTMatch2Shifts(SADs, dxArr, dyArr, bestBlks.dims());
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

    stackedBlks = af::constant(0, blkSize, blkSize, nBest);

//    af::sync();
//    af::timer timerFor = af::timer::start();

    for (int i = 0; i < nBest; i++)
    {
//        af::array xCropi = xCropRange(af::span, i);
//        af::array yCropi = yCropRange(af::span, i);
//        stackedBlks(af::span, af::span, i) = arfSeries(xCropi, yCropi, bestIndices[i]);
        int x1 = x - dxHost[i] ;
        int x2 = x1 + blkSize - 1;
        int y1 = y - dyHost[i];
        int y2 = y1 + blkSize - 1;
        stackedBlks(af::span, af::span, i) = arfSeries(af::seq(x1, x2), af::seq(y1, y2), bestIndices[i]);
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
    af::array onesAr = af::constant(1, img.dims());
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

void RProcessing::phaseCorrelate(af::array &refBlk, af::array &shiftedArray, const af::array &arDim, af::array &shifts)
{
    /// Calculate the shifts necessary to to align the shiftedArray onto array.
    /// See equivalent code at https://github.com/scikit-image/scikit-image/blob/master/skimage/feature/register_translation.py#L109

    af::array arrayProduct = fft2(refBlk) * conjg(fft2(shiftedArray));
    af::array cc = af::abs(ifft2(arrayProduct));
    findMaxLoc(cc, shifts);
    af::array mask = shifts > arDim / 2;
    shifts(mask) -= arDim(mask);
}

void RProcessing::phaseCorrelate2(af::array &stackedBlks, af::array &shifts, const af::array &arDim)
{
    /// Same as phaseCorrelate but use a GPU-batched Fourier-transform.
    /// instead of doing it in a cpu-loop.

//    af::array (*rMatMulPtr)(const af::array&, const af::array&);
//    rMatMulPtr = RProcessing::rMatMul;

    af::array refBlk = stackedBlks.slice(0);
    af::array blks = stackedBlks.slices(1, stackedBlks.dims(2)-1);
    af::array refBlkF = fft2(refBlk);
    af::array blksF = conjg(fft2(blks));
    af::array arrayProduct = batchFunc(blksF, refBlkF, &RProcessing::rMatMul);
    af::array cc = af::abs(ifft2(arrayProduct));

    findMaxLoc2(cc, shifts);
    /// Fix the shifts, considering square blocks. shifts(0, af::span) are x-shifts, shifts(1, af::span) are y-shifts
    af::array mask = shifts > arDim / 2;
    shifts(mask) -= arDim(mask);
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
    af::array val, idx0;
    af::min(val, idx0, af::moddims(a, a.dims(0)*a.dims(1), a.dims(2)));
    dx = idx0 % a.dims(0);
    dy = idx0 / a.dims(0);
}

void RProcessing::findMaxLoc(af::array &a, af::array &locxy)
{
    /// Get the location of the maximum in each image or 2D array in a series.
    /// Works in 2D.
    af::array val, idx;
    af::max(val, idx, af::moddims(a, a.dims(0)*a.dims(1), a.dims(2)));
    locxy(0) = idx % a.dims(0);
    locxy(1) = idx / a.dims(0);
}

void RProcessing::findMaxLoc2(af::array &a, af::array &locxy)
{
    /// Get the location of the minimum in each image or 2D array in a series.
    /// Works in 3D.
    af::array val, idx0;
    af::max(val, idx0, af::moddims(a, a.dims(0)*a.dims(1), a.dims(2)));
    locxy(0, af::span) = idx0 % a.dims(0);
    locxy(1, af::span) = idx0 / a.dims(0);
}



void RProcessing::fetchTMatch2Shifts(af::array & a, int &dx, int &dy, af::dim4 dims)
{
    /// Used with matchTemplate2 to the get correct the location of the minimum
    findMinLoc(a, dx, dy);
    dx -= (dims[0]/2 -1);
    dy -= (dims[0]/2 -1);
}

void RProcessing::fetchTMatch2Shifts(af::array & a, af::array &dx, af::array &dy, af::dim4 dims)
{
    findMinLoc(a, dx, dy);
    dx -= (dims[0]/2 -1);
    dy -= (dims[0]/2 -1);
}

void RProcessing::fetchTMatchShiftsGfor(af::array & ar, af::array &dxAr, af::array &dyAr)
{
    float minVal = 0;
    unsigned minLoc = 0;

    gfor(af::seq i, ar.dims(2))
    {
        af::min<float>(&minVal, &minLoc, ar(af::span, af::span, i));
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

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    std::cout << "Loading master Bias..." << std::endl;
    loadMasterBias();


    std::cout << "Loading master Dark..." << std::endl;
    loadMasterDark();

    std::cout << "Loading master Flat..." << std::endl;
    loadMasterFlat();

    /// Flat fielding needs also to have at least the bias removed.
    if (masterBias !=NULL)
    {
        std::cout << "Subtracking Bias to Flat..." << std::endl;
        cv::subtract(masterFlat->matImage, masterBias->matImage, masterFlat->matImage);
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
    std::cout << "RProcessing::calibrate()" << std::endl;
    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    std::cout << "Calibrating Lights..." << std::endl;
    for(int i = 0; i < treeWidget->getLightUrls().size(); i++)
    {
        ImageManager lightManager(treeWidget->getLightUrls().at(i));
        //cv::Mat matResult;
        cv::Mat lightMat = lightManager.getRMatImage()->matImage;

        if (lightMat.type() != CV_32F)
        {
            lightMat.convertTo(lightMat, CV_32F);
            std::cout << "Converted light to CV_32F" << std::endl;
        }

        if (masterDark!=NULL )
        {
            std::cout << "Subtracting Dark..." << std::endl;
            cv::subtract(lightMat, masterDark->matImage, lightMat);
        }
        else if (masterBias != NULL)
        {
            cv::subtract(lightMat, masterBias->matImage, lightMat);
        }

        if (masterFlat != NULL)
        {
            std::cout << "Flat fielding..." << std::endl;
            /// A flat field completely modify the range of the image.
            /// We could restore that range and continue treating the image with an integral type
            /// or we could treat it as float as long as we stay consistent.
            /// In case we treat as float, demosaicing with openCV needs image to range from 0 and 1
            cv::Scalar oldMean = cv::mean(lightMat);
            std::cout << "Dividing by flat-field..." << std::endl;
            cv::divide(lightMat, masterFlat->matImage, lightMat);
            cv::Scalar newMean = cv::mean(lightMat);
            cv::Scalar scale = oldMean/newMean;
            lightMat = lightMat * (float)scale.val[0];


        }

        cv::threshold(lightMat, lightMat, 0, 0, cv::THRESH_TOZERO);

        /// Copy the result into the list.
        /// This may be an overkill. We could output directly into the elements of the list
        resultList << new RMat(lightMat, lightManager.getRMatImage()->isBayer(), lightManager.getRMatImage()->getInstrument(), lightManager.getRMatImage()->getXPOSURE(), lightManager.getRMatImage()->getTEMP());
        resultList.at(i)->flipUD = lightManager.getRMatImage()->flipUD;
        resultList.at(i)->calcMinMax();
        resultList.at(i)->calcStats();
        resultList.at(i)->setFileInfo(lightManager.getRMatImage()->getFileInfo());
        showMinMax(resultList.at(i)->matImageGray);
    }

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

bool RProcessing::prepRegistration()
{
    // If the treeWidget is not empty, it can also be from a currentROpenGLWidget
    // It makes sense to use the displayFirstCheckBox button to determine whether
    // we load it from the treeWidget or from that currentROpenGLWidget.
    if (!treeWidget->rMatLightList.isEmpty() && !useUrlsFromTreeWidget)
    {
        rMatLightList = treeWidget->rMatLightList;
    }
    else if(!treeWidget->getLightUrls().empty())
    {
        emit tempMessageSignal(QString("Batch processing %1 images").arg(rMatLightList.size()), 10000);
        loadRMatLightList(treeWidget->getLightUrls());
    }
    else
    {
        emit tempMessageSignal(QString("No lights to register"), 10000);
        return false;
    }
    if (cvRectROIList.isEmpty())
    {
        if (!cvRectROI.empty())
        {
            cvRectROIList.append(cvRectROI);
        }
        else
        {
            emit tempMessageSignal(QString("ROI not defined"), 10000);
            return false;
        }
    }
    // Clear the resultList if not empty
    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    if (rMatLightList.at(0)->isBayer())
    {
        resultList << new RMat(rMatLightList.at(0)->matImageRGB, false, rMatLightList.at(0)->getInstrument(), rMatLightList.at(0)->getXPOSURE(), rMatLightList.at(0)->getTEMP());
        resultList.at(0)->setFileInfo(rMatLightList.at(0)->getFileInfo());
    }
    else
    {
        resultList << new RMat(rMatLightList.at(0)->matImage, false, rMatLightList.at(0)->getInstrument(), rMatLightList.at(0)->getXPOSURE(), rMatLightList.at(0)->getTEMP());
        resultList.at(0)->setFileInfo(rMatLightList.at(0)->getFileInfo());
    }


    return true;

}


void RProcessing::registerSeries()
{
    /// Here we register rMatLightList. It is assigned in two ways:
    /// 1) By Drag and Drop in the QMdiArea
    /// 2) After a calibration like in calibrate()
    /// We use resultList as the (temporary?) output list.

    bool status = prepRegistration();
    if (!status) {return;}

    int nFrames = rMatLightList.size();

    // Define the motion model
    const int warp_mode_1 = cv::MOTION_TRANSLATION;
    const int warp_mode_2 = cv::MOTION_TRANSLATION;

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
    cv::Mat refMat;
    // Normalized version
    cv::Mat refMatN;

    rMatLightList.at(0)->matImageGray.convertTo(refMat, CV_32F);
    // Normalize to a multiple of the exposure time * median?
    //float lowThresh = rMatLightList.at(0)->getIntensityLow();
    //float highThresh = rMatLightList.at(0)->getIntensityHigh();
    float normFactor = 1.0f / rMatLightList.at(0)->getXPOSURE();

    if(applyMask & (maskCircleRadius !=0))
    {
        refMat = circleMaskMat(refMat, maskCircleX, maskCircleY, maskCircleRadius);
        //emit resultSignal(refMat, false, instruments::DSLR);
    }

    // Normalize to the normFactor (e.g: the mean, high threshold from rMat->calcStats(), ...)
    refMatN = refMat * normFactor;
    //cv::threshold(refMatN, registeredMatN, 0.9, 0.9, cv::THRESH_TRUNC);


    // If ROI is used
    if (useROI)
    {
        refMatN = refMatN(cvRectROI);
        emit resultSignal(refMatN, false, rMatLightList.at(0)->getInstrument());
    }

    // Get a resampled version. 1/4 on each axis.
    cv::Mat refMatR;
    cv::resize(refMatN, refMatR, cv::Size(), 0.25, 0.25, CV_INTER_AREA);


    for (int i = 1 ; i < nFrames; ++i)
    {
        qDebug("Registering image #%i/%i", i, nFrames);

        cv::Mat registeredMat;
        // Normalized mat image
        cv::Mat registeredMatN;
        // Rebinned mat Image
        cv::Mat registeredMatR;

        // Registration requires floating points
        rMatLightList.at(i)->matImageGray.convertTo(registeredMat, CV_32F);

        if(applyMask & (maskCircleRadius !=0))
        {
            registeredMat = circleMaskMat(registeredMat, maskCircleX, maskCircleY, maskCircleRadius);
            //emit resultSignal(refMat, false, instruments::DSLR);
        }

        registeredMatN = registeredMat * normFactor;
        // If ROI is used
        if (useROI)
        {
            registeredMatN = registeredMatN(cvRectROI);
        }

        //cv::threshold(registeredMatN, registeredMatN, 0.9, 0.9, cv::THRESH_TRUNC);



        // Get a resampled version. 1/4 on each axis;
        cv::resize(registeredMatN, registeredMatR, cv::Size(), 0.25, 0.25, CV_INTER_AREA);

        cv::Mat warp_matrix_1 = cv::Mat::eye(2, 3, CV_32F);
        double eccEps = 0;
        // 1st pass of the ECC algorithm on the decimated image. The results are stored in warp_matrix.
        eccEps = cv::findTransformECC(
                    refMatR,
                    registeredMatR,
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
                    refMatN,
                    registeredMatN,
                    warp_matrix_1,
                    warp_mode_2,
                    criteria2
                    );

        qDebug() << "eccEps 2 =" << eccEps;
        std::cout << "result warp_matrix 2 =" << std::endl << warp_matrix_1 << std::endl << std::endl;

        // To do the alignment, use warpAffine. Needs to be applied on 3 channels separately and reassemble the channels?
        // That should be ok as the shifts will be applied rigidly on the 3 channels.

        if (rMatLightList.at(0)->isBayer())
        {
            // RGB array for splitting channels
            cv::Mat tempMatRGB[3];
            // Split the channels
            cv::split(rMatLightList.at(i)->matImageRGB, tempMatRGB);
            cv::warpAffine(tempMatRGB[0], tempMatRGB[0], warp_matrix_1, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
            cv::warpAffine(tempMatRGB[1], tempMatRGB[1], warp_matrix_1, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
            cv::warpAffine(tempMatRGB[2], tempMatRGB[2], warp_matrix_1, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);

            tempMatRGB[0].convertTo(tempMatRGB[0], CV_16U);
            tempMatRGB[1].convertTo(tempMatRGB[1], CV_16U);
            tempMatRGB[2].convertTo(tempMatRGB[2], CV_16U);
            vector<cv::Mat> channels;
            // Watch the order. 1st array will be last channel.
            channels.push_back(tempMatRGB[0]);
            channels.push_back(tempMatRGB[1]);
            channels.push_back(tempMatRGB[2]);
            cv::Mat registeredMatRGB;
            cv::merge(channels, registeredMatRGB);
            // registeredMat is necessarily non-bayer.
            resultList << new RMat(registeredMatRGB, false, rMatLightList.at(i)->getInstrument());
        }
        else
        {
            cv::warpAffine(registeredMat, registeredMat, warp_matrix_1, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
            // registeredMat is necessarily non-bayer.
            resultList << new RMat(registeredMat, false, rMatLightList.at(i)->getInstrument());
            resultList.at(i)->setBscale(normFactor);

        }

    }

}

void RProcessing::registerSeriesXCorrPropagate()
{
    /// Register by pairs of nearest frames (in time) and propagate up to the first so that the series
    /// is coaligned with the first image.

    // Check that registration is properly setup, including normalization.
    bool status = prepRegistration();

    if (!status)
    {
        return;
    }

    cv::Mat refMat;
    cv::Mat currentMatImage;
    cv::Mat warpMatrixTotal = cv::Mat::eye( 2, 3, CV_32FC1 );
    for (int i=0; i < rMatLightList.size()-1; i++)
    {
        rMatLightList.at(i)->matImageGray.convertTo(refMat, CV_32F);
        rMatLightList.at(i+1)->matImageGray.convertTo(currentMatImage, CV_32F);

        cv::Mat refMatN = refMat / rMatLightList.at(i)->getXPOSURE();
        cv::Mat currentMatImageN = currentMatImage / rMatLightList.at(i+1)->getXPOSURE();

        //shift = shift + calculateSADShift(refMatN, currentMatImageN, cvRectROIList, 50);
        cv::Mat warpMat = calculateXCorrShift(refMatN, currentMatImageN, cvRectROIList);
        warpMatrixTotal.at<float>(0, 2) += warpMat.at<float>(0, 2);
        warpMatrixTotal.at<float>(1, 2) += warpMat.at<float>(1, 2);

        std::cout << "RProcessing::registerSeriesXCorrPropagate() frame # " << i << std::endl;
        std::cout << "RProcessing::registerSeriesXCorrPropagate() file: " << rMatLightList.at(i)->getFileInfo().baseName().toStdString() << std::endl;
        std::cout << "RProcessing::registerSeriesXCorrPropagate() ShiftX = " << warpMat.at<float>(0, 2) << std::endl;
        std::cout << "RProcessing::registerSeriesXCorrPropagate() ShiftY = " << warpMat.at<float>(1, 2) << std::endl;

        cv::Mat shiftedMat = shiftImage(rMatLightList.at(i+1), warpMatrixTotal);
        resultList << new RMat(shiftedMat, false, rMatLightList.at(i+1)->getInstrument(), rMatLightList.at(i+1)->getXPOSURE(), rMatLightList.at(i+1)->getTEMP());
        resultList.at(i+1)->setFileInfo(rMatLightList.at(i+1)->getFileInfo());
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
    //const int warp_mode_2 = cv::MOTION_EUCLIDEAN;

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
    /// Here we register rMatLightList. It is assigned in two ways:
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
    // Clear the resultList if not empty
    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    // Get the 1st image of the rMatLightList as the reference image (and put as 1st element of resultList)
    cv::Mat refMat;
    // Normalized version
    cv::Mat refMatN;

    if (rMatLightList.at(0)->isBayer())
    {
        resultList << new RMat(rMatLightList.at(0)->matImageRGB, false, rMatLightList.at(0)->getInstrument());
    }
    else
    {
        resultList << new RMat(rMatLightList.at(0)->matImage, false, rMatLightList.at(0)->getInstrument());
    }


    rMatLightList.at(0)->matImageGray.convertTo(refMat, CV_32F);
    // Normalize to a multiple of the exposure time * median?
    float normFactor = 1.0f / rMatLightList.at(0)->getXPOSURE();
    // Normalize to the normFactor (e.g: the mean, high threshold from rMat->calcStats(), ...)
    refMatN = refMat * normFactor;


    for (int i=1; i < rMatLightList.size(); i++)
    {
        cv::Mat currentMatImage;
        rMatLightList.at(i)->matImageGray.convertTo(currentMatImage, CV_32F);

        cv::Mat currentMatImageN = currentMatImage * normFactor;

        cv::Point2d shift = cv::phaseCorrelate(refMatN, currentMatImageN);
        std::cout << "Shifts = " << shift << std::endl;

        cv::Mat registeredMat;
        cv::Mat warpMat = cv::Mat::eye(2, 3, CV_32F);
        warpMat.at<float>(0, 2) = shift.x;
        warpMat.at<float>(1, 2) = shift.y;

        if (rMatLightList.at(0)->isBayer())
        {
            // RGB array for splitting channels
            cv::Mat tempMatRGB[3];
            // Split the channels
            cv::split(rMatLightList.at(i)->matImageRGB, tempMatRGB);
            cv::warpAffine(tempMatRGB[0], tempMatRGB[0], warpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
            cv::warpAffine(tempMatRGB[1], tempMatRGB[1], warpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
            cv::warpAffine(tempMatRGB[2], tempMatRGB[2], warpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);

            tempMatRGB[0].convertTo(tempMatRGB[0], CV_16U);
            tempMatRGB[1].convertTo(tempMatRGB[1], CV_16U);
            tempMatRGB[2].convertTo(tempMatRGB[2], CV_16U);
            vector<cv::Mat> channels;
            // Watch the order. 1st array will be last channel.
            channels.push_back(tempMatRGB[0]);
            channels.push_back(tempMatRGB[1]);
            channels.push_back(tempMatRGB[2]);
            cv::Mat registeredMatRGB;
            cv::merge(channels, registeredMatRGB);
            // registeredMat is necessarily non-bayer.
            resultList << new RMat(registeredMatRGB, false, rMatLightList.at(i)->getInstrument());
        }
        else
        {
            cv::warpAffine(rMatLightList.at(i)->matImage, registeredMat, warpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
            // registeredMat is necessarily non-bayer.
            resultList << new RMat(registeredMat, false, rMatLightList.at(i)->getInstrument());
            resultList.at(i)->setBscale(normFactor);

        }
    }

}

void RProcessing::registerSeriesCustom()
{
    if (cvRectROI.empty())
    {
        emit tempMessageSignal(QString("ROI not defined"));
        return;
    }

    if(applyMask & (maskCircleRadius ==0))
    {
        emit tempMessageSignal(QString("circular mask not defined"));
        return;
    }

    if (!treeWidget->rMatLightList.isEmpty())
    {
        rMatLightList = treeWidget->rMatLightList;
    }
    else
    {
        emit tempMessageSignal(QString("No lights to register"));
        return;
    }
    // Clear the resultList if not empty
    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    // Get the 1st image of the rMatLightList as the reference image (and put as 1st element of resultList)
    cv::Mat refMat;
    // Normalized version
    cv::Mat refMatN;

    if (rMatLightList.at(0)->isBayer())
    {
        resultList << new RMat(rMatLightList.at(0)->matImageRGB, false, rMatLightList.at(0)->getInstrument());
    }
    else
    {
        resultList << new RMat(rMatLightList.at(0)->matImage, false, rMatLightList.at(0)->getInstrument());
    }

    rMatLightList.at(0)->matImageGray.convertTo(refMat, CV_32F);

    float normFactor = 1.0f / rMatLightList.at(0)->getXPOSURE();
    refMatN = refMat * normFactor;

    std::cout << "cvRectROI = " << cvRectROI << std::endl;

    // Originally I chose to use the first image in the timeline as a reference
    // During the total eclipse, statistical properties change rapidly, and normalization by exposure time is not enough
    // (it's also the case with clouds, fire haze passing quickly etc...)
    // Thus I should opt for a propagating approach, and make a by-pair coalignment, and propagate the shift with respect to the first image
    // This is the purpose of the function registerSeriesCustomPropagate()
    for (int i=1; i < rMatLightList.size(); i++)
    {
        cv::Mat currentMatImage;
        rMatLightList.at(i)->matImageGray.convertTo(currentMatImage, CV_32F);

        cv::Mat currentMatImageN = currentMatImage * normFactor;
        cv::Point shift = calculateSADShift(refMatN, currentMatImageN, cvRectROI, 50);
        std::cout << "Shifts = " << shift << std::endl;

        cv::Mat shiftedMat = shiftImage(rMatLightList.at(i), shift);
        resultList << new RMat(shiftedMat, false, rMatLightList.at(i)->getInstrument());
    }

}

void RProcessing::registerSeriesCustomPropagate()
{
    // Purpose:
    // Originally, in registerSeriesCustom(), I used the first image in the timeline as a reference
    // During the total eclipse, statistical properties change rapidly, and normalization by exposure time is not enough
    // (it's also the case with clouds, fire haze passing quickly etc...)
    // Thus I should opt for a propagating approach, and make a by-pair coalignment, and propagate the shift with respect to the first image

    std::cout << "Starting registerSeriesCustomPropagate()... " << std::endl;

    if (cvRectROIList.empty())
    {
        emit tempMessageSignal(QString("ROI not defined"));
        return;
    }

    if(applyMask & (maskCircleRadius ==0))
    {
        emit tempMessageSignal(QString("circular mask not defined"));
        return;
    }

    if (!treeWidget->rMatLightList.isEmpty())
    {
        rMatLightList = treeWidget->rMatLightList;
    }
    else
    {
        emit tempMessageSignal(QString("No lights to register"));
        return;
    }
    // Clear the resultList if not empty
    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    std::cout << "registerSeriesCustomPropagate() check passed. " << std::endl;

    cv::Mat refMat;
    // Normalized version
    cv::Mat refMatN;
    cv::Mat currentMatImage;
    cv::Mat currentMatImageN;


    if (rMatLightList.at(0)->isBayer())
    {
        resultList << new RMat(rMatLightList.at(0)->matImageRGB, false, rMatLightList.at(0)->getInstrument());
    }
    else
    {
        resultList << new RMat(rMatLightList.at(0)->matImage, false, rMatLightList.at(0)->getInstrument());
    }

    float normFactor = 1.0f / rMatLightList.at(0)->getXPOSURE();

    std::cout << "cvRectROI = " << cvRectROI << std::endl;
    std::cout << "cvRectROIList(0) = " << cvRectROIList.at(0) << std::endl;

    cv::Point shift(0, 0);
    for (int i=0; i < rMatLightList.size()-1; i++)
    {
        rMatLightList.at(i)->matImageGray.convertTo(refMat, CV_32F);
        refMatN = refMat * normFactor;

        rMatLightList.at(i+1)->matImageGray.convertTo(currentMatImage, CV_32F);

        currentMatImageN = currentMatImage * normFactor;
        //shift = shift + calculateSADShift(refMatN, currentMatImageN, cvRectROI, 50);
        shift = shift + calculateSADShift(refMatN, currentMatImageN, cvRectROIList, 50);
        std::cout << "Shifts = " << shift << std::endl;

        cv::Mat shiftedMat = shiftImage(rMatLightList.at(i+1), shift);
        resultList << new RMat(shiftedMat, false, rMatLightList.at(i)->getInstrument());
    }
}

cv::Point RProcessing::calculateSADShift(cv::Mat refMat, cv::Mat matImage, cv::Rect fov, int maxLength)
{

    cv::Mat mask;
    if(applyMask)
    {
        mask = circleMask(refMat, maskCircleX, maskCircleY, maskCircleRadius);
        //emit resultSignal(refMat, false, instruments::DSLR);
    }
    else
    {
        mask = cv::Mat::ones(refMat.size(), CV_8U);
    }

    cv::GaussianBlur(refMat, refMat, cv::Size(0, 0), 3);
    cv::GaussianBlur(matImage, matImage, cv::Size(0, 0), 3);

    cv::Mat diffMat;
    cv::Mat absDiffMat;
    cv::Mat sMask = mask(fov);
    cv::Mat sMask2;
    cv::Mat sRefMat = refMat(fov);
    cv::Mat sMatImage;
    cv::Mat SAD = cv::Mat::zeros(maxLength, maxLength, CV_32F);
    cv::Rect fov2 = fov;

    // Loop through all test shifts using sum of absolute difference
    // Beware, SAD is initially populated with zeros. So I must not "miss a spot".
    for (int i=0; i < maxLength; ++i)
    {
        int dy = -maxLength / 2 + i;
        fov2.y = fov.y + dy;

        for (int j=0; j < maxLength; ++j)
        {
            int dx = -maxLength/2 + j;
            fov2.x = fov.x + dx;
            //std::cout << "[dx, dy] " << "[" << dx << "," << dy << "]" << std::endl;

            sMatImage = matImage(fov2);
            // Combine the shifted mask and the reference mask into one mask, with zeros wherever one of them is zeros
            cv::bitwise_and(sMask, mask(fov2), sMask2);
            cv::subtract(sRefMat, sMatImage, diffMat, sMask2);
            absDiffMat = cv::abs(diffMat);
            cv::Scalar sad = cv::sum(absDiffMat) / cv::sum(sMask2);
            SAD.at<float>(i,j) = (float) sad.val[0];
        }
    }

    cv::Point minLoc;
    cv::Point maxLoc;
    double minValue;
    double maxValue;
    cv::minMaxLoc(SAD, &minValue, &maxValue, &minLoc, &maxLoc);
    cv::Point shift;
    shift.x = minLoc.x - maxLength/2;
    shift.y = minLoc.y - maxLength/2;

    return shift;
}

cv::Point RProcessing::calculateSADShift(cv::Mat refMat, cv::Mat matImage, QList<cv::Rect> fovList, int maxLength)
{
    cv::Point shift;
    for (int i=0; i < fovList.size(); i++)
    {
        cv::Rect fov = fovList.at(i);
        shift += calculateSADShift(refMat, matImage, fov, maxLength);
    }

    // Average the shifts
    shift /= fovList.size();
    std::cout << "calculateSADShift:: shift = "<< std::endl;
    std::cout << shift << std::endl;

    return shift;
}

cv::Mat RProcessing::calculateXCorrShift(cv::Mat refMat, cv::Mat matImage, cv::Rect fov)
{
    cv::GaussianBlur(refMat, refMat, cv::Size(0, 0), 3);
    cv::GaussianBlur(matImage, matImage, cv::Size(0, 0), 3);

    cv::Mat sRefMat = refMat(fov);
    cv::Mat sMatImage = matImage(fov);

    double eccEps = 0;
    // Define the motion mode
    const int warpMode = cv::MOTION_TRANSLATION;
    // Specify the number of iterations.
    //int number_of_iterations = 50; // stellar
    int number_of_iterations = 1000; // stellar
    // Specify the threshold of the increment
    // in the correlation coefficient between two iterations
    double termination_eps = 1e-2; // stellar


    // Define termination criteria
    cv::TermCriteria criteria(cv::TermCriteria::COUNT, number_of_iterations, termination_eps);
    //cv::TermCriteria criteria(cv::TermCriteria::COUNT+cv::TermCriteria::EPS, number_of_iterations, termination_eps);

    cv::Mat warpMatrix = cv::Mat::eye(2, 3, CV_32F);

    eccEps = cv::findTransformECC(
                sRefMat,
                sMatImage,
                warpMatrix,
                warpMode,
                criteria
                );

    return warpMatrix;
}


cv::Mat RProcessing::calculateXCorrShift(cv::Mat refMat, cv::Mat matImage, QList<cv::Rect> fovList)
{

    cv::Mat warpMatrixTotal = cv::Mat::eye( 2, 3, CV_32FC1 );

    for (int i=0; i < fovList.size(); i++)
    {
        cv::Rect fov = fovList.at(i);
        cv::Mat warpMatrix_i = calculateXCorrShift(refMat, matImage, fov);

        // Translation in 1st dimension (rows? y-axis?)
        warpMatrixTotal.at<float>(0,2) += warpMatrix_i.at<float>(0, 2);
        // Translation in 2nd dimension (cols? x-axis?)
        warpMatrixTotal.at<float>(1,2) += warpMatrix_i.at<float>(1, 2);
    }

    // Average the shifts
    warpMatrixTotal.at<float>(0,2) /= fovList.size();
    warpMatrixTotal.at<float>(1,2) /= fovList.size();

    std::cout << "RProcessing::calculateXCorrShift:: shift X = " << warpMatrixTotal.at<float>(0,2) << std::endl;
    std::cout << "RProcessing::calculateXCorrShift:: shift Y = " << warpMatrixTotal.at<float>(1,2) << std::endl;

    return warpMatrixTotal;

}

cv::Mat RProcessing::shiftToWarp(cv::Point shift)
{
    cv::Mat warpMat = cv::Mat::eye(2, 3, CV_32F);
    warpMat.at<float>(0, 2) = shift.x;
    warpMat.at<float>(1, 2) = shift.y;

    return warpMat;
}

void RProcessing::registerSeriesByTemplateMatching()
{
    /// The reason of trying this method in addition to the X-correlation rose from saturated image during the eclipse.
    /// The longer the exposure, the farther the off-limb coronal features saturate.
    /// Thus it becomes necessary to exclude a larger saturated part. Experience shows that enhanced X-Correlation fails
    /// when the ROI(s) become(s) small, and the other similarity metrics offered by template matching offer alternatives.

    /// Work with a single ROI for now. Check if it exists.
    if (cvRectROI.empty())
    {
        emit tempMessageSignal(QString("ROI not defined"));
        return;
    }
    resultList << new RMat(rMatLightList.at(0)->matImageRGB, false, rMatLightList.at(0)->getInstrument());
    // Reference image
    cv::Mat refMatN;
    rMatLightList.at(0)->matImageGray.convertTo(refMatN, CV_32F);
    refMatN = refMatN / rMatLightList.at(0)->getXPOSURE();

    for (int i = 1; i < rMatLightList.size(); ++i)
    {
        cv::Mat currentMatImageN;
        rMatLightList.at(i)->matImageGray.convertTo(currentMatImageN, CV_32F);
        currentMatImageN = currentMatImageN / rMatLightList.at(i)->getXPOSURE();

        cv::Mat warpMat = calculateTemplateMatchShift(refMatN, currentMatImageN, cvRectROI);
        cv::Mat shiftedMat = shiftImage(rMatLightList.at(i), warpMat);
        resultList << new RMat(shiftedMat, false, rMatLightList.at(i)->getInstrument());
    }
}

void RProcessing::registerSeriesByTemplateMatchingPropagate()
{
    /// Work if a single ROI for now. Check that it exists
    if (cvRectROI.empty())
    {
        emit tempMessageSignal(QString("ROI not defined"));
        return;
    }

    if (rMatLightList.at(0)->isBayer())
    {
        resultList << new RMat(rMatLightList.at(0)->matImageRGB, false, rMatLightList.at(0)->getInstrument(), rMatLightList.at(0)->getXPOSURE(), rMatLightList.at(0)->getTEMP());
        resultList.at(0)->setFileInfo(rMatLightList.at(0)->getFileInfo());
    }
    else
    {
        resultList << new RMat(rMatLightList.at(0)->matImage, false, rMatLightList.at(0)->getInstrument(), rMatLightList.at(0)->getXPOSURE(), rMatLightList.at(0)->getTEMP());
        resultList.at(0)->setFileInfo(rMatLightList.at(0)->getFileInfo());
    }

    cv::Mat refMat;
    cv::Mat currentMatImage;
    cv::Mat warpMatrixTotal = cv::Mat::eye( 2, 3, CV_32FC1 );

    for (int i=0; i < rMatLightList.size()-1; i++)
    {
        rMatLightList.at(i)->matImageGray.convertTo(refMat, CV_32F);
        rMatLightList.at(i+1)->matImageGray.convertTo(currentMatImage, CV_32F);

        cv::Mat refMatN = refMat / rMatLightList.at(i)->getXPOSURE();
        cv::Mat currentMatImageN = currentMatImage / rMatLightList.at(i+1)->getXPOSURE();

        cv::Mat warpMat = calculateTemplateMatchShift(refMatN, currentMatImageN, cvRectROI);
        warpMatrixTotal.at<float>(0, 2) += warpMat.at<float>(0, 2);
        warpMatrixTotal.at<float>(1, 2) += warpMat.at<float>(1, 2);

        cv::Mat shiftedMat = shiftImage(rMatLightList.at(i+1), warpMatrixTotal);
        resultList << new RMat(shiftedMat, false, rMatLightList.at(i+1)->getInstrument(), rMatLightList.at(i+1)->getXPOSURE(), rMatLightList.at(i+1)->getTEMP());
        resultList.at(i+1)->setFileInfo(rMatLightList.at(i+1)->getFileInfo());
    }

}



cv::Point RProcessing::templateMatch(cv::Mat img, cv::Mat templ, int matchMethod)
{
    cv::Mat result;
    /// Create the result matrix
    int result_cols =  img.cols - templ.cols + 1;
    int result_rows = img.rows - templ.rows + 1;

    result.create( result_rows, result_cols, CV_32FC1 );

    /// Do the Matching and Normalize
    cv::matchTemplate( img, templ, result, matchMethod );
    cv::normalize( result, result, 0, 1, cv::NORM_MINMAX, -1, cv::Mat() );
    /// Localizing the best match with minMaxLoc
    double minVal; double maxVal; cv::Point minLoc; cv::Point maxLoc;
    cv::Point matchLoc;

    cv::minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, cv::Mat() );

    /// For SQDIFF and SQDIFF_NORMED, the best matches are lower values. For all the other methods, the higher the better
    if( matchMethod  == CV_TM_SQDIFF || matchMethod == CV_TM_SQDIFF_NORMED )
    { matchLoc = minLoc; }
    else
    { matchLoc = maxLoc; }

    return matchLoc;
}

cv::Mat RProcessing::calculateTemplateMatchShift(cv::Mat refMat, cv::Mat matImage, cv::Rect fov)
{
    // Template image. Extract patch with the cv::Rect fov.
    // Testing with the reference image itself
    cv::Mat templ = matImage(fov);
    templ.convertTo(templ, CV_32F);

    cv::Point matchLoc = templateMatch(refMat, templ, CV_TM_SQDIFF);

//        std::cout << "cvRectROI = " << cvRectROI << std::endl;
//        std::cout << "cv::Point matchLoc = " << matchLoc << std::endl;

    // Convert that location into a shift with respect to the original image
    // If the current has moved by a, the algorithm gives a shift of -a.
    // So we need to invert the result to know by how much the current image is shifted with respect to the reference image
    // Finally, by using WARP_INVERSE when warping the image in shiftImage(), the shift given need not to be inverted again.
    cv::Point shift;
    shift.x = -(matchLoc.x - cvRectROI.x);
    shift.y = -(matchLoc.y - cvRectROI.y);
    std::cout << "registerSeriesByTemplateMatching():: shift = " << shift << std::endl;

    // Shift the image
    cv::Mat warpMatrix = shiftToWarp(shift);

    return warpMatrix;
}




cv::Mat RProcessing::shiftImage(RMat * rMatImage, cv::Mat warpMat)
{
    cv::Mat registeredMat;
    if (rMatImage->isBayer())
    {
        // RGB array for splitting channels
        cv::Mat tempMatRGB[3];
        // Split the channels
        cv::split(rMatImage->matImageRGB, tempMatRGB);
        cv::warpAffine(tempMatRGB[0], tempMatRGB[0], warpMat, tempMatRGB[0].size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
        cv::warpAffine(tempMatRGB[1], tempMatRGB[1], warpMat, tempMatRGB[0].size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
        cv::warpAffine(tempMatRGB[2], tempMatRGB[2], warpMat, tempMatRGB[0].size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);

        tempMatRGB[0].convertTo(tempMatRGB[0], CV_16U);
        tempMatRGB[1].convertTo(tempMatRGB[1], CV_16U);
        tempMatRGB[2].convertTo(tempMatRGB[2], CV_16U);
        vector<cv::Mat> channels;
        // Watch the order. 1st array will be last channel.
        channels.push_back(tempMatRGB[0]);
        channels.push_back(tempMatRGB[1]);
        channels.push_back(tempMatRGB[2]);

        cv::merge(channels, registeredMat);
        // registeredMat is necessarily non-bayer.
    }
    else
    {
        cv::warpAffine(rMatImage->matImage, registeredMat, warpMat, rMatImage->matImage.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
        // registeredMat is necessarily non-bayer.
    }
    return registeredMat;
}

cv::Mat RProcessing::shiftImage(RMat *rMatImage, cv::Point shift)
{
    cv::Mat warpMat = cv::Mat::eye(2, 3, CV_32F);
    warpMat.at<float>(0, 2) = shift.x;
    warpMat.at<float>(1, 2) = shift.y;
    cv::Mat shiftedMat = shiftImage(rMatImage, warpMat);
    return shiftedMat;
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
        ImageManager *newImageManager = new ImageManager(treeWidget->getLightUrls().at(i));
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

        int numDots = 128;

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
        cv::Mat matDistances(distances, false);
        cv::Scalar mean, stddev;
        cv::meanStdDev(matDistances, mean, stddev);
        double meand = mean.val[0];
        qDebug("mean = %f", meand);

        float median = calcMedian(distances, 0.1);

        std::vector<cv::Point> newPoints;
        for (int j = 0; j < wernerPoints.n; j++)
        {
//            if ( abs(distances[j] - mean[0] ) <  stddev[0] )
            if ( abs(distances[j] - median ) <  stddev.val[0] )
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
        qDebug("[xc2, yc2, Rm2] = [%.2f, %.2f, %.2f]", circleOut2.a, circleOut2.b, circleOut2.r);

        /// 3rd pass
        cv::Point2f circleCenter2(circleOut2.a, circleOut2.b);
        std::vector<float> distances2(cleanDataPoints.n);
        for (int j = 0; j < cleanDataPoints.n; j++)
        {
            cv::Point2f point( (float) cleanDataPoints.X[j], (float) cleanDataPoints.Y[j]);
            float distance = cv::norm(point-circleCenter2);
            distances2[j] = distance;
        }
        cv::Mat matDistances2(distances2, false);
        cv::meanStdDev(matDistances2, mean, stddev);
        median = calcMedian(distances2, 0.1);

        std::vector<cv::Point> newPoints2;
        for (int j = 0; j < cleanDataPoints.n; j++)
        {
            //if ( abs(distances2[j] - mean[0] ) <  stddev[0] )
            if ( abs(distances2[j] - median ) <  stddev.val[0] )
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
            if (newPoints2.empty())
            {
                // draw the 1st fitted circle in green
                cv::circle(contoursMat, circleCenter, circleOut.r, green, 2, 8);
            }
            else
            {
                // draw the 3rd fitted circle in orange
                cv::Point2f circleCenter3(circleOut3.a, circleOut3.b);
                cv::circle(contoursMat, circleCenter3, circleOut3.r, red, 1, 8);
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
    cv::Mat matRadii(radii, false);
    cv::Scalar tempRadius = cv::mean(matRadii);
    meanRadius = (float) tempRadius.val[0];
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
    cv::Mat matImageF;
    matImage.convertTo(matImageF, CV_32F);
    cv::Mat matSlice;
    cv::Mat gradXLeft, gradXRight, gradYBottom, gradYTop;
    cv::Mat gradX, gradY;
    cv::Point minLoc, maxLoc;
    double minVal, maxVal;
    int naxis1 = matImage.cols;
    int naxis2 = matImage.rows;
    // 1D smoothing kernels
    cv::Mat smoothKerX = cv::Mat::ones(1, 3, CV_32F);
    smoothKerX /= smoothSize;
    cv::Mat smoothKerY = cv::Mat::ones(3, 1, CV_32F);
    smoothKerY /= smoothSize;
    /// Kernels for central derivatives
    cv::Mat kernelXCentDeriv = (cv::Mat_<float>(1,3)<<-0.5, 0, 0.5);
    /// Over y-axis, forward and backward direction
    cv::Mat kernelYCentDeriv = (cv::Mat_<float>(3,1)<<-0.5, 0, 0.5);
    //cv::Size kernelSize(smoothSize, smoothSize);
    if (smoothSize == 0) { smooth = false; }

    for (int ii = 0; ii < numDots; ii++)
    {
        int X = naxis1/4 + ii*naxis1/(2*numDots);
        int Y = naxis2/4 + ii*naxis2/(2*numDots);

        /// Left-hand slices (no copy)
        matSlice = matImageF.rowRange(Y, Y+1);

        if (ii == 0)
        {
            std::cout << "x , y =" << X << ", " << Y << std::endl;
            std::cout << "matSlice = " << std::endl;
            std::cout << matSlice(cv::Range(0, 1), cv::Range(0, 10)) << std::endl;
        }

        if (smooth)
        {
            cv::filter2D(matSlice, matSlice, -1, smoothKerX, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        }

        if (ii == 0)
        {
            std::cout << "matSlice (smoothed)= " << std::endl;
            std::cout << matSlice(cv::Range(0, 1), cv::Range(0, 10)) << std::endl;
        }

        cv::filter2D(matSlice, gradX, -1, kernelXCentDeriv, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        gradXLeft = cv::abs(gradX(cv::Range(0, 1), cv::Range(0, naxis1/4)));

        if (ii == 0)
        {
            std::cout << "gradXLeft = " << endl;
            std::cout << gradXLeft(cv::Range(0, 1), cv::Range(0, 10)) << std::endl;
        }


        cv::minMaxLoc(gradXLeft, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii] = maxLoc.x;
        dat->Y[ii] = Y;

        /// Right-hand slices (no copy)
        gradXRight = cv::abs(gradX(cv::Range(0, 1), cv::Range(3*naxis1/4, naxis1)));
        cv::minMaxLoc(gradXRight, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + numDots] = maxLoc.x + 3*naxis1/4;
        dat->Y[ii + numDots] = Y;

        ///Bottom slices (no copy)
        matSlice = matImageF.colRange(X, X+1);

        if (smooth)
        {
            cv::filter2D(matSlice, matSlice, -1, smoothKerX, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        }
        cv::filter2D(matSlice, gradY, -1, kernelYCentDeriv, cv::Point(-1, -1), 0, cv::BORDER_REPLICATE);
        gradYBottom = cv::abs(gradY(cv::Range(0, naxis2/4), cv::Range(0, 1)));
        cv::minMaxLoc(gradYBottom, &minVal, &maxVal, &minLoc, &maxLoc);
        dat->X[ii + 2*numDots] = X;
        dat->Y[ii + 2*numDots] = maxLoc.y;

        ///Top slices (no copy)
        gradYTop = cv::abs(gradY(cv::Range(3*naxis2/4, naxis2), cv::Range(0, 1)));
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
        normalizedRMatImageList.at(i)->setInstrument(rMatImageList.at(i)->getInstrument());
        normalizedRMatImageList.at(i)->setSOLAR_R(rMatImageList.at(i)->getSOLAR_R());
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

RMat *RProcessing::normalizeToXposure(RMat *rMat)
{
    cv::Mat normalizedMatImage;
    rMat->matImage.convertTo(normalizedMatImage, CV_32F);
    normalizedMatImage = normalizedMatImage / rMat->getXPOSURE();

    RMat* normalizedRMat = new RMat(normalizedMatImage, rMat->isBayer());
    normalizedRMat->setImageTitle(QString("Normalized image "));
    normalizedRMat->setDate_time(rMat->getDate_time());
    normalizedRMat->setSOLAR_R(rMat->getSOLAR_R());
    normalizedRMat->flipUD = rMat->flipUD;

    return normalizedRMat;
}

QList<RMat *> RProcessing::normalizeSeriesToXposure(QList<RMat *> rMatImageList)
{
    QList<RMat*> normalizedRMatImageList;

    for (int i =0 ; i < rMatImageList.size() ; ++i)
    {
        normalizedRMatImageList << normalizeToXposure(rMatImageList.at(i));
        //normalizedRMatImageList.at(i)->setInstrument(rMatImageList.at(i)->getInstrument());
        normalizedRMatImageList.at(i)->setSOLAR_R(rMatImageList.at(i)->getSOLAR_R());
        normalizedRMatImageList.at(i)->setImageTitle(normalizedRMatImageList.at(i)->getImageTitle() + QString("# %1").arg(i));
        normalizedRMatImageList.at(i)->setDate_time(rMatImageList.at(i)->getDate_time());

    }

    return normalizedRMatImageList;
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
    cv::threshold(normalizedMatImage, normalizedMatImage, newRange, newRange, cv::THRESH_TRUNC);

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

    if (matImage.type() == CV_16U)
    {
        for (int x = 0; x < matImage.cols; x++)
        {   /// 1st rows
            matImage.at<ushort>(0, x) = matImage.at<ushort>(1, x);
        }
        for (int y = 0; y < matImage.rows; y++)
        {   /// 1st rows
            matImage.at<ushort>(y, 0) = matImage.at<ushort>(y, 1);
        }
    }
    else if (matImage.type() == CV_32F)
    {
        for (int x = 0; x < matImage.cols; x++)
        {   /// 1st rows
            matImage.at<float>(0, x) = matImage.at<float>(1, x);
        }
        for (int y = 0; y < matImage.rows; y++)
        {   /// 1st rows
            matImage.at<float>(y, 0) = matImage.at<float>(y, 1);
        }
    }

}

void RProcessing::setupMaskingCircle(int circleX, int circleY, int radius)
{
    this->maskCircleX = circleX;
    this->maskCircleY = circleY;
    this->maskCircleRadius = radius;
}

void RProcessing::clearROIs()
{
    cvRectROIList.clear();
}

void RProcessing::appendROIList(QRect qRect)
{
    cv::Rect cvRect(qRect.x(), qRect.y(), qRect.width(), qRect.height());
    cvRectROIList << cvRect;
    for (int i=0; i < cvRectROIList.size(); i++)
    {
        std::cout << "RProcessing::appendROIList  cvRect = " << cvRectROIList.at(i) << std::endl;
    }

}

cv::Mat RProcessing::circleMaskMat(cv::Mat matImage, int circleX, int circleY, int radius)
{
    cv::Mat mask = cv::Mat::ones(matImage.size(), CV_8U);
    cv::Point circleCenter(circleX, circleY);
    cv::circle(mask, circleCenter, radius, cv::Scalar::all(0), -1);
    cv::Mat maskedMat;
    matImage.copyTo(maskedMat, mask);
    return maskedMat;
}

cv::Mat RProcessing::circleMask(cv::Mat matImage, int circleX, int circleY, int radius)
{
    cv::Mat mask = cv::Mat::ones(matImage.size(), CV_8U);
    cv::Point circleCenter(circleX, circleY);
    cv::circle(mask, circleCenter, radius, cv::Scalar::all(0), -1);
    return mask;
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
    masterFlatN->matImage = masterFlatN->matImage / meanValueF;
    masterFlatN->calcMinMax();
    masterFlatN->setImageTitle(QString("Master Flat normalized"));
    float bscale = (float) masterFlat->getDataMax() / meanValueF;
    masterFlatN->setBscale(bscale);
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

void RProcessing::setBinning(int binning)
{
    this->binning = binning;
}

void RProcessing::setBlkSize(int blkSize)
{
    this->blkSize = blkSize;
}

void RProcessing::setNBest(int nBest)
{
    this->nBest = nBest;
}

void RProcessing::setQualityMetric(QString qualityMetric)
{
    this->qualityMetric = qualityMetric;
}

void RProcessing::setApplyMask(bool status)
{
    this->applyMask = status;
}

void RProcessing::setCvRectROIList(QList<cv::Rect> cvRectList)
{
    this->cvRectROIList = cvRectList;
}

void RProcessing::setMaskCircleX(int circleX)
{
    this->maskCircleX = circleX;
}

void RProcessing::setMaskCircleY(int circleY)
{
    this->maskCircleY = circleY;
}

void RProcessing::setMaskCircleRadius(int circleRadius)
{
    this->maskCircleRadius = circleRadius;
}

void RProcessing::setUseUrlsFromTreeWidget(bool status)
{
    this->useUrlsFromTreeWidget = status;
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
    qDebug("ar.dims() = [%d, %d, %d, %d]", ar.dims(0), ar.dims(1), ar.dims(2), ar.dims(3));
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




void RProcessing::makeAlignedStack2(af::array &stackedBlks, const af::array &arfSeries, const af::array &qualityBinnedSeries, const af::array &arDim, const af::array &xRange, const af::array &yRange, const int nBest, const int &blkSize, const int &binnedBlkSize, const int &binning, int &x, int &y)
{
    int xB = x/binning;
    int yB = y/binning;

    af::array binnedCube = qualityBinnedSeries(af::seq(xB, xB + binnedBlkSize -1), af::seq(yB, yB + binnedBlkSize -1), af::span);
    af::array gradSum = flat(sum(sum(binnedCube, 0), 1));

    /// Sort the sum of the gradient-norm. Sorted Array is the array ordered in decreasing quality
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, gradSum, 0, false);
    af::array bestInds = sortIndices(af::seq(0, nBest-1));
    unsigned int *bestIndsH = bestInds.host<unsigned int>();

    stackedBlks = arfSeries(xRange, yRange, bestInds);
    af::array refBlk = stackedBlks.slice(0);
    af::array shifts = af::constant(0, 2, f32);

    qDebug("makeAlignedStack2:");
    for (int i = 1; i < nBest; i++)
    {
        af::array blk = stackedBlks.slice(i);
        phaseCorrelate(refBlk, blk, arDim, shifts);
        af::array xShift = tile(shifts(0), blkSize);
        af::array yShift = tile(shifts(1), blkSize);
        qDebug("Shift at Frame %d:", i);
        af_print(shifts);

        /// Need to check these values... some of them are not physically accurate
        /// although they are precise in a correlation sense.
//        af_print(xShift);

        af::array xr = xRange - xShift;
        af::array yr = yRange - yShift;
        stackedBlks(af::span, af::span, i) = arfSeries(xr, yr, bestIndsH[i]);
    }

}

void RProcessing::makeAlignedStackGradient(af::array &stackedBlks, const af::array &arfSeries, const af::array &qualityBinnedSeries, const af::array & arDim, const af::array &xRange, const af::array &yRange, int &x, int &y)
{
    /// Assumes qualitySeries as the sum of absolute gradient components |dx| + |dy| from blockProcessingGradient()
    /// Need to check if valid for blockProcessingSobel() as well...

    int xB = x/binning;
    int yB = y/binning;
    int binnedBlkSize = blkSize/binning;

    af::array qualityBlk = qualityBinnedSeries(af::seq(xB, xB + binnedBlkSize -1), af::seq(yB, yB + binnedBlkSize -1), af::span);
    af::array qualityBlk2 = af::moddims(qualityBlk, qualityBlk.dims(0)*qualityBlk.dims(1), qualityBlk.dims(2));
    // Quality set to total sum of (|dx| + |dy|)
    //af::array quality = af::flat(af::sum(qualityBlk2, 0));

    // Quality set to max of (|dx| + |dy|)
    af::array quality = af::flat(af::max(qualityBlk2, 0));

    /// Sort the sum of the gradient-norm. Sorted Array is the array ordered in decreasing quality
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, quality, 0, false);
    af::array bestInds = sortIndices(af::seq(0, nBest-1));
    unsigned int *bestIndsH = bestInds.host<unsigned int>();

    stackedBlks = arfSeries(xRange, yRange, bestInds);
//    std::cout << "Block values" << std::endl;
//    af_print(stackedBlks(af::seq(0, 9), af::seq(0,9), 1))

    if (nBest == 1)
    {
        return;
    }

    af::array shifts = af::constant(0, 2, nBest-1);
    phaseCorrelate2(stackedBlks, shifts, arDim);


    for (int i = 1; i < nBest; i++)
    {
        af::array xShift = tile(shifts(0, i-1), blkSize);
        af::array yShift = tile(shifts(1, i-1), blkSize);
        /// Need to check these values... some of them are not physically accurate
        /// although they are precise in a correlation sense.

        af::array xr = xRange - xShift;
        af::array yr = yRange - yShift;
        stackedBlks(af::span, af::span, i) = arfSeries(xr, yr, bestIndsH[i]);
    }
}

void RProcessing::makeAlignedStackLaplace(af::array &stackedBlks, const af::array &arfSeries, const af::array &qualitySeries, const af::array &arDim, const af::array &xRange, const af::array &yRange, int &x, int &y)
{
    /// Same as makeAlignedStackGradient but assumes the Laplacian metric in the qualitySeries.
    int xB = x/binning;
    int yB = y/binning;
    int binnedBlkSize = blkSize/binning;

    af::array qualityBlk = qualitySeries(af::seq(xB, xB + binnedBlkSize -1), af::seq(yB, yB + binnedBlkSize -1), af::span);
    af::array qualityBlk2 = af::abs(af::moddims(qualityBlk, qualityBlk.dims(0)*qualityBlk.dims(1), qualityBlk.dims(2)));
    af::array quality = af::flat(af::var(qualityBlk2, false, 0));
    //af::array quality = af::flat(af::max(qualityBlk2, 0));
    //af::array quality = af::flat(af::sum(qualityBlk2, 0));

    /// Sort the quality in descending order
    af::array sortedArray;
    af::array sortIndices;
    af::sort(sortedArray, sortIndices, quality, 0, false);
    af::array bestInds = sortIndices(af::seq(0, nBest-1));
    af_print(bestInds);
    af_print(xRange);
    af_print(yRange);
    unsigned int *bestIndsH = bestInds.host<unsigned int>();

    stackedBlks = arfSeries(xRange, yRange, bestInds);

    if (nBest == 1)
    {
        return;
    }

    af::array shifts = af::constant(0, 2, nBest-1);
    phaseCorrelate2(stackedBlks, shifts, arDim);
    //af_print(shifts);

    for (int i = 1; i < nBest; i++)
    {
        af::array xShift = tile(shifts(0, i-1), blkSize);
        af::array yShift = tile(shifts(1, i-1), blkSize);

        af::array xr = xRange - xShift;
        af::array yr = yRange - yShift;
        stackedBlks(af::span, af::span, i) = arfSeries(xr, yr, bestIndsH[i]);
    }
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



