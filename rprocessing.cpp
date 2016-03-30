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


RProcessing::RProcessing(QObject *parent): QObject(parent),
    masterBias(NULL), masterDark(NULL), masterFlat(NULL), masterFlatN(NULL)
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
    if (!masterDark != NULL && masterDark != masterBias)
    {
        delete masterDark;
    }

    if (!masterFlat != NULL)
    {
        delete masterFlat;
    }

    if (!masterFlatN != NULL)
    {
        delete masterFlatN;
    }

    if (!resultList.empty())
    {
        qDeleteAll(resultList);
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

    qDebug() << "RProcessing::exportToFits() bayer =" << bayer;

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
        if (masterBias != NULL)
        {
            cv::Mat tempMatFlat;
            cv::Mat tempMatBias;
            masterBias->matImage.copyTo(tempMatBias);
            masterFlat->matImage.copyTo(tempMatFlat);

            tempMatFlat.convertTo(tempMatFlat, CV_32F);
            tempMatBias.convertTo(tempMatBias, CV_32F);

            cv::subtract(tempMatFlat, tempMatBias, masterFlat->matImage);
        }

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

    int naxis1 = rMatList.at(0)->matImage.cols;
    int naxis2 = rMatList.at(0)->matImage.rows;

    // Need to create a Mat image that will host the result of the average.
    // It needs to be the same type, and has the number of channels as the Mat images of the series.
    cv::Mat avgImg = cv::Mat::zeros(naxis2, naxis1, rMatList.at(0)->matImage.type());
    avgImg.convertTo(avgImg, CV_32F);

     // The accumulate function in the for-loop below can only work with up to 3 color channels

    for(int i = 0; i < rMatList.size(); i++)
    {
        cv::Mat tempImage;
//        cv::Mat imageClone = rMatList->operator [](i).matImage.clone();

//        imageClone.convertTo(tempImage, CV_32F);

        rMatList.at(i)->matImage.convertTo(tempImage, CV_32F);
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
        qDebug("RProcessing::calibrate():: image #%i", i);
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

    emit listResultSignal(resultList);
}

void RProcessing::cannyEdgeDetection(int thresh)
{
    // Check if data exist
    if (!treeWidget->rMatLightList.isEmpty())
    {
        rMatLightList = treeWidget->rMatLightList;
    }
    else
    {
        emit tempMessageSignal(QString("No lights for Canny edge detection"));
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

    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        setupCannyDetection(i);
        cannyDetect(thresh);

        /// Draw contours of all the edges
        contoursRMat = new RMat(contoursMat.clone(), false);
        contoursRMat->setImageTitle(QString("All canny edges contours"));
        contoursRMatList << contoursRMat;

        qDebug("RProcessing::cannyEdgeDetection():: [centers.at(i).x ; centers.at(i).y] = [%f ; %f]", centers.at(i).x, centers.at(i).y);
        qDebug("RProcessing::cannyEdgeDetection():: ellRect.boundingRect().center = [%f ; %f]", (float) ellRectList.at(i).boundingRect().width/2.0f, (float) ellRectList.at(i).boundingRect().height/2.0f);
        qDebug("RProcessing::cannyEdgeDetection():: ellRect.size = [%f ; %f]", ellRectList.at(i).size.width, ellRectList.at(i).size.height);
        qDebug("RProcessing::cannyEdgeDetection():: ellRect.angle = %f", ellRectList.at(i).angle);

    }

}

void RProcessing::setupCannyDetection(int i)
{

    // Deep-copy of one sample
    cv::Mat sampleMat = rMatLightList.at(i)->matImage.clone();
    sampleMat.convertTo(sampleMat, CV_32F);
    sampleMat = sampleMat / rMatLightList.at(i)->getMean();
    //sampleMat = sampleMat * currentROpenGLWidget->getAlpha() + currentROpenGLWidget->getBeta();
    sampleMatN = sampleMat.clone();


//    cv::threshold(sampleMatN, sampleMatN, 100, 100, cv::THRESH_TRUNC);
//    cv::threshold(sampleMatN, sampleMatN, 60, 60, cv::THRESH_TOZERO);

    cv::normalize(sampleMat, sampleMat8, 0, 255, cv::NORM_MINMAX);
    cv::normalize(sampleMatN, sampleMatN, 0, 255, cv::NORM_MINMAX);

    sampleMatN = 255 - sampleMatN;

//    float newDataRange = 100.0f - 40.0f;
//    float alpha = 255.0f / newDataRange;
//    float beta = -40.0f * 255.0f /newDataRange;

    float newMin = 150;
    float newMax = 200;
    float newDataRange = newMax - newMin;
    float alpha = 255.0f / newDataRange;
    float beta = -newMin * 255.0f /newDataRange;

    // Convert to 8 bit
    sampleMatN.convertTo(sampleMatN, CV_8U, alpha, beta);
    sampleMat8.convertTo(sampleMat8, CV_8U);

//    RMat *rMatN = new RMat(sampleMatN.clone(), false);
//    rMatN->setImageTitle(QString("Normalized-rescaled 8-bit image"));
//    emit resultSignal(rMatN);

//    cv::threshold(sampleMatN, sampleMatN, 225, 225, cv::THRESH_TRUNC);    // for a negative image
//    cv::threshold(sampleMatN, sampleMatN, 175, 175, cv::THRESH_TOZERO); // for a negative image

}

void RProcessing::cannyDetect(int thresh)
{
    qDebug("updating Canny with thresh = %i", thresh);

    /// Detect edges using cannycompareContourAreas
    double thresh1 = ((double) thresh) / 2.0;
    double thresh2 = (double) thresh;

    // Canny edge detection works best when blurring a bit.
    cv::blur(sampleMatN, sampleMatN, cv::Size(3,3));
    cv::Canny(sampleMatN, contoursMat, thresh1, thresh2, 3);
/// void adaptiveThreshold(InputArray src, OutputArray dst, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C)
    //cv::adaptiveThreshold(sampleMatN, cannyMat, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 5, -5);

//    cannyRMat = new RMat(cannyMat, false);
//    cannyRMat->setImageTitle(QString("Canny edges"));
//    emit resultSignal(cannyRMat);


    // Find contours
    std::vector<std::vector<cv::Point> > contours;
    std::vector<cv::Vec4i> hierarchy;
    cv::findContours(contoursMat.clone(), contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE, cv::Point(0, 0));

    cv::cvtColor(contoursMat, contoursMat, CV_GRAY2RGB);

    // sort contours
    std::sort(contours.begin(), contours.end(), compareContourAreas);
    // grab contours
    std::vector<cv::Point> biggestContour = contours[contours.size()-1];
    //std::vector<cv::Point> smallestContour = contours[0];

    std::vector<cv::Point> allContours;
    // gather points of all contours in one big vector
    for (int i = 0 ; i < contours.size() ; ++i)
    {
        for (int j = 0 ; j < contours.at(i).size() ; ++j)
        {
            allContours.push_back(contours.at(i).at(j));
        }

    }

    std::vector< std::vector<cv::Point> > singleBiggestContour;
    singleBiggestContour.push_back(biggestContour);

    singleBiggestContourList << singleBiggestContour;


    if (showContours)
    {

//                cv::Scalar color = cv::Scalar( 0, 255, 0);
//                cv::drawContours( contoursMat, singleBiggestContour, 0, color, 2, 8);
        cv::RNG rng(12345);

        for( int i = 0; i< contours.size(); i++ )
        {
            cv::Scalar color = cv::Scalar( rng.uniform(0, 255), rng.uniform(0,255), rng.uniform(0,255) ); // random rgb value
            //cv::Scalar color = cv::Scalar( 0, 255, 0);
            cv::drawContours( contoursMat, contours, i, color, 2, 8, hierarchy, 0, cv::Point() );
        }
    }


    //cv::RotatedRect ellRect = cv::fitEllipse(biggestContour);
    cv::RotatedRect ellRect = cv::fitEllipse(allContours);

    ellRectList << ellRect;
    centers << ellRect.center;

    ellRect.size.width -= 2;
    ellRect.size.height -= 2;

    if (showLimb)
    {
        cv::Scalar red = cv::Scalar(255, 0, 0);
        cv::ellipse(contoursMat, ellRect, red, 2, 8);
    }


    // Make another ellipse, 2px wider, red.
//    cv::Scalar red = cv::Scalar( 255, 0, 0);
//    ellRect2.size.width += 4;
//    ellRect2.size.height += 4;
//    cv::ellipse(contoursMat, ellRect2, red, 2, 8);

//    cannyQImage = new QImage(contoursMat.data, contoursMat.cols, contoursMat.rows, QImage::Format_RGB888);
//    emit resultQImageSignal(ellQImage);
//    emit ellipseSignal(ellRect);


}

void RProcessing::updateCannyDetection(int thresh)
{
    cannyDetect(thresh);

    /// Wrap the results of cannyDetect in RMat objects.
    /// 1) Create the Canny RMat image, with all the canny edges
//    cannyRMat = new RMat(cannyMat, false);
//    cannyRMat->setImageTitle(QString("Canny image"));

    /// 2) Draw contours of all the edges
    contoursRMat = new RMat(contoursMat, false);
    contoursRMat->setImageTitle(QString("All canny edges contours"));

    /// 3) Create the RMat of of the ellipse fitting the the biggest contours of the canny image.
//    ellipseRMat = new RMat(ellipseMatRGB888, false);
//    ellipseRMat->setImageTitle(QString("Fitted ellipse"));

}

void RProcessing::cannyRegisterSeries()
{
    if (centers.isEmpty())
    {
        tempMessageSignal(QString("Run limb fitting first."));
        return;
    }

    if (!resultList.isEmpty())
    {
        qDeleteAll(resultList);
        resultList.clear();
    }

    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        /// Register series
        cv::Mat warpMat = cv::Mat::eye( 2, 3, CV_32FC1 );
        cv::Point2f origin(rMatLightList.at(i)->matImage.cols / 2.0f, rMatLightList.at(i)->matImage.rows / 2.0f);

        cv::Point2f delta = origin - centers.at(i);
        warpMat.at<float>(0, 2) = delta.x;
        warpMat.at<float>(1, 2) = delta.y;
        std::cout << "warpMat = " << std::endl << " " << warpMat << std::endl << std::endl;

        cv::Mat registeredMat = rMatLightList.at(i)->matImage.clone();
        registeredMat.convertTo(registeredMat, CV_32F);
        //cv::warpAffine(registeredMat, registeredMat, warpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
        cv::warpAffine(registeredMat, registeredMat, warpMat, registeredMat.size(), cv::INTER_LANCZOS4);
        registeredMat = registeredMat / rMatLightList.at(i)->getMean();
        cv::normalize(registeredMat, registeredMat, 0, 255, cv::NORM_MINMAX);
        registeredMat.convertTo(registeredMat, CV_8U);
        cv::cvtColor(registeredMat, registeredMat, CV_GRAY2RGB);

        if (showContours)
        {
            cv::drawContours(registeredMat, singleBiggestContourList.at(i), 0, cv::Scalar(0, 255, 0), 2, 8, cv::noArray(), 0 , delta);
        }

        if (showLimb)
        {
            ellRectList[i].center = origin;
            cv::ellipse(registeredMat, ellRectList.at(i), cv::Scalar(255, 0, 0), 2, 8);
        }

        resultList << new RMat(registeredMat.clone(), false);
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
