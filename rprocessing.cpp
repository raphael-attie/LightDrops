#include "rprocessing.h"

#include <QFileDialog>

// cifitsio
#include <cfitsio/fitsio.h>

//opencv
#include <opencv2/world.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/video.hpp>
//#include <opencv2/core/ocl.hpp>

#include "imagemanager.h"
#include "parallelcalibration.h"
#include "typedefs.h"


RProcessing::RProcessing(QObject *parent): QObject(parent),
    masterBias(NULL), masterDark(NULL), masterFlat(NULL), masterFlatN(NULL), useXCorr(false),
    radius(0), radius1(0), radius2(0), radius3(0)
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
        qDebug("RProcessing:: cannyEdgeDetection() on image #%i", i);
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

    emit listResultSignal(resultList);
}

void RProcessing::registerSeries(QList<RMat *> rMatList)
{

    int nFrames = rMatList.size();

    resultList2.reserve(nFrames);

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
    rMatList.at(0)->matImage.convertTo(refMat, CV_16U);


    RMat *refRMat = new RMat(refMat, false, rMatLightList.at(0)->getInstrument());
    resultList2 << refRMat;
    resultList2.at(0)->setImageTitle(QString("X-corr registered image # 1"));

    // Get the new reference image (at i = 0).
    cv::Mat refMat1 = normalizeClipByThresh(rMatList.at(0), 1000.0f, 2000.0f);
    // The output is already converted to CV_32F, so no need to convert
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

        // Get the image to co-align with respect to the reference image
        cv::Mat registeredMat = normalizeClipByThresh(rMatList.at(i), 1000.0f, 2000.0f);

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
                    warp_mode_2,
                    criteria2
                    );

        qDebug() << "eccEps 2 =" << eccEps;
        std::cout << "result warp_matrix 2 =" << std::endl << warp_matrix_1 << std::endl << std::endl;

        cv::warpAffine(rMatList.at(i)->matImage, registeredMat, warp_matrix_1, registeredMat.size(), cv::INTER_CUBIC + CV_WARP_INVERSE_MAP);
        registeredMat.convertTo(registeredMat, CV_16U);
        resultList2 << new RMat(registeredMat, false, rMatLightList.at(i)->getInstrument()); // This RMat is necessarily non-bayer.
        resultList2.at(i)->setImageTitle(QString("X-corr registered image # %1").arg(i));
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

    if (!circleOutList.isEmpty())
    {
        circleOutList.clear();
    }

    radius = 0;

    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        qDebug("RProcessing:: cannyEdgeDetection() on image #%i", i);
        setupCannyDetection(i);
        cannyDetect(thresh);
        limbFit(i);

        /// Get results showing contours of all the edges
        contoursRMat = new RMat(contoursMat.clone(), false);
        contoursRMat->setImageTitle(QString("Canny edges: Image # %1").arg(i));
        contoursRMatList << contoursRMat;

//        qDebug("RProcessing::cannyEdgeDetection():: [centers.at(i).x ; centers.at(i).y] = [%f ; %f]", centers.at(i).x, centers.at(i).y);
//        qDebug("RProcessing::cannyEdgeDetection():: ellRect.boundingRect().center = [%f ; %f]", (float) ellRectList.at(i).boundingRect().width/2.0f, (float) ellRectList.at(i).boundingRect().height/2.0f);
//        qDebug("RProcessing::cannyEdgeDetection():: ellRect.size = [%f ; %f]", ellRectList.at(i).size.width, ellRectList.at(i).size.height);
//        qDebug("RProcessing::cannyEdgeDetection():: ellRect.angle = %f", ellRectList.at(i).angle);
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

}

void RProcessing::setupCannyDetection(int i)
{

    // Normalize the data to have equivalent statistics / histograms
    //normalizeByStats(rMatLightList.at(i));

    cv::Mat normalizedMat = normalizeClipByThresh(rMatLightList.at(i), rMatLightList.at(i)->getIntensityLow(), rMatLightList.at(i)->getIntensityHigh());

    normalizedMat.convertTo(sampleMat8, CV_8U, 256.0f / rMatLightList.at(i)->getDataRange());
    normalizedMat.convertTo(sampleMatN, CV_8U, 256.0f / rMatLightList.at(i)->getDataRange());

    // Setting for contrast stretching to boost contrast between limb and off-limb
    float newMin = 0;
    float newMax = 100;
    float newDataRange = newMax - newMin;
    float alpha = 256.0f / newDataRange;
    float beta = -newMin * 256.0f /newDataRange;

    // Convert to 8 bit with contrast stretching to boost contrast between limb and off-limb
    sampleMatN.convertTo(sampleMatN, CV_8U, alpha, beta);


}

void RProcessing::cannyDetect(int thresh)
{
    /// Detect edges using cannycompareContourAreas
    double thresh1 = ((double) thresh) / 3.0;
    double thresh2 = (double) thresh;

    // Canny edge detection works best when blurring a bit.
    cv::blur(sampleMatN, sampleMatN, cv::Size(3,3));
    cv::Canny(sampleMatN, contoursMat, thresh1, thresh2, 3);

/// void adaptiveThreshold(InputArray src, OutputArray dst, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C)
    //cv::adaptiveThreshold(sampleMatN, contoursMat, 255, CV_ADAPTIVE_THRESH_MEAN_C, CV_THRESH_BINARY, 5, -5);

}

void RProcessing::limbFit(int i)
{
    // Find contours
    vector< vector <cv::Point> > contours;
    vector< cv::Vec4i > hierarchy;
    cv::findContours(contoursMat.clone(), contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_NONE, cv::Point(0, 0));

    cv::cvtColor(contoursMat, contoursMat, CV_GRAY2RGB);

    if (contours.size() == 0)
    {
        qDebug("No contours found at image %i", i+1);
        exit(1);
    }

    // sort contours
    std::sort(contours.begin(), contours.end(), compareContourAreas);
    // grab contours
    //vector< cv::Point > biggestContour = contours[contours.size()-1];
    //std::vector<cv::Point> smallestContour = contours[0];

    // gather points of all contours in one big vector

    vector< vector<cv::Point> > biggestContours;
    size_t nSelectedContours = std::min(contours.size(), (size_t) 20);

    for (size_t ii = 1 ; ii <= nSelectedContours ; ++ii)
    {
        biggestContours.push_back(contours[contours.size() - ii]);
    }

    biggestContoursList << biggestContours;

    size_t nContourPoints = 0;

       for (int ii = 0; ii < biggestContours.size(); ++ii)
       {
            nContourPoints += biggestContours[ii].size();
       }

    vector< cv::Point > contours1D;
    contours1D.reserve(nContourPoints);

    for (int ii = 0; ii < biggestContours.size(); ++ii)
      {
        const vector< cv::Point > & v = biggestContours[ii];
        contours1D.insert( contours1D.end() , v.begin() , v.end() );
      }

    if (showContours)
    {

                cv::Scalar color = cv::Scalar( 0, 255, 0);
//                cv::drawContours( contoursMat, biggestContours, 0, color, 2, 8);

                for( int ii = 0; ii < biggestContours.size(); ++ii )
                {
                    cv::drawContours( contoursMat, biggestContours, ii, color, 2, 8, hierarchy, 0, cv::Point() );
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
//    qDebug("Hyper fit output:");
//    qDebug("Center at [%f ; %f]", circleOut1.a, circleOut1.b);
//    qDebug("Radius R = %f", circleOut1.r);
    radius1 += circleOut1.r;

//    reals circleX = ellRect.center.x;
//    reals circleY = ellRect.center.y;
//    reals circleR = 913.0f;
//    Circle circleInit(circleX, circleY, circleR);
    Circle circleInit = circleOut1;
    //circleInit.r = 917.0f;

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

//    Circle circleInitR(circleX, circleY, circleR);
//    Circle circleInitR = circleOut1;
//    circleInitR.r = circleR;
//    Circle circleOut3;
//    CircleFitByLevenbergMarquardtReduced(contourData, circleInit, lambdaIni, circleOut3);
//    qDebug("LMA-reduced output:");
//    qDebug("status = %i", status);
//    qDebug("Center at [%f ; %f]", circleOut3.a, circleOut3.b);
//    qDebug("Radius R = %f", circleOut3.r);
//    qDebug("");

//    centers.append(cv::Point2f((float) circleOut3.a, (float) circleOut3.b));
    //centers.append(ellRect.center);
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

    if (!limbFitPreviewList.isEmpty())
    {
        qDeleteAll(limbFitPreviewList);
        limbFitPreviewList.clear();
    }



    for (int i = 0 ; i < rMatLightList.size() ; ++i)
    {
        qDebug("Canny-registering image # %i ", i+1);
        /// Register series
        cv::Mat warpMat = cv::Mat::eye( 2, 3, CV_32FC1 );
        cv::Point2f origin(rMatLightList.at(i)->matImage.cols / 2.0f, rMatLightList.at(i)->matImage.rows / 2.0f);

        cv::Point2f delta = origin - centers.at(i);
        warpMat.at<float>(0, 2) = delta.x;
        warpMat.at<float>(1, 2) = delta.y;
        //std::cout << "warpMat = " << std::endl << " " << warpMat << std::endl << std::endl;

        cv::Mat registeredMat = rMatLightList.at(i)->matImage.clone();
        registeredMat.convertTo(registeredMat, CV_32F);
        //cv::warpAffine(registeredMat, registeredMat, warpMat, registeredMat.size(), cv::INTER_LANCZOS4 + CV_WARP_INVERSE_MAP);
        cv::warpAffine(registeredMat, registeredMat, warpMat, registeredMat.size(), cv::INTER_LANCZOS4);
        registeredMat.convertTo(registeredMat, CV_16U);

        RMat *resultMat = new RMat(registeredMat, false, rMatLightList.at(i)->getInstrument());
        resultMat->setImageTitle(QString("Registered image # ") + QString::number(i));
        resultMat->setDate_time(rMatLightList.at(i)->getDate_time());
        resultList << resultMat;

        cv::Mat previewMat = registeredMat;
        if (rMatLightList.at(i)->getInstrument() == instruments::USET)
        {
            registeredMat.convertTo(previewMat, CV_8U, 1.0f/16.0f, 0.0f);

        }
        else
        {
            registeredMat.convertTo(previewMat, CV_8U, 1.0f/255.0f, 0.0f);
        }
        cv::cvtColor(previewMat, previewMat, CV_GRAY2RGB);

        RMat *previewRMat = new RMat(previewMat, false);
        normalizeByStatsInPlace(previewRMat);
        previewRMat->setImageTitle(QString("(Preview) Registered image # ") + QString::number(i+1));
        previewRMat->setDate_time(rMatLightList.at(i)->getDate_time());

//        if (showContours)
//        {
//            cv::drawContours(registeredMat8, biggestContoursList.at(i), 0, cv::Scalar(0, 255, 0), 2, 8, cv::noArray(), 0 , delta);
//        }

        if (showLimb)
        {
            Circle circle = circleOutList.at(i);
            cv::Point2f circleCenter(circle.a, circle.b);
            circleCenter += delta;
            cv::Scalar red = cv::Scalar(255, 0, 0);
            cv::circle(previewRMat->matImage, circleCenter, circle.r, red, 2, 8);


//            ellRectList[i].center = origin;
//            cv::ellipse(registeredMat8, ellRectList.at(i), cv::Scalar(255, 0, 0), 2, 8);
        }

        /// void putText(Mat& img, const string& text, Point org, int fontFace, double fontScale, Scalar color, int thickness=1, int lineType=8, bool bottomLeftOrigin=false )
        //putText(registeredMat8, date_time, cv::Point(0, 20), cv::FONT_HERSHEY_SIMPLEX, 1, cv::Scalar(255, 255, 255), true);


        limbFitPreviewList << previewRMat;

    }

    if (useXCorr)
    {
        qDebug("Canny registration:: refinement with cross-correlation");
        registerSeries(resultList);
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

QList<RMat *> RProcessing::normalizeSeriesByStats(QList<RMat*> rMatImageList)
{
    QList<RMat*> normalizedRMatImageList;

    for (int i =0 ; i < rMatImageList.size() ; ++i)
    {
        normalizedRMatImageList << normalizeByStats(rMatImageList.at(i));
        normalizedRMatImageList.at(i)->setImageTitle(normalizedRMatImageList.at(i)->getImageTitle() + QString("# %1").arg(i));
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
    /// This is like normalizeClipByThresh() using newMin and newMax as intensityLow and
    /// intensityHigh from RMat::calcStats()

    cv::Mat matImage;

    rMat->matImage.convertTo(matImage, CV_32F);
    qDebug("RProcessing::normalizeByStats()   rMat->getDataRange() = %f", rMat->getDataRange());

    float newMin = rMat->getIntensityLow();
    float newMax = rMat->getIntensityHigh();
    float newDataRange = newMax - newMin;
    float alpha = rMat->getDataRange() / newDataRange;
    float beta = -newMin * rMat->getDataRange() /newDataRange;

    matImage.convertTo(matImage, rMat->matImage.type(), alpha, beta);

    RMat* normalizeRMat = new RMat(matImage, rMat->isBayer(), rMat->getInstrument());
    normalizeRMat->setImageTitle(QString("Normalized image "));
    qDebug("RProcessing::normalizeByStats()   [newMin , newMax] = [%f , %f]", newMin, newMax);
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



cv::Mat RProcessing::normalizeClipByThresh(RMat *rMat, float newMin, float newMax)
{
    /// This function do contrast stretching and clips the intensity between newMin and newMax.

    // Change the series "in situ". No copy.
    cv::Mat matImage = rMat->matImage.clone();
    matImage.convertTo(matImage, CV_32F);

    float newDataRange = newMax - newMin;
    float alpha = rMat->getDataRange() / newDataRange;
    float beta = -newMin * rMat->getDataRange() /newDataRange;
    matImage = matImage * alpha + beta;
    // Now we need to clip the image between the max and min of the extrema of the instrument data type range.
    cv::threshold(matImage, matImage, newMax, newMax, cv::THRESH_TRUNC);
    matImage = newMax - matImage;
    cv::threshold(matImage, matImage, newMax - newMin, newMax - newMin, cv::THRESH_TRUNC);
    matImage = newMax - matImage;

    matImage.convertTo(matImage, rMat->matImage.type());



    return matImage;
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

QList<RMat *> RProcessing::getLimbFitPreviewList()
{
    return limbFitPreviewList;
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
