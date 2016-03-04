#include "imagemanager.h"

#include <QUrl>
#include <QFileInfo>

ImageManager::ImageManager(QString filePathQStr) :
    newFitsImage(NULL), fitsSeries(NULL), newRawImage(NULL)
{

    if (filePathQStr.isEmpty())
    {
        qDebug("filePath is empty");
        return;
    }

    fitsExtList = QList<QString>();
    fitsExtList << "fits" << "fts" << "fit";

    rawExtList = QList<QString>();
    rawExtList << "cr2";


    // Instantiate an object that hold information on the user's monitor.

    // Process the list of filenames of the dragged file.
    // Here it processes the name of 1st of the list to know the file type,
    // but then it actually loads the whole lot! This may change in the future.

    this->filePathQStr = filePathQStr;
    QFileInfo fileInfo(filePathQStr);
    fileName = fileInfo.fileName();
    fileExt = fileInfo.suffix().toLower();

    if (fitsExtList.contains(fileExt))
    {
        loadFits();
    }
    else if (rawExtList.contains(fileExt))
    {
        loadRaw();
    }
    else
    {
        qDebug("Unknown file extension");
        return;
    }

    qDebug() << "ImageManager:: dataMin = " << dataMin;
    qDebug() << "ImageManager:: dataMax = " << dataMax;

}

ImageManager::~ImageManager()
{
    delete newFitsImage;
    delete newRawImage;
}


void ImageManager::loadFits()
{
    newFitsImage = new MyFitsImage(filePathQStr);

    cv::minMaxLoc(newFitsImage->matFits, &dataMin, &dataMax);

    instruments instrument;
    // USET data?
    if (newFitsImage->getKeyValues().contains(QString("USET")))
    {
        qDebug("ImageManager::loadFits():: loading USET data");
        instrument = instruments::USET;
    }
    else if (dataMin < -100.0)
    {
        instrument = instruments::MAG;
    }
    else
    {
        instrument = instruments::generic;
    }

    rMatImage = RMat(newFitsImage->matFits, newFitsImage->isBayer(), instrument);
    rMatImage.setBscale(newFitsImage->getBscale());
    rMatImage.setBzero(newFitsImage->getBzero());

    rMatImage.setDataMin((float) dataMin);
    rMatImage.setDataMax((float) dataMax);

    rMatImage.setWbRed(1.0);
    rMatImage.setWbGreen(1.0);
    rMatImage.setWbBlue(1.0);
}

void ImageManager::loadRaw()
{
    newRawImage = new RawImage(filePathQStr);
    cv::minMaxLoc(newRawImage->matCFA, &dataMin, &dataMax);

    rMatImage = RMat(newRawImage->matCFA, true, instruments::DSLR);

    rMatImage.setDataMin((float) dataMin);
    rMatImage.setDataMax((float) dataMax);

    rMatImage.setWbRed(newRawImage->getWbRed());
    rMatImage.setWbGreen(newRawImage->getWbGreen());
    rMatImage.setWbBlue(newRawImage->getWbBlue());
}

MyFitsImage* ImageManager::getNewFitsImage()
{
    return newFitsImage;
}

QString ImageManager::getFileName()
{
    return fileName;
}

float ImageManager::getWbRed() const
{
    return wbRed;
}

float ImageManager::getWbGreen() const
{
    return wbGreen;
}

float ImageManager::getWbBlue() const
{
    return wbBlue;
}
