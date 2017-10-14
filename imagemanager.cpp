#include "imagemanager.h"

#include <QUrl>
#include <QFileInfo>
#include <QHeaderView>

ImageManager::ImageManager(QUrl url) :
    error(0), url(url), newFitsImage(NULL), fitsSeries(NULL), newRawImage(NULL), rMatImage(NULL), nKeys(1)
{


    if (url.isEmpty())
    {
        qDebug("filePath is empty");
        return;
    }

    this->filePathQStr = url.toLocalFile();

    fitsExtList = QList<QString>();
    fitsExtList << "fits" << "fts" << "fit";

    rawExtList = QList<QString>();
    rawExtList << "cr2" << "CR2";


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
    else if (fileExt == QString("tiff"))
    {
        loadTiff();
    }
    else
    {
        std::cout << "ImageManager:: filePathQStr=" << filePathQStr.toStdString() << std::endl;
        qDebug("Unknown file extension");
        error = 1;
        return;
    }
    rMatImage->setFileInfo(fileInfo);
    rMatImage->setUrl(url);
    createTableWidget();

    rMatImage->calcStats();
    rMatImage->setImageTitle(fileName);



}

ImageManager::~ImageManager()
{
    qDebug("ImageManager::~ImageManager() calling ImageManager destructor");
//    delete tableWidget;
    delete rMatImage;
//    delete newFitsImage;
//    delete newRawImage;

}

QUrl ImageManager::getUrl() const
{
    return url;
}

bool ImageManager::getError() const
{
    return error;
}


void ImageManager::loadFits()
{
    newFitsImage = new MyFitsImage(filePathQStr);

    rMatImage = new RMat(newFitsImage->getMatFits(), newFitsImage->isBayer());
    rMatImage->setBscale(newFitsImage->getBscale());
    rMatImage->setBzero(newFitsImage->getBzero());
    rMatImage->setExpTime(newFitsImage->getExpTime());

    rMatImage->setWbRed(1.0);
    rMatImage->setWbGreen(1.0);
    rMatImage->setWbBlue(1.0);

    instruments instrument;
    // USET data?
    if (newFitsImage->getKeyValues().contains(QString("USET")))
    {
        instrument = instruments::USET;
        //fixUset();
    }
    else if (newFitsImage->getKeyValues().contains(QString("DSLR")))
    {
        instrument = instruments::DSLR;
    }
    else if (rMatImage->getDataMin() < -100.0)
    {
        instrument = instruments::MAG;
    }
    else
    {
        instrument = instruments::generic;
    }
    rMatImage->setInstrument(instrument);

    if (newFitsImage->getKeyNames().contains(QString("SOLAR_R")))
    {
        int keyInd = newFitsImage->getKeyNames().indexOf("SOLAR_R");
        rMatImage->setSOLAR_R(std::stof(newFitsImage->getKeyValues().at(keyInd).toStdString()));
    }

    if (newFitsImage->getKeyNames().contains(QString("EXPTIME")))
    {
        int keyInd = newFitsImage->getKeyNames().indexOf("EXPTIME");
        rMatImage->setXPOSURE(std::stof(newFitsImage->getKeyValues().at(keyInd).toStdString()));
        rMatImage->setExpTime(std::stof(newFitsImage->getKeyValues().at(keyInd).toStdString()));
    }


    if (newFitsImage->getKeyNames().contains(QString("XPOSURE")))
    {
        int keyInd = newFitsImage->getKeyNames().indexOf("XPOSURE");
        rMatImage->setXPOSURE(std::stof(newFitsImage->getKeyValues().at(keyInd).toStdString()));
    }

    if (newFitsImage->getKeyNames().contains(QString("DATE-OBS")))
    {

        int dateIndex = newFitsImage->getKeyNames().indexOf(QString("DATE-OBS"));
        date_obs = newFitsImage->getKeyValues().at(dateIndex);
    }
    else
    {
        date_obs = QString("Date???");
    }

    if (newFitsImage->getKeyNames().contains(QString("TIME")))
    {

        int timeIndex = newFitsImage->getKeyNames().indexOf(QString("TIME"));
        time_obs = newFitsImage->getKeyValues().at(timeIndex);
    }
    else
    {
        time_obs = QString("Time???");
    }

    rMatImage->setDate_obs(date_obs);
    rMatImage->setTime_obs(time_obs);
    rMatImage->setDate_time(date_obs + QString(" ") + time_obs);
    nKeys = newFitsImage->getNKeys();

}

void ImageManager::loadRaw()
{
    newRawImage = new RawImage(filePathQStr);
    //newRawImage = new RawImage2();
    rMatImage = new RMat(newRawImage->getMatCFA(), true, instruments::DSLR);
    rMatImage->setTEMP(newRawImage->getTEMP());
    rMatImage->setXPOSURE(newRawImage->getXPOSURE());
    nKeys = newRawImage->getKeyNames().size();

//    rMatImage->setWbRed(newRawImage->getWbRed());
//    rMatImage->setWbGreen(newRawImage->getWbGreen());
//    rMatImage->setWbBlue(newRawImage->getWbBlue());


}

void ImageManager::loadTiff()
{
    cv::Mat tiffImage = cv::imread(filePathQStr.toStdString(), CV_LOAD_IMAGE_COLOR);
    rMatImage = new RMat(tiffImage, false, instruments::TIFF);

}

void ImageManager::createTableWidget()
{
    // The table widget displays the header keywords/value pairs.

    tableWidget = new QTableWidget(nKeys, 3);
    tableWidget->verticalHeader()->setVisible(false);
    tableWidget->horizontalHeader()->setVisible(true);
    tableWidget->setMinimumWidth(400);
    tableWidget->setMaximumWidth(800);

    tableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);

    tableWidget->setHorizontalHeaderItem(0, new QTableWidgetItem(QString("Keywords")));
    tableWidget->setHorizontalHeaderItem(1, new QTableWidgetItem(QString("Values")));
    tableWidget->setHorizontalHeaderItem(2, new QTableWidgetItem(QString("Comments")));

    tableWidget->horizontalHeaderItem(0)->setTextAlignment(Qt::AlignLeft);
    tableWidget->horizontalHeaderItem(1)->setTextAlignment(Qt::AlignLeft);
    tableWidget->horizontalHeaderItem(2)->setTextAlignment(Qt::AlignLeft);

    tableWidget->setColumnWidth(1, 150);
    tableWidget->setColumnWidth(2, 400);


    if (fitsExtList.contains(fileExt))
    {
        setupFitsTableWidget();
    }
    else if (rawExtList.contains(fileExt))
    {
        setupRawTableWidget();
    }
    else if (fileExt == QString("tiff"))
    {
        setupTiffTableWidget();
    }
    else
    {
        qDebug("Unsupported file extension");
        error = 1;
        return;
    }

    tableWidget->adjustSize();
}

void ImageManager::setupFitsTableWidget()
{
    QFont tFont;
    tFont.setPointSize(12);
    tFont.setFamily("Arial");
    int rowHeight = 20;

    for (int ii=0; ii<newFitsImage->getNKeys(); ii++)
    {
       tableWidget->setItem(ii, 0,  new QTableWidgetItem(newFitsImage->getKeyNames().at(ii)));
       tableWidget->setItem(ii, 1,  new QTableWidgetItem(newFitsImage->getKeyValues().at(ii)));
       tableWidget->setItem(ii, 2,  new QTableWidgetItem(newFitsImage->getKeyComments().at(ii)));
       tableWidget->item(ii, 0)->setFont(tFont);
       tableWidget->item(ii, 1)->setFont(tFont);
       tableWidget->item(ii, 2)->setFont(tFont);
       tableWidget->setRowHeight(ii, rowHeight);
    }
    tableWidget->setMinimumHeight((newFitsImage->getNKeys() + 3) * rowHeight);
}

void ImageManager::setupRawTableWidget()
{
    QFont tFont;
    tFont.setPointSize(12);
    tFont.setFamily("Arial");
    int rowHeight = 20;

    for (int ii=0; ii < newRawImage->getKeyNames().size(); ii++)
    {
       tableWidget->setItem(ii, 0,  new QTableWidgetItem(newRawImage->getKeyNames().at(ii)));
       tableWidget->setItem(ii, 1,  new QTableWidgetItem(newRawImage->getKeyValues().at(ii)));
       tableWidget->setItem(ii, 2,  new QTableWidgetItem(newRawImage->getKeyComments().at(ii)));
       tableWidget->item(ii, 0)->setFont(tFont);
       tableWidget->item(ii, 1)->setFont(tFont);
       tableWidget->item(ii, 2)->setFont(tFont);
       tableWidget->setRowHeight(ii, rowHeight);
    }
    tableWidget->setMinimumHeight((nKeys + 3) * rowHeight);
}

void ImageManager::setupTiffTableWidget()
{

}

void ImageManager::fixUset()
{

    for (int x = 0; x < 2048; x++)
    {
        rMatImage->matImage.at<short>(0, x) = 4095;
//        rMatImage->matImage.at<short>(1, x) = 4095;
//        rMatImage->matImage.at<short>(2, x) = 4095;
//        rMatImage->matImage.at<short>(3, x) = 4095;
//        rMatImage->matImage.at<short>(4, x) = 4095;
//        rMatImage->matImage.at<short>(5, x) = 4095;
    }

}

MyFitsImage* ImageManager::getNewFitsImage() const
{
    return newFitsImage;
}

QString ImageManager::getFileName() const
{
    return fileName;
}

QTableWidget* ImageManager::getTableWidget() const
{
    return tableWidget;
}

RMat *ImageManager::getRMatImage() const
{
    return rMatImage;
}

QString ImageManager::getDate_obs() const
{
    return date_obs;
}

QString ImageManager::getTime_obs() const
{
    return time_obs;
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
