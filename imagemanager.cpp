#include "imagemanager.h"

#include <QUrl>
#include <QFileInfo>
#include <QHeaderView>

ImageManager::ImageManager(QString filePathQStr) :
    error(0), newFitsImage(NULL), fitsSeries(NULL), newRawImage(NULL), rMatImage(NULL)
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
        createTableWidget();
    }
    else if (rawExtList.contains(fileExt))
    {
        loadRaw();
    }
    else
    {
        qDebug("Unknown file extension");
        error = 1;
        return;
    }

    rMatImage->setImageTitle(fileName);
    rMatImage->setFileInfo(fileInfo);

    qDebug() << "ImageManager:: dataMin = " << rMatImage->getDataMin();
    qDebug() << "ImageManager:: dataMax = " << rMatImage->getDataMax();

}

ImageManager::~ImageManager()
{
    qDebug("ImageManager::~ImageManager() calling ImageManager destructor");
    delete tableWidget;
    delete rMatImage;
    delete newFitsImage;
    delete newRawImage;
}

bool ImageManager::getError()
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
        qDebug("ImageManager::loadFits():: loading USET data");
        instrument = instruments::USET;
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
}

void ImageManager::loadRaw()
{
    newRawImage = new RawImage(filePathQStr);

    rMatImage = new RMat(newRawImage->matCFA, true, instruments::DSLR);

    rMatImage->setWbRed(newRawImage->getWbRed());
    rMatImage->setWbGreen(newRawImage->getWbGreen());
    rMatImage->setWbBlue(newRawImage->getWbBlue());
}

void ImageManager::createTableWidget()
{
    tableWidget = new QTableWidget(newFitsImage->getNKeys(), 3);
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

    int rowHeight = 20;
    QFont tFont;
    tFont.setPointSize(12);
    tFont.setFamily("Arial");

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
    tableWidget->adjustSize();
}

MyFitsImage* ImageManager::getNewFitsImage()
{
    return newFitsImage;
}

QString ImageManager::getFileName()
{
    return fileName;
}

QTableWidget* ImageManager::getTableWidget() const
{
    return tableWidget;
}

RMat *ImageManager::getRMatImage()
{
    return rMatImage;
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
