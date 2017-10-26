#include "rlistimagemanager.h"


RListImageManager::RListImageManager()
{

}

RListImageManager::RListImageManager(QList<QUrl> urlList)
{
    this->urlList = urlList;

    for (long ii = 0; ii < urlList.size(); ii++)
    {
        ImageManager* newImageManager = new ImageManager(urlList.at(ii));

        if (!newImageManager->getError())
        {
            imageManagerList << newImageManager;
            rMatImageList << newImageManager->getRMatImage();
            tableWidgetList << newImageManager->getTableWidget();
        }
        else
        {
            qDebug("RListImageManager::RListImageManager:: Could not load file");
            return;
        }
    }
}



RListImageManager::~RListImageManager()
{
    qDebug("RListImageManager::~RListImageManager() calling RListImageManager destructor");

    if (!imageManagerList.empty())
    {
        qDebug("ListImageManager:: Cleaning up imageManagerList");
        qDeleteAll(imageManagerList);
        imageManagerList.clear();
    }

       rMatImageList.clear();

 if (!tableWidgetList.empty())
       {
           qDebug("Cleaning up tableWidgetList");
           tableWidgetList.clear();
       }
}

void RListImageManager::loadData(QList<QUrl> urlList)
{
    this->urlList = urlList;

    cleanLists();

    for (long ii = 0; ii < urlList.size(); ii++)
    {
        QString filePathQStr = urlList.at(ii).toLocalFile();
        qDebug() << "RListImageManager::RListImageManager():: File: " << filePathQStr;

        ImageManager* newImageManager = new ImageManager(urlList.at(ii));

        if (!newImageManager->getError())
        {
            imageManagerList << newImageManager;
            rMatImageList << newImageManager->getRMatImage();
            tableWidgetList << newImageManager->getTableWidget();
        }
        else
        {
            qDebug("RListImageManager::RListImageManager:: Could not load file");
            return;
        }
    }
}


//getters

QList<QUrl> RListImageManager::getUrlList()
{
    return urlList;
}

QList<ImageManager*> RListImageManager::getImageManagerList() const
{
    return imageManagerList;
}

QList<QTableWidget*> RListImageManager::getTableWidgetList() const
{
    return tableWidgetList;
}

QList<RMat *> RListImageManager::getRMatImageList()
{
    return rMatImageList;
}

QVector<double> RListImageManager::getMeanSeries() const
{
    std::cout << "RListImageManager::getMeanSeries()" << std::endl;

    QVector<double> meanSeries;
    for (int i =0; i<rMatImageList.size(); i++)
    {
        meanSeries << (double) rMatImageList.at(i)->getMean();
    }
    return meanSeries;
}

void RListImageManager::cleanLists()
{
    if (!imageManagerList.empty())
    {
        qDebug("Cleaning up imageManagerList");
        qDeleteAll(imageManagerList);
        imageManagerList.clear();
    }

    if (!rMatImageList.empty())
    {
        qDebug("Cleaning up rMatImageList");
        rMatImageList.clear();
    }

    if (!tableWidgetList.empty())
    {
        qDebug("Cleaning up tableWidgetList");
        tableWidgetList.clear();
    }

}

