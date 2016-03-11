#include "rlistimagemanager.h"

RListImageManager::RListImageManager(QList<QUrl> urlList)
{
    this->urlList = urlList;

    for (long ii = 0; ii < urlList.size(); ii++)
    {
        QString filePathQStr = urlList.at(ii).toLocalFile();
        qDebug() << "File: " << filePathQStr;
        ImageManager *imageManager = new ImageManager(filePathQStr);
        if (imageManager->rMatImage.empty())
        {
            return;
        }
        imageManagerList << imageManager;
        rMatImageList.push_back(imageManager->rMatImage);
        tableWidgetList << imageManager->getTableWidget();
    }
}


void* RListImageManager::fetchData(int frameIndex)
{
    return imageManagerList.at(frameIndex)->rMatImage.matImage.data;
}

ImageManager* RListImageManager::fetchImageManager(int frameIndex)
{
    return imageManagerList.at(frameIndex);
}

//getters

QList<ImageManager*> RListImageManager::getImageManagerList()
{
    return imageManagerList;
}

QList<QUrl> RListImageManager::getUrlList()
{
    return urlList;
}

QList<QTableWidget*> RListImageManager::getTableWidgetList() const
{
    return tableWidgetList;
}

