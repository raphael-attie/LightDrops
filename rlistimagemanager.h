#ifndef RLISTIMAGEMANAGER_H
#define RLISTIMAGEMANAGER_H

#include <QtCore>

#include "imagemanager.h"
#include "rmat.h"

class RListImageManager
{
public:
    RListImageManager(QList<QUrl> urlList);

    // return buffer of current image.
    void* fetchData(int frameIndex);
    ImageManager* fetchImageManager(int frameIndex);

    //public member (which shall be passed by reference)
    QList<RMat> rMatImageList;

    //getters
    QList<QUrl> getUrlList();
    QList<ImageManager*> getImageManagerList();


private:
    QList<QUrl> urlList;
    QList<ImageManager*> imageManagerList;

};

#endif // RLISTIMAGEMANAGER_H
