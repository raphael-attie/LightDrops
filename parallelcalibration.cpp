#include "parallelcalibration.h"

ParallelCalibration::ParallelCalibration(QString dir, QList<QUrl> lightFiles, RMat masterDark, RMat masterFlat, QList<RMat> &rMatList):
    exportDir(dir), files(lightFiles), dark(masterDark), flat(masterFlat), resultList(rMatList)
{

}

void ParallelCalibration::operator ()(const cv::Range& range) const
{
    for(int i = range.start; i < range.end; i++)
    {
        QString filePathQStr = files.at(i).toLocalFile();
        ImageManager *lightManager = new ImageManager(filePathQStr);
        RMat tempResult;
        cv::subtract(lightManager->rMatImage.matImage, dark.matImage, tempResult.matImage);
        cv::divide(tempResult.matImage, flat.matImage, tempResult.matImage);
        tempResult.matImage.copyTo(resultList.at(i).matImage);
    }
}
