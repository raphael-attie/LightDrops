#include "parallelcalibration.h"

ParallelCalibration::ParallelCalibration(QString dir, QList<QUrl> lightFiles, RMat* masterDark, RMat* masterFlat, QList<RMat*> rMatList):
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
        cv::subtract(lightManager->getRMatImage()->getMatImage(), dark->getMatImage(), tempResult.getMatImage());
        cv::divide(tempResult.getMatImage(), flat->getMatImage(), tempResult.getMatImage());
        tempResult.getMatImage().copyTo(resultList.at(i)->getMatImage());
    }
}
