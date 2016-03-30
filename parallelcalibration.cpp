#include "parallelcalibration.h"

ParallelCalibration::ParallelCalibration(QList<QUrl> lightFiles, RMat* masterDark, RMat* masterFlat, QList<RMat*> rMatList)
{
    files = lightFiles;
    // Use the RMat copy constructor for dark and flat. This is cheap.
    dark = RMat(*masterDark);
    if (masterFlat == NULL)
    {
        flat = RMat();
    }
    else
    {
        flat = RMat(*masterFlat);
    }

    resultList = rMatList;
}


void ParallelCalibration::operator ()(const cv::Range& range) const
{
    for(int i = range.start; i < range.end; i++)
    {
        QString filePathQStr = files.at(i).toLocalFile();
        ImageManager *lightManager = new ImageManager(filePathQStr);
        cv::Mat matResult;
        cv::Mat lightMat = lightManager->getRMatImage()->matImage;

        if (lightMat.type() != CV_32F)
        {
            lightMat.convertTo(lightMat, CV_32F);
        }
        if (dark.matImage.type() != CV_32F)
        {
            dark.matImage.convertTo(dark.matImage, CV_32F);
        }

        cv::subtract(lightMat, dark.matImage, matResult);

        if (!flat.matImage.empty())
        {
            if (flat.matImage.type() != CV_32F)
            {
                flat.matImage.convertTo(flat.matImage, CV_32F);
            }
            cv::divide(matResult, flat.matImage, matResult);
        }
        /// Copy the result into the list.
        /// This may be an overkill. We could output directly into the elements of the list
        matResult.copyTo(resultList.at(i)->matImage);
        //delete lightManager;
    }
}
