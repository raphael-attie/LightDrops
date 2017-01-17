#include "RawImage.h"

// libraw
#include <libraw/libraw.h>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/core/ocl.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

// Qt
#include <iostream>
#include <QObject>
#include <QString>
#include <QDebug>
#include <QColor>
#include <QImage>
#include <QElapsedTimer>
#include <QDateTime>

using namespace cv;
using namespace std;

RawImage::RawImage()
{

}

RawImage::RawImage(QString filePath):filePath_(filePath), imageProcessed(NULL)
{
    QElapsedTimer timer1, timer2;

    int error = 0;

//    qDebug() << "Raw file path:" + filePath.toLatin1();

//    unsigned char x;
//    ifstream inputF(filePathStr, ios::binary);
//    inputF.seekg(0, ios::beg);
//    //qDebug() << "inputF.tellg() = " << inputF.tellg();
//    inputF >> noskipws;
//    istream_iterator<unsigned char> it(inputF);

//    unsigned char *array = new unsigned char[100];
//    int i = 0;
//    while(inputF.tellg() < 60)
//    {
//        qDebug() << "inputF.tellg() = " << inputF.tellg();
//        inputF >> hex >> setw(2) >> array[i];
//        i++;
//    }

//    for( int i = 0; i < 60; i++ )
//    {
//          cout << uppercase << showbase << hex << array[i] << endl;
//    }

//    int position = 0;
//    char * array = new char[100];

//    while(inputF.tellg() < 60)
//    {
//        qDebug() << "inputF.tellg() = " << inputF.tellg();
//        inputF.get(array[position]);
//        position++;
//    }
//    qDebug() << "position = " << position;
//    array[position-1] = '\0';
//    ///displaying the data stored in the array
//    qDebug() << "array[0] = " << array[0];
//    qDebug() << "array[1] = " << array[1];
//    qDebug() << "array[2] = " << array[2];




    /// Sample program to print the Exif metadata of an image



//    for (Exiv2::ExifData::const_iterator i = exifData.begin(); i != exifData.end() ; i++)
//    {
//        qDebug() << "i->key() = " << QString::fromStdString((*i).key());
//    }

//    for (int i = 0 ; i < strList.size() ; i++)
//    {

//        Exiv2::Exifdatum data = exifData[strList[i]];

//        QString dataName = QString::fromStdString(strList[i]);
//        QString dataValue = QString::fromStdString(data.toString());

//        qDebug() << dataName << " = " << dataValue;
//    }





//   Exiv2::ExifData::const_iterator end = exifData.end();
//   for (Exiv2::ExifData::const_iterator i = exifData.begin(); i != end; ++i)
//   {
//       QString keyName = QString::fromStdString(i->key());
//       QString keyValue = QString::fromStdString(i->value().toString());
//       qDebug() << keyName << " : ";
//    }


    error = rawProcess.open_file(filePath_.toStdString().c_str());

    if (error !=0)
    {
        qDebug() << "error on open_file = " << error;
    }

    //rawProcess.imgdata.params.use_camera_matrix = 0;
    rawProcess.imgdata.params.use_camera_wb = 0;
    rawProcess.imgdata.params.no_auto_bright = 1;
    //rawProcess.imgdata.params.user_qual = -1;

    extractExif();

    timer1.start();

    error = rawProcess.unpack();
    if (error!=0)
    {
        qDebug() << "unpack() error =" << error;
    }
    qDebug()<<"Time for rawProcess.unpack()=" <<  timer1.elapsed() << "ms";
    qDebug() << "Bayer Pattern: rawProcess.imgdata.idata.cdesc=" << rawProcess.imgdata.idata.cdesc;

    naxis1 = rawProcess.imgdata.sizes.iwidth;
    naxis2 = rawProcess.imgdata.sizes.iheight;
    nPixels = naxis1 * naxis2;

    //    rawProcess.subtract_black();

    // Advised by lexa (Author of Libraw):
    //    Each row contains non-visible pixels on both ends:
    //    - on left (0.... left_margin-1)
    //    - and on right (from visible area end to row end)
    //    Also, rows may be aligned for efficient SSE/AVX access. LibRaw internal code do not align rows, but if you use LibRaw+RawSpeed, RawSpeed will align rows on 16 bytes.
    //    So, it is better to use imgdata.sizes.raw_pitch (it is in bytes, so divide /2 for bayer data) instead of raw_width.

    //    int raw_width = (int) rawProcess.imgdata.sizes.raw_width;
    int raw_width = (int) rawProcess.imgdata.sizes.raw_pitch/2;
    int top_margin = (int) rawProcess.imgdata.sizes.top_margin;
    int left_margin = (int) rawProcess.imgdata.sizes.left_margin;
    int first_visible_pixel = (int) (raw_width * top_margin + left_margin);

    qDebug() << "first_visible_pixel =" << first_visible_pixel;


    matCFA.create(naxis2, naxis1, CV_16UC1);
    ushort *rawP = matCFA.ptr<ushort>(0);

    // The following block re-assign the data in a more contiguous form,
    // where the "dark" margins (top and left) have been removed from the raw buffer (rawdata.raw_image).
    // Hence the dimension of this new buffer (rawP) match the final image size of the camera.
    // We convert them afterwards to float.
    // These steps cost about 150 ms at most on a 22 MP DSLR, core i5 macbook Pro Retina 2013.
    // Without parallelization (QThread or alike), this is about 5-10% of the total time to load and display the image.
    long j=0;
    for(int r = 0; r < naxis2; r++)
    {
        for(int c=0; c < naxis1; c++)
        {
            rawP[j] = rawProcess.imgdata.rawdata.raw_image[(r + top_margin)*raw_width + left_margin + c];
            j++;
        }
    }
   //matCFA.convertTo(matCFA, CV_32F);
    matCFA.convertTo(matCFA, CV_16U);


    // White balance
    wbRed = rawProcess.imgdata.color.cam_mul[0] /1000.0f;
    wbGreen = rawProcess.imgdata.color.cam_mul[1]/1000.0f;
    wbBlue = rawProcess.imgdata.color.cam_mul[2]/1000.0f;

    qDebug("white balance cam_mul =%f , %f , %f" , wbRed, wbGreen, wbBlue);

    //timer2.restart();
    //std::vector<Mat> rgbChannels(3);
    //split(bayerMat32, rgbChannels);
    //imRed = rgbChannels.at(0);
    //imGreen = rgbChannels.at(1);
    //imBlue = rgbChannels.at(2);

    //imRed = imRed * redBalance;
    //imGreen = imGreen * greenBalance;
    //imBlue = imBlue * blueBalance;

    //std::vector<cv::Mat> matVector;
    //matVector.push_back(imRed);
    //matVector.push_back(imGreen);
    //matVector.push_back(imBlue);
    //cv::merge(matVector, bayerMat32);


    //qDebug()<<"Time for openCV white balancing:" <<  timer2.elapsed() << "ms";


    // raw2image()
    //    timer2.start();
    //    error = rawProcess.raw2image();
    //    if (error !=0)
    //    {
    //        qDebug() << "error on raw2image()";
    //    }
    //    qDebug()<<"Time for rawProcess.raw2image()" <<  timer2.elapsed() << "ms";


    //    for(int i = 0;i < nPixels; i++)
    //    {
    //        qDebug("i=%d R=%d G=%d B=%d G2=%d\n", i,
    //                     rawProcess.imgdata.image[i][0],
    //                     rawProcess.imgdata.image[i][1],
    //                     rawProcess.imgdata.image[i][2],
    //                     rawProcess.imgdata.image[i][3]
    //             );
    //    }

    //    timer2.restart();
    //    error = rawProcess.dcraw_process();
    //    if (error !=0)
    //    {
    //        qDebug() << "error on dcraw_process()";
    //    }
    ////    qDebug()<<"Time for rawProcess.dcraw_process()" <<  timer2.elapsed() << "ms";

    //    ushort *redP = imRed.ptr<ushort>(0);
    //    ushort *greenP = imGreen.ptr<ushort>(0);
    //    ushort *blueP = imBlue.ptr<ushort>(0);

    //    for (int i = 0; i < nPixels; i++)
    //    {

    //        redP[i] = rawProcess.imgdata.image[i][0];
    //        greenP[i] = rawProcess.imgdata.image[i][1];
    //        blueP[i] = rawProcess.imgdata.image[i][2];
    //    }


    //    std::vector<cv::Mat> matVector;
    //    matVector.push_back(imRed);
    //    matVector.push_back(imGreen);
    //    matVector.push_back(imBlue);
    ////    //matVector.push_back(imGreen2);

    //    cv::Mat bayerMat32b;
    //    cv::merge(matVector, bayerMat32b);
    //    bayerMat32b.convertTo(bayerMat32b, CV_32FC3);

    //    qDebug()<<"Time for dcraw_process() and channel assignment: " <<  timer2.elapsed() << "ms";

    //    Mat bayerMat16, rgbMat16;
    //    bayerMat32.convertTo(bayerMat16, CV_16UC1);
    //    cvtColor(bayerMat16, rgbMat16, CV_BayerBG2RGB);
    //    rgbMat16.convertTo(bayerMat32, CV_32FC3);


}

RawImage::~RawImage()
{
    imRed.release();
    imGreen.release();
    imBlue.release();
    delete[] imageProcessed;
}


void RawImage::extractExif()
{
    /// Extract exif data
//    image = Exiv2::ImageFactory::open(filePath_.toStdString());
//    assert(image.get() != 0);
//    image->readMetadata();
//    Exiv2::ExifData &exifData = image->exifData();
//    if (exifData.empty())
//    {
//        std::string error(filePath_.toStdString());
//        error += ": No Exif data found in the file";
//        throw Exiv2::Error(1, error);
//    }
//    Exiv2::Exifdatum temperatureData = exifData["Exif.CanonSi.CameraTemperature"];

    QDateTime dt = QDateTime::fromTime_t( rawProcess.imgdata.other.timestamp );

    dispatchMetaDatum("Brand", QString::fromUtf8(rawProcess.imgdata.idata.make));
    dispatchMetaDatum("Model", QString::fromUtf8(rawProcess.imgdata.idata.model));
    dispatchMetaDatum("Colors", QString::number(rawProcess.imgdata.idata.colors), "Number of color channels");
    dispatchMetaDatum("Bayer pattern", QString::fromUtf8(rawProcess.imgdata.idata.cdesc), "Color pattern of the Bayer matrix");
    dispatchMetaDatum("DATE", dt.toString("MMMM d yyyy hh:mm:ss t"));
    dispatchMetaDatum("NAXIS1", QString::number(rawProcess.imgdata.sizes.width), "Image width (px)");
    dispatchMetaDatum("NAXIS2", QString::number(rawProcess.imgdata.sizes.height), "Image height (px)");
    dispatchMetaDatum("Flip", QString::number(rawProcess.imgdata.sizes.flip), "0: 0; 3: 180 deg; 5: 90 deg CCW; 6: 90 deg CW");
    dispatchMetaDatum("Order", QString::number(rawProcess.imgdata.other.shot_order), "Shot ordered number");
    dispatchMetaDatum("ISO", QString::number(rawProcess.imgdata.other.iso_speed));
    dispatchMetaDatum("XPOSURE", QString::number(rawProcess.imgdata.other.shutter), "Exposure time (s)");
    //dispatchMetaDatum("Temperature", QString::fromStdString(temperatureData.print()), "Temperature of camera sensor in degrees Celsius");

}

LibRaw RawImage::getRawProcess() const
{
    return rawProcess;
}

Mat RawImage::getImRed() const
{
    return imRed;
}

Mat RawImage::getImGreen() const
{
    return imGreen;
}

Mat RawImage::getImBlue() const
{
    return imBlue;
}


qint32 RawImage::getNaxis1() const
{
    return naxis1;
}

qint32 RawImage::getNaxis2() const
{
    return naxis2;
}

int RawImage::getNPixels() const
{
    return nPixels;
}


float RawImage::getDataMaxRed() const
{
    return dataMaxRed;
}

float RawImage::getDataMaxGreen() const
{
    return dataMaxGreen;
}

float RawImage::getDataMaxBlue() const
{
    return dataMaxBlue;
}

float RawImage::getDataMinRed() const
{
    return dataMinRed;
}

float RawImage::getDataMinGreen() const
{
    return dataMinGreen;
}

float RawImage::getDataMinBlue() const
{
    return dataMinBlue;
}

float RawImage::getWbRed() const
{
    return wbRed;
}

float RawImage::getWbGreen() const
{
    return wbGreen;
}

float RawImage::getWbBlue() const
{
    return wbBlue;
}

QVector<QString> RawImage::getKeyNames() const
{
    return keyNames;
}

QVector<QString> RawImage::getKeyValues() const
{
    return keyValues;
}

QVector<QString> RawImage::getKeyComments() const
{
    return keyComments;
}

void RawImage::dispatchMetaDatum(const char* keyName, QString keyValue, const char* comment)
{
    keyNames << QObject::tr(keyName);
    keyValues << keyValue;
    keyComments << QObject::tr(comment);
}
