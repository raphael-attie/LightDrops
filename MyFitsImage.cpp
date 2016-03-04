#include "MyFitsImage.h"

//#include <cfitsio/fitsio.h> // Windows
#include <fitsio.h>
#include <cmath>
#include <algorithm>

// openmp
#include <QtDebug>

// opencv
#include <opencv2/core.hpp>
#include <opencv2/core/ocl.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

// ArrayFire
#include <arrayfire.h>


using namespace cv;

MyFitsImage::MyFitsImage(QString filePath) :
hduType(0), naxis1(0), naxis2(0), bayer(false), nPixels(0)
{

	fitsfile *fptr;
    nKeys=0;
    int status = 0, morekeys=0, nfound, anynul;
    long naxes[2], firstPixel, ii;
    double nullval;
    char comment[FLEN_COMMENT], keyString[FLEN_VALUE], card[FLEN_CARD];
    char keyword[FLEN_KEYWORD], keyValue[FLEN_VALUE];
    bscale = 1;
    bzero = 0;

	std::string filePathStr(filePath.toStdString());

	if (fits_open_data(&fptr, filePathStr.c_str(), READONLY, &status))
	{
        std::cout << "MyFitsImage:: Error opening FITS file" << std::endl;
		printerror(status);
	}
		

	if (fits_get_hdu_type(fptr, &hduType, &status))
		printerror(status);
	
    printHDUType(hduType);

	// If ZCMPTYPE does not exist, data are uncompressed -> Read NAXIS
	// else, data are compressed, read ZNAXIS. 
    if (fits_read_key(fptr, TSTRING, "ZCMPTYPE", keyString, NULL, &status))
	{   
        qDebug() << "MyFitsImage:: No ZCMPTYPE keyword, reading UNCOMPRESSED data" ;
		status = 0;
		/* read the NAXIS1 and NAXIS2 keyword to get image size */
		if (fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status))
		{
            qDebug() << "MyFitsImage:: Error reading NAXIS keyword";
			printerror(status);
		}

        if (fits_read_key(fptr, TINT, "BITPIX", &bitpix, NULL, &status))
        {
            qDebug() << "MyFitsImage:: Error reading BITPIX keyword";
            printerror(status);
        }
        qDebug() << "MyFitsImage:: BITPIX =" << bitpix;

	}
	else
	{
        qDebug() << "MyFitsImage:: ZCMPTYPE = " << keyString;
        qDebug() << "MyFitsImage:: Reading COMPRESSED data";
		/* read the NAXIS1 and NAXIS2 keyword to get image size */
		if (fits_read_keys_lng(fptr, "ZNAXIS", 1, 2, naxes, &nfound, &status))
		{
            qDebug() << "MyFitsImage:: Error reading ZNAXIS keyword";
			printerror(status);
		}
        if (fits_read_key(fptr, TINT, "ZBITPIX", &bitpix, NULL, &status))
        {
            qDebug() << "MyFitsImage:: Error reading ZBITPIX keyword";
            printerror(status);
        }
        qDebug() << "MyFitsImage:: ZBITPIX =" << bitpix;

	}

    if (fits_read_key(fptr, TINT, "BSCALE", &bscale, NULL, &status))
    {
        qDebug() << "MyFitsImage:: No BSCALE keyword, setting default BSCALE = 1";
        status = 0;
    }
    qDebug() << "MyFitsImage:: BSCALE =" << bscale;

    if (fits_read_key(fptr, TINT, "BZERO", &bzero, NULL, &status))
    {
        qDebug() << "MyFitsImage:: No BZERO keyword, setting default BZERO = 0";
        status = 0;
    }
    qDebug() << "MyFitsImage:: BZERO =" << bzero;


    if (fits_read_key(fptr, TLOGICAL, "BAYER", &bayer, NULL, &status))
    {
        qDebug() << "MyFitsImage:: No BAYER keyword" ;
        status = 0;
    }
    qDebug() << "MyFitsImage:: BAYER =" << bayer;



	naxis1 = naxes[0];
	naxis2 = naxes[1];
	nPixels = naxis1 * naxis2; // Total number of pixels in the image

    fits_get_hdrspace(fptr, &nKeys, &morekeys, &status);

     for (ii=1; ii<=nKeys; ii++)
     {
         fits_read_keyn(fptr, ii, keyword, keyValue, comment, &status);
         //qDebug() << "keyn Keyword"<<ii<<":" << keyword << keyValue << comment;
         keyNames<< QString(keyword);
         QString tempValue(keyValue);
         tempValue = tempValue.simplified();
         tempValue.replace(" ", "");
         tempValue.replace("'", "");
         keyValues<< tempValue;
         keyComments<< QString(comment);
     }


    // anynul is set to 1 if there any undefined pixel value.
    anynul = 0;
    // undefined (e.g., blank) pixels are set to NaN.
    nullval = NAN;

	firstPixel	= 1;

    if (bitpix == USHORT_IMG)
    {
        ushort* image1D = new ushort[nPixels]();
        if (fits_read_img(fptr, TUSHORT, firstPixel, nPixels, &nullval, image1D, &anynul, &status))
            printerror(status);

        matFits = Mat(naxis2, naxis1, CV_16U, image1D);

        delete[] image1D;
    }
    else if (bitpix == SHORT_IMG && bzero == 0)
    {
        short int* image1D = new short int[nPixels]();
        if (fits_read_img(fptr, TSHORT, firstPixel, nPixels, &nullval, image1D, &anynul, &status))
            printerror(status);

        matFits = Mat(naxis2, naxis1, CV_16S, image1D);

        delete[] image1D;
    }
    else if (bitpix == SHORT_IMG && bzero == 32768)
    {
        ushort* image1D = new ushort[nPixels]();
        if (fits_read_img(fptr, TUSHORT, firstPixel, nPixels, &nullval, image1D, &anynul, &status))
            printerror(status);

        matFits = Mat(naxis2, naxis1, CV_16U, image1D);

        delete[] image1D;
    }
    else if (bitpix == SHORT_IMG  || bitpix == FLOAT_IMG || bitpix == LONG_IMG)
    {
        float* image1D = new float[nPixels]();
        if (fits_read_img(fptr, TFLOAT, firstPixel, nPixels, &nullval, image1D, &anynul, &status))
            printerror(status);
        matFits = Mat(naxis2, naxis1, CV_32F, image1D);

        delete[] image1D;
    }
    else
    {
        if (fits_close_file(fptr, &status))
            printerror(status);
        qDebug("MyFitsImage:: data type not supported, returning.");
        return;
    }

    if (fits_close_file(fptr, &status))
        printerror(status);

    //matFits.convertTo(matFits, CV_32F);
    if (matFits.type() != CV_32F)
    {
        matFits.convertTo(matFits, CV_16U);
    }

}

MyFitsImage::~MyFitsImage()
{

}

// getters

qint32 MyFitsImage::getNaxis1() const
{
	return naxis1;
}


qint32 MyFitsImage::getNaxis2() const
{
	return naxis2;
}

bool MyFitsImage::isBayer() const
{
    return bayer;
}

int MyFitsImage::getBscale() const
{
    return bscale;
}

int MyFitsImage::getBzero() const
{
    return bzero;
}

void MyFitsImage::printHDUType(int hduType)
{
	switch (hduType)
	{
    case IMAGE_HDU: std::cout << "HDUType = Primary Array or IMAGE HDU" << std::endl;
		break;
	case ASCII_TBL: std::cout << "HDUType = ASCII table HDU" << std::endl;
		break;
	case BINARY_TBL: std::cout << "HDUType = Binary table HDU" << std::endl;
		break;
	case ANY_HDU: std::cout << "HDUType = matches any type of HDU" << std::endl;
		break;
	default:
		break;
	}
	return;
}

void MyFitsImage::printerror(int status)
{
	if (status)
	{
		fits_report_error(stderr, status); /* print error report */

		exit(status);    /* terminate the program, returning error status */
	}
	return;
}


int MyFitsImage::getNKeys() const
{
    return nKeys;
}

QVector<QString> MyFitsImage::getKeyNames() const
{
    return keyNames;
}

QVector<QString> MyFitsImage::getKeyValues() const
{
    return keyValues;
}

QVector<QString> MyFitsImage::getKeyComments() const
{
    return keyComments;
}
