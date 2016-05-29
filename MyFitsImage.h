#ifndef MYFITSIMAGE_H
#define MYFITSIMAGE_H

#include "winsockwrapper.h"
#include <QDataStream>
#include <QtCore>

#include <cmath>

#include <iostream>
#include <vector>
#include <cstring>

// opencv
#include <opencv2/world.hpp>



// this will save performance by disbling any range checking
//#define BOOST_DISABLE_ASSERTS 1


class MyFitsImage

{
public:
	MyFitsImage(QString filePath);
	~MyFitsImage();

	// getters

	qint32 getNaxis1() const;
	qint32 getNaxis2() const;
    bool isBayer() const;
    float getBscale() const;
    float getExpTime() const;
    int getBzero() const;
    cv::Mat getMatFits() const;

    int getNKeys() const;
    QVector<QString> getKeyNames() const;
    QVector<QString> getKeyValues() const;
    QVector<QString> getKeyComments() const;

    // static
	static void printerror(int status);
	static void printHDUType(int hduType);

private:

	//valarray<float> *fitsValArray;
	//array_2D_double *myImage;
	int hduType;
    // dimensions
    int naxis1;
    int naxis2;
    int nPixels;
    int nKeys;
    float bscale, expTime;
    int bzero;
    bool bayer;
    int bitpix;

    ushort* image1D_ushort;
    float* image1D_float;
    short int* image1D_shortint;


    cv::Mat matFits;

    QVector<QString> keyNames;
    QVector<QString> keyValues;
    QVector<QString> keyComments;

	
};

#endif

