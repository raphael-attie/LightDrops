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
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>



// this will save performance by disbling any range checking
//#define BOOST_DISABLE_ASSERTS 1

class MyFitsImage

{
public:
	MyFitsImage(QString filePath);
	~MyFitsImage();

    cv::Mat matFits;

	// getters

	qint32 getNaxis1() const;
	qint32 getNaxis2() const;
    bool isBayer() const;
    int getBscale() const;
    int getBzero() const;

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
    int bitpix;
    int bscale;
    int bzero;
    bool bayer;


	// dimensions
	qint32 naxis1; 
	qint32 naxis2;

	long nPixels;


    int nKeys;
    QVector<QString> keyNames;
    QVector<QString> keyValues;
    QVector<QString> keyComments;

	
};

#endif

