#include "RFrame.h"
#include <iostream>

#include <QtGui>
#include <QApplication>
#include <QScreen>
#include <QMouseEvent>
#include <QDebug>
#include <QLabel>
#include <QHBoxLayout>

RFrame::RFrame(QWidget *parent) : QWidget(parent), image(NULL)
{
    initialize();
}

RFrame::RFrame(QImage *image, uint magnification) : QWidget()
{
    // Default displayed size of the widget
    // This default should be greater for retina display
    this->image = image;
    this->magnification = magnification;

    initialize();
}

RFrame::~RFrame()
{

}

void RFrame::initialize()
{
    sx = 1;
    sy = -1;

    naxis1 = 60;
    naxis2 = 60;
    imageCoordX = 0;
    imageCoordY = 0;
    drawCross = false;
    magnification = 2;

    setCursor(Qt::CrossCursor);

    linePenBlack = QPen(Qt::black);
    linePenBlack.setWidthF(0);

    linePenRed = QPen(QColor(240, 20, 20, 255));
    linePenRed.setWidthF(1);

    linePenGreen = QPen(QColor(0, 255, 0, 255));
    linePenGreen.setWidthF(1);


    valueFont = QFont("Times", 12, QFont::Normal);

    // Set the widget size

    pixelDensity = 1;//myScreen->devicePixelRatio();
    resizeFac  = 1 / pixelDensity * magnification;

    // Set the widget size to the same dimension as the painted Image.
    widgetWidth  = naxis1 * resizeFac;
    widgetHeight = naxis2 * resizeFac;

    // Set the scaling factor for the paintEvent
    sx = resizeFac;
    sy = -sx;

    // One need to size the widget here, otherwise it doesn't show at all
    //resize(fitsWidgetWidth , fitsWidgetHeight);
    setFixedSize(widgetWidth , widgetHeight);

    int infoHeight = 15;
    int infoWidth = this->width();
    rectText1 = QRect(1, this->height()- infoHeight, infoWidth, infoHeight);
    rectText2 = QRect(1, 0, this->width(), infoHeight);

    posText1 = QPointF(1, this->height() - 3);
    posText2 = QPointF(1, 12);

    grayTransparentBrush = QBrush(QColor(128, 128, 128, 100));

    setupCrossSquare();
}

void RFrame::paintEvent(QPaintEvent *)
{

    if (image == NULL)
    {
        return;
    }
//    QElapsedTimer timer1;
//    timer1.start();

	QPainter painter(this);
    painter.scale(sx, sy);

    QPointF origin(0, -naxis2);
    painter.drawImage(origin, *image);

    //qDebug()<<"Time for paintEvent (display): " <<  timer1.elapsed() << "ms";


    if (drawCross)
    {
        //setCursor(Qt::BlankCursor);

        // Green
        painter.setPen(linePenGreen);
//        painter.drawLine(hLineLeft);
//        painter.drawLine(vLineBot);
//        painter.drawLine(hLineRight);
//        painter.drawLine(vLineTop);
//        painter.drawRect(subRect1);

        painter.drawLine(hLineCross);
        painter.drawLine(vLineCross);

        // Black
        painter.setPen(linePenBlack);
//        painter.drawLine(hLineLeft);
//        painter.drawLine(vLineBot);
//        painter.drawLine(hLineRight);
//        painter.drawLine(vLineTop);
//        painter.drawRect(subRect2);

        painter.drawLine(hLineCross);
        painter.drawLine(vLineCross);


        painter.scale(1/sx, 1/sy);
        // Display intensity and coordinates where pointed by cursor
        painter.setPen(linePenGreen);
        painter.setFont(valueFont);

        painter.fillRect(rectText1, grayTransparentBrush);
        QString text1(QString("I(X,Y) = ")+cursorText);
        painter.drawText(posText1, text1);

        painter.fillRect(rectText2, grayTransparentBrush);
        QString text2 = QString("X: %1 ; Y: %2").arg(imageCoordX).arg(imageCoordY);
        painter.drawText(posText2, text2);

    }

}

void RFrame::mouseMoveEvent(QMouseEvent *event)
{
    event->ignore();
}

void RFrame::mousePressEvent(QMouseEvent *event)
{
    event->ignore();
}

void RFrame::mouseReleaseEvent(QMouseEvent *event)
{
    setCursor(Qt::CrossCursor);
    event->ignore();
}

void RFrame::setImage(QImage *image)
{
	this->image = image;
}

qint32 RFrame::getNaxis1()
{
	return naxis1;
}

qint32 RFrame::getNaxis2()
{
	return naxis2;
}

quint32 RFrame::getWidgetWidth() const
{
    return widgetWidth;
}

quint32 RFrame::getWidgetHeight() const
{
    return widgetHeight;
}

qreal RFrame::getResizeFac() const
{
    return resizeFac;
}

void RFrame::setImageCoordX(const int coordX)
{
    imageCoordX = coordX;
}

void RFrame::setImageCoordY(const int coordY)
{
    imageCoordY = coordY;
}

void RFrame::setPaintRectangle(bool enabled)
{
    paintRectangle = enabled;
}


void RFrame::setDrawCross(bool enabled)
{
    drawCross = enabled;
}

void RFrame::setupCrossSquare()
{
    // Setup the square-cross dimensions
    qreal squareSize = 4;
    qreal lineLength = 0.5;

    qreal vX = naxis1/2 ;
    qreal vY = -naxis2/2 - 0.5 + 1.5;

    qreal y1 = vY - squareSize/2;
    qreal y2 = y1 + lineLength;
    vLineBot = QLineF(vX, y1, vX, y2);
    vLineTop = QLineF(vX, y1 + squareSize - lineLength, vX, y1 + squareSize);

    qreal x1 = vX - squareSize/2 ;
    qreal x2 = x1 + lineLength;
    hLineLeft = QLineF(x1, vY, x2, vY);
    hLineRight = QLineF(x1 + squareSize -lineLength, vY, x1 + squareSize, vY);

    subRect1 = QRectF(x1, y1, squareSize, squareSize);
    subRect2 = QRectF(x1, y1, squareSize, squareSize);

    qreal hCrossX1 = naxis1/2 - 2;
    qreal hCrossX2 = hCrossX1 + 4;
    qreal hCrossY = - naxis2/2 + 1;
    hLineCross = QLineF(hCrossX1, hCrossY, hCrossX2, hCrossY);

    qreal vCrossX = naxis1/2;
    qreal vCrossY1 = -naxis2/2 - 1;
    qreal vCrossY2 = vCrossY1 + 4;
    vLineCross = QLineF(vCrossX, vCrossY1, vCrossX, vCrossY2);

}

void RFrame::setCursorPoint(QPointF cursorPoint)
{
    this->cursorPoint = cursorPoint;
}

void RFrame::setCursorText(QString cursorText)
{
    this->cursorText = cursorText;
}
