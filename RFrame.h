#ifndef RFrame_H
#define RFrame_H

#include "winsockwrapper.h"
#include <QWidget>
#include <QScrollArea>
#include <qpainter.h>
#include <QScreen>
#include <QMouseEvent>

class RFrame : public QWidget // QWidget used as the QPaintDevice.
{
	Q_OBJECT

public:

    RFrame(QWidget *parent =0);
    RFrame(QImage *image, uint magnification=1);
    ~RFrame();

    void initialize();
	void setImage(QImage *image);
	qint32 getNaxis1();
	qint32 getNaxis2();
    quint32 getWidgetWidth() const;
    quint32 getWidgetHeight() const;
    qreal getResizeFac() const;
    void setImageCoordX(int coordX);
    void setImageCoordY(int coordY);
    void setPaintRectangle(bool enabled);
    void setDrawCross(bool enabled);
    void setCursorPoint(QPointF cursorPoint);
    void setCursorText(QString cursorText);
    void setupCrossSquare();

protected:

    void paintEvent(QPaintEvent *);
    void mouseMoveEvent(QMouseEvent *event);
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

private:



	QImage *image;
	qreal sx;
	qreal sy;
    qreal resizeFac;
	qint32 naxis1;
	qint32 naxis2;
    qreal pixelDensity;

    long imageCoordX;
    long imageCoordY;
    bool paintRectangle;
    bool drawCross;

    QLineF vLineBot;
    QLineF hLineLeft;
    QLineF vLineTop;
    QLineF hLineRight;
    QRectF subRect1;
    QRectF subRect2;
    QPointF cursorPoint;
    QString cursorText;

    QPen linePenBlack;
    QPen linePenRed;
    QPen linePenGreen;
    QFont valueFont;

    quint32 widgetWidth;
    quint32 widgetHeight;
    uint magnification;

    QRect rectText1;
    QRect rectText2;
    QPointF posText1;
    QPointF posText2;

    QBrush grayTransparentBrush;

};

#endif
