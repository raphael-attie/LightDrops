#ifndef RopenGLWidget_H
#define RopenGLWidget_H

#include "winsockwrapper.h"
#include <QtCore>
#include <QWidget>
#include <QOpenGLWidget>
#include <QtGui/QWindow>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLTexture>
#include <QMdiSubWindow>

//opencv
//#include <opencv2/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/world.hpp>

#include "rlistimagemanager.h"
//#include "RawImage2.h"
#include "rmat.h"
#include "rsubwindow.h"

//QCustomPlot
#include <qcustomplot/qcustomplot.h>


class QPainter;
class QOpenGLContext;
class QOpenGLPaintDevice;

class ROpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    ROpenGLWidget(QWidget *parent = 0);
    ROpenGLWidget(RListImageManager *rListImageManager, QWidget *parent = 0);
    ROpenGLWidget(QList<RMat*> rMatImageList, QWidget *parent = 0);
    ROpenGLWidget(RMat *rMatImage, QWidget *parent = 0);
    ~ROpenGLWidget();

    void initialize();
    void initSubQImage();
    void updateSubQImage();
    void setupHistoPlots();
    void updateCustomPlotLineItems();

    QList<cv::Mat> matImageListRGB;

    //events
    void focusInEvent(QFocusEvent *event);
    void focusOutEvent(QFocusEvent *event);
    void wheelEvent(QWheelEvent *wheelEvent);

    //getters
    quint32 getImageCoordX();
    quint32 getImageCoordY();
    RListImageManager *getRListImageManager();
    QList<RMat*> getRMatImageList();
    float getResizeFac();

    int getFrameIndex();
    QString getCalibrationType();
    QList<QString> getWindowTitleList();
    QString getWindowTitle();
    RSubWindow *getTableRSubWindow();
    QList<QCustomPlot *> getCustomPlotList();
    QCustomPlot* fetchCurrentCustomPlot();

    float getGamma();
    float getNewMax();
    float getNewMin();
    float getLimbNewMax();
    float getLimbNewMin();
    float getLimbGamma();
    float getAlpha();
    float getBeta();

    float getWbRed();
    float getWbGreen();
    float getWbBlue();

    float getRadius();
    float getScaleLimb();

    QRect getFovRect() const;

    //setters
    void setRMatImageList(QList<RMat*> rMatImageList);
    void setNewMax(float newMax);
    void setNewMin(float newMin);
    void setLimbNewMax(float limbNewMax);
    void setLimbNewMin(float limbNewMin);
    void setLimbGamma(float limbGamma);
    void setAlpha(float newAlpha);
    void setBeta(float newBeta);
    void setAlphaLimb(float alphaLimb);
    void setBetaLimb(float betaLimb);
    void setGamma(float newGamma);
    void setImax(float iMax);
    void setLambda(float lambda);
    void setMu(float mu);
    void setApplyToneMapping(bool state);
    void setUseInverseGausssian(bool state);
    void setScaleLimb(bool state);
    void setImageCoordX(quint32 x);
    void setImageCoordY(quint32 y);

    void setFrameIndex(int newFrameIndex);
    void setCalibrationType(QString calibrationType);

    void setWindowTitleList(QList<QString> newWindowTitleList);
    void setWindowTitle(QString newWindowTitle);
    void setResizeFac(float resizeFac);

    void setWbRed(float wbRed);
    void setWbGreen(float wbGreen);
    void setWbBlue(float wbBlue);

    void setTableSize(QSize size);
    void setRadius(float radius);
    void setDisplayOnlyFirst(bool status);

    //ROI
    void setRoiSelected(bool isSelected);
    void setCircleSelected(bool isSelected);
    void setUseMultiROI(bool isSelected);
    void clearROIs();



signals:

    void mousePressed();
    void gotSelected(ROpenGLWidget* ROpenGLWidget);
    void sendSubQImage(QImage *image, float intensity, int x, int y);
    void sendTableWidget(QTableWidget* tableWidget);
    void sendNewTitle(QString title);
    void roiSignal(QRect rect);
    void roiSignal(QList<QRect> rectList);
    void circleSignal(int circleX, int circleY, int circleRadius);

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

public slots:
    void setupTableWidget(int value);
    void changeFrame(int frameIndex);


private:
    void prepImage();

    bool prepareShaderProgram( const QString vertexShaderPath,
                               const QString fragmentShaderPath );

    QList<RMat*> rMatImageList;

    void loadGLTexture();
    void calculateDefaultSize();


    int nFrames;
    uint counter;
    double dataMax;
    double dataMin;
    float dataRange;
    float intensity;
    int sliderValueHigh;
    int sliderValueLow;
    int sliderValueGamma;

    float newMin, newMax;
    float limbNewMin, limbNewMax, limbGamma;
    float alpha, beta, gamma;
    /// Limb
    float alphaLimb, betaLimb;
    /// ToneMapping
    float iMax, lambda, mu, radius;
    bool applyToneMapping;
    bool useInverseGaussian;
    bool scaleLimb;


    float wbRed;
    float wbGreen;
    float wbBlue;

    float resizeFac;


    RListImageManager *rListImageManager;
    cv::Mat matImageRGB;
    cv::Mat matFrame;
    cv::Mat subMatImage;
    cv::Rect FOV;
    QImage *subQImage;

//    GLuint shaderProgram;
//    GLuint vertexShader;
//    GLuint fragmentShader;
//    GLuint vao;
//    GLuint vbo;
//    GLuint ebo;
//    GLuint texture;

    bool displayOnlyFirst;
    QOpenGLShaderProgram m_shader;
    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vertexBuffer;
    QOpenGLBuffer m_indexBuffer;

    QOpenGLTexture *texture;
    QVector<QOpenGLTexture*> textureVector;

    QList<QString> windowTitleList;
    QString windowTitle;
    QString calibrationType;
    uint frameIndex;

    qint32 naxis1;
    qint32 naxis2;
    uint subNaxis;



    RSubWindow *tableRSubWindow;
    QTableWidget *tableWidget;
    QSize tableSize;

    QList<QCustomPlot*> customPlotList;
    QList<QCPItemLine*> vertLineHighList;
    QList<QCPItemLine*> vertLineLowList;
    QList<QCPItemText*> itemTextHighList;
    QList<QCPItemText*> itemTextLowList;

    QCustomPlot* customPlot;
    QCPItemLine *vertLineHigh;
    QCPItemLine *vertLineLow;
    QCPItemText *itemTextHigh;
    QCPItemText *itemTextLow;

    // ROI QPainter-related
    bool roiSelected;
    bool useMultiROI;
    bool circleSelected;
    QString painterString;

    QPoint initCursorPos;
    int initCursorX;
    int initCursorY;
    int currentCursorX;
    int currentCursorY;

    QList<QRect> painterRectList;
    QRect painterRect;
    QRect circleRect;



    // Field of view in the image reference frame corresponding to the selection rectangle in QOpenGLWidget
    QList<QRect> fovRectList;
    QRect fovRect;
    QRect fovCircleRect;
    int circleX;
    int circleY;
    quint32 circleDiameter;
    quint32 circleRadius;

    quint32 imageCoordX;
    quint32 imageCoordY;
    quint32 imageCircleX;
    quint32 imageCircleY;
    quint32 imageCircleDiameter;
    quint32 imageCircleRadius;


};

#endif // RopenGLWidget_H
