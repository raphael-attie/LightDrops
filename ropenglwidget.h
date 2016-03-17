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

//opencv
#include <opencv2/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include "rlistimagemanager.h"
#include "RawImage.h"
#include "rmat.h"

class QPainter;
class QOpenGLContext;
class QOpenGLPaintDevice;

class ROpenGLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
    Q_OBJECT
public:
    ROpenGLWidget(QWidget *parent = 0);
    ROpenGLWidget(RListImageManager *rListImageManager, QWidget *parent = 0);
    ROpenGLWidget(QList<RMat> & rMatImageList, QWidget *parent = 0);
    ROpenGLWidget(RMat &rImage, QWidget *parent = 0);
    ~ROpenGLWidget();

    void initSubQImage();
    void updateSubQImage();
    void updateInfo();

    QList<RMat> rMatImageList;
    QList<cv::Mat> matImageListRGB;

    //events
    void focusInEvent(QFocusEvent *event);
    void focusOutEvent(QFocusEvent *event);
    void wheelEvent(QWheelEvent *wheelEvent);

    //getters
    QList<RMat> getRMatImageList();
    QSize getOglDefaultSize();
    quint32 getImageCoordX();
    quint32 getImageCoordY();
    RListImageManager *getRListImageManager();
    int getFrameIndex();
    QString getCalibrationType();
    QList<QString> getWindowTitleList();
    QString getWindowTitle();

    float getGamma();
    float getNewMax();
    float getNewMin();

    float getWbRed();
    float getWbGreen();
    float getWbBlue();


    //setters
    void setNewMax(float newMax);
    void setNewMin(float newMin);
    void setAlpha(float newAlpha);
    void setBeta(float newBeta);
    void setGamma(float newGamma);

    void setFrameIndex(int newFrameIndex);
    void setRListImageManager(RListImageManager *rListImageManager);
    void setCalibrationType(QString calibrationType);


    void setWindowTitleList(QList<QString> newWindowTitleList);
    void setWindowTitle(QString newWindowTitle);

    void setWbRed(float wbRed);
    void setWbGreen(float wbGreen);
    void setWbBlue(float wbBlue);

signals:

    void mousePressed();
    void gotSelected(ROpenGLWidget* ROpenGLWidget);
    void sendSubQImage(QImage *image, float intensity, int x, int y);

protected:
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);

private:
    void initialize();
    cv::Mat prepImage(RMat rMatImage);

    bool prepareShaderProgram( const QString vertexShaderPath,
                               const QString fragmentShaderPath );

    bool isBayer;
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
    float alpha, beta, gamma;

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

    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vertexBuffer;
    QOpenGLBuffer m_indexBuffer;
    QOpenGLShaderProgram m_shader;
    QOpenGLTexture *texture;
    QVector<QOpenGLTexture*> textureVector;

    QList<QString> windowTitleList;
    QString windowTitle;
    QString calibrationType;
    uint frameIndex;

    QSize oglDefaultSize;
    qint32 naxis1;
    qint32 naxis2;
    uint subNaxis;

    quint32 imageCoordX;
    quint32 imageCoordY;

};

#endif // RopenGLWidget_H
