#include "ropenglwidget.h"

#include <QApplication>
#include <QOpenGLWidget>
#include <QtGui/QWindow>
#include <QtGui/QOpenGLFunctions>
#include <QtGui/QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLTexture>
#include <QScreen>

#include <QMouseEvent>

#include <QElapsedTimer>
#include <QTimer>

//opencv
#include <opencv2/core.hpp>
#include <opencv2/core/ocl.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

#include <unistd.h>

using namespace cv;

// Shader sources
const GLchar* vertexSource =
    "#version 330 core\n"
    "layout (location = 0) in vec3 position;"
    "layout (location = 1) in vec2 texCoord;"
    "out vec2 TexCoord;"
    "void main()"
    "{"
    "    gl_Position = vec4(position, 1.0f);"
    "    TexCoord = texCoord;"
    "}";
const GLchar* fragmentSource =
    "#version 330 core\n"
    "in vec2 TexCoord;"
    "out vec4 color;"
    "uniform usampler2D ourTexture;"
    "void main()"
    "{"
    "    vec3 textureColor = texture(ourTexture, TexCoord).rgb /65535.0 * 20.0;"
    "    color = vec4(textureColor.rgb, 1.0);"
    "}";


ROpenGLWidget::ROpenGLWidget(QWidget *parent) :
    QOpenGLWidget(parent)
{

}

ROpenGLWidget::ROpenGLWidget(RListImageManager *rListImageManager, QWidget *parent)
    :QOpenGLWidget(parent), tableWidget(NULL), tableRSubWindow(NULL)
    , m_shader(0)
    , m_vertexBuffer(QOpenGLBuffer::VertexBuffer)
    , m_indexBuffer(QOpenGLBuffer::IndexBuffer)
{
    this->rListImageManager = rListImageManager;
    this->rMatImageList = rListImageManager->getRMatImageList();
    nFrames = rMatImageList.size();

    initialize();
    //setupTableWidget();
    tableRSubWindow = new RSubWindow();
    tableSize = rListImageManager->getTableWidgetList().at(0)->size();

}

ROpenGLWidget::ROpenGLWidget(QList<RMat *> rMatImageList, QWidget *parent)
    :QOpenGLWidget(parent), rListImageManager(NULL), tableWidget(NULL), tableRSubWindow(NULL)
    , m_shader(0)
    , m_vertexBuffer(QOpenGLBuffer::VertexBuffer)
    , m_indexBuffer(QOpenGLBuffer::IndexBuffer)
{
    this->rMatImageList = rMatImageList;
    nFrames = rMatImageList.size();

    initialize();
}

ROpenGLWidget::ROpenGLWidget(RMat *rMatImage, QWidget *parent)
    :QOpenGLWidget(parent), rListImageManager(NULL), tableWidget(NULL), tableRSubWindow(NULL)
    , m_shader(0)
    , m_vertexBuffer(QOpenGLBuffer::VertexBuffer)
    , m_indexBuffer(QOpenGLBuffer::IndexBuffer)
{
    this->rMatImageList =  QList<RMat*>() << rMatImage;
    //this->rMatImageList.push_back(rMatImage);
    nFrames = 1;

    initialize();
}

ROpenGLWidget::~ROpenGLWidget()
{
    qDebug("ROpenGLWidget::~ROpenGLWidget() destructor");
    m_indexBuffer.destroy();
    m_vertexBuffer.destroy();
    m_vao.destroy();

    delete rListImageManager;

//    glDeleteTextures(1, &texture);

//    glDeleteProgram(shaderProgram);
//    glDeleteShader(fragmentShader);
//    glDeleteShader(vertexShader);

//    glDeleteBuffers(1, &ebo);
//    glDeleteBuffers(1, &vbo);
//    glDeleteVertexArrays(1, &vao);

}


void ROpenGLWidget::initialize()
{
    setFocusPolicy(Qt::StrongFocus);
    frameIndex = 0;
    counter = 0;
    dataMax = rMatImageList.at(0)->getDataMax();
    dataMin = rMatImageList.at(0)->getDataMin();
    dataRange = (float) (dataMax - dataMin);
    calibrationType = QString("");
    naxis1 = rMatImageList.at(0)->getMatImage().cols;
    naxis2 = rMatImageList.at(0)->getMatImage().rows;

    double min;
    double max;

    prepImage();

    cv::minMaxLoc(rMatImageList.at(0)->getMatImage(), &min, &max);
    qDebug("ROpenGLWidget::initialize()::rMatImageList.at(0)->getMatImage()   min = %f , max = %f", min, max );


    matImageRGB = matImageListRGB.at(0);


    intensity = 0;
    newMin = dataMin;
    newMax = dataMax;
    alpha = 1.0f/dataRange;
    beta = (float) (-dataMin / dataRange);
    gamma = 1.0;

    // White balance default
    wbRed = 1.0f;
    wbGreen = 1.0f;
    wbBlue = 1.0f;

    imageCoordX = naxis1/2;
    imageCoordY = naxis2/2;

    initSubQImage();

    /// Create series of image titles
    /// If the frames do not come from disk files, need to create one from the default frame titles
    if (rListImageManager == NULL)
    {
        for (int i = 0; i < nFrames; i++)
        {
            windowTitleList << rMatImageList.at(i)->getImageTitle() + QString::number(i+1);
        }
    }
    else
    {
        for (int i = 0; i < nFrames; i++)
        {
            windowTitleList << rMatImageList.at(i)->getImageTitle();
        }
    }

}


void ROpenGLWidget::prepImage()
{
    for (int ii = 0; ii < nFrames; ii++)
    {
        cv::Mat tempMatRGB = cv::Mat::zeros(naxis2, naxis1, rMatImageList.at(ii)->getMatImage().type());

        if (rMatImageList.at(ii)->isBayer())
        {
            cv::Mat tempMat16;
            rMatImageList.at(ii)->getMatImage().convertTo(tempMat16, CV_16U);
            cv::cvtColor(tempMat16, tempMatRGB, CV_BayerBG2RGB);
            //tempMatRGB.convertTo(tempMatRGB, CV_32FC3);
        }
        else if (rMatImageList.at(ii)->getMatImage().channels() == 1)
        {
            cv::Mat tempMat = rMatImageList.at(ii)->getMatImage();
            cv::cvtColor(tempMat, tempMatRGB, CV_GRAY2RGB);
        }
        else
        {
            rMatImageList.at(ii)->getMatImage().copyTo(tempMatRGB);
        }        
        matImageListRGB.append(tempMatRGB);
    }

}

void ROpenGLWidget::initializeGL()
{

    GLfloat vertices[] = {
        // Positions          // Texture Coords
         1.0f,  1.0f, 0.0f,   1.0f, 1.0f, // Top Right
         1.0f, -1.0f, 0.0f,   1.0f, 0.0f, // Bottom Right
        -1.0f, -1.0f, 0.0f,   0.0f, 0.0f, // Bottom Left
        -1.0f,  1.0f, 0.0f,   0.0f, 1.0f  // Top Left
    };

    GLuint elements[] = {  // Note that we start from 0!
        0, 1, 3, // First Triangle (first half of the rectangle)
        1, 2, 3  // Second Triangle (second half of the rectangle)
    };

    qDebug("initializeGL()");
    initializeOpenGLFunctions();

    glClearColor( 0.0f, 0.0f, 0.0f, 1.0f );

    if (matImageRGB.type() == CV_32F || matImageRGB.type() == CV_32FC3)
    {
        if ( !prepareShaderProgram( ":/shaders/vertexShader.vert", ":/shaders/fragment32.frag") )
         return;
    }
    else
    {
        if ( !prepareShaderProgram( ":/shaders/vertexShader.vert", ":/shaders/fragment16.frag") )
         return;
    }


    // Create Vertex Array Object
    m_vao.create();
    // Create Vertex Buffer Object (vbo)
    m_vertexBuffer.create();
    // Create the "element" or "index" buffer object
    m_indexBuffer.create();

    // Bind VAO
    m_vao.bind();

    // Bind Vertex Buffer + check
    if ( !m_vertexBuffer.bind() )
    {
        qWarning() << "Could not bind vertex buffer to the context";
        return;
    }

    // Setup Vertex and index Buffer
    m_vertexBuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);

    // Send data to vertex buffer, equivalent to  glBufferData()
    m_vertexBuffer.allocate(vertices, sizeof(vertices));

    // Bind index Buffer + check
    if ( !m_indexBuffer.bind() )
    {
        qWarning() << "Could not bind vertex buffer to the context";
        return;
    }
    m_indexBuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
    // Send data to index buffer, equivalent to  glBufferData()
    m_indexBuffer.allocate(elements, sizeof(elements));

    ///----- Pure openGL commands -------///

//    glGenVertexArrays(1, &vao);
//    glBindVertexArray(vao);

//    glGenBuffers(1, &vbo);

//    glBindBuffer(GL_ARRAY_BUFFER, vbo);
//    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

//    glGenBuffers(1, &ebo);


//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
//    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(elements), elements, GL_STATIC_DRAW);

        // Bind the shader program so that we can associate variables from
        // our application to the shaders
        if ( !m_shader.bind() )
        {
            qWarning() << "Could not bind shader program to context";
            return;
        }

        m_shader.setAttributeBuffer(0, GL_FLOAT, 0, 3, 5*sizeof(GLfloat));
        m_shader.enableAttributeArray(0);

        m_shader.setAttributeBuffer(1, GL_FLOAT, 3*sizeof(GLfloat), 2, 5*sizeof(GLfloat));
        m_shader.enableAttributeArray(1);

    m_vertexBuffer.release();
    m_vao.release();

//    glBindVertexArray(0);
//    glBindBuffer(GL_ARRAY_BUFFER, 0);
//    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Prepare and load texture

    loadGLTexture();
    m_shader.release();

    GLenum glErr = GL_NO_ERROR;

    while((glErr = glGetError()) != GL_NO_ERROR)
    {
      qDebug("ROpenGLWidget::initializeGL():: ERROR INITIALISING OPENGL");
      exit(1);
    }


}

void ROpenGLWidget::loadGLTexture()
{

    for (int ii =0; ii < nFrames; ii++)
    {
        cv::Mat tempImageRGB = matImageListRGB.at(ii);
        //tempImageRGB.convertTo(tempImageRGB, CV_16UC3);

        cv::Mat tempImageGray;
        cv::cvtColor(tempImageRGB, tempImageGray, CV_RGB2GRAY);

        QOpenGLTexture *oglt = new QOpenGLTexture(QOpenGLTexture::Target2D);
        oglt->setSize(naxis1, naxis2);


        oglt->setMagnificationFilter(QOpenGLTexture::Nearest);
        //oglt->setMagnificationFilter(QOpenGLTexture::NearestMipMapNearest);
        oglt->setMinificationFilter(QOpenGLTexture::NearestMipMapNearest);


        if (tempImageRGB.type() == CV_32FC3)
        {
            qDebug("ROpenGLWidget::loadGLTexture():: image is 32 bits float, 3 channels");
            oglt->setFormat(QOpenGLTexture::RGB32F);
            oglt->allocateStorage(QOpenGLTexture::RGB, QOpenGLTexture::Float32);
            oglt->setData(QOpenGLTexture::RGB, QOpenGLTexture::Float32, tempImageRGB.data);
        }
        else if (tempImageRGB.type() == CV_16UC3)
        {
            qDebug("ROpenGLWidget::loadGLTexture():: image is 16-bit unsigned integer, 3 channels");
            oglt->setFormat(QOpenGLTexture::RGB16U);
            oglt->allocateStorage(QOpenGLTexture::RGB_Integer, QOpenGLTexture::UInt16);
            oglt->setData(QOpenGLTexture::RGB_Integer, QOpenGLTexture::UInt16, tempImageRGB.data);
        }
        else if (tempImageRGB.type() == CV_16SC3)
        {
            qDebug("ROpenGLWidget::loadGLTexture():: image is 16-bit SIGNED integer, 3 channels");
            oglt->setFormat(QOpenGLTexture::RGB16I);
            oglt->allocateStorage(QOpenGLTexture::RGB_Integer, QOpenGLTexture::Int16);
            oglt->setData(QOpenGLTexture::RGB_Integer, QOpenGLTexture::Int16, tempImageRGB.data);
        }
        else if (tempImageRGB.type() == CV_16U)
        {
            qDebug("ROpenGLWidget::loadGLTexture():: image is 16-bit unsigned integer, 1 channel");
            oglt->setFormat(QOpenGLTexture::R16U);
            oglt->allocateStorage(QOpenGLTexture::Red_Integer, QOpenGLTexture::UInt16);
            oglt->setData(QOpenGLTexture::Red_Integer, QOpenGLTexture::UInt16, tempImageRGB.data);
        }
        else
        {
            qDebug("ROpenGLWidget::loadGLTexture():: Unsupported image format");
        }
//        qDebug() << "QOpenGLTexture::isAutoMipMapGenerationEnabled() =" << oglt->isAutoMipMapGenerationEnabled();
//        qDebug() << "QOpenGLTexture::mipLevels() =" << oglt->mipLevels();

        textureVector << oglt;
    }

    GLenum glErr = GL_NO_ERROR;

    while((glErr = glGetError()) != GL_NO_ERROR)
    {
      qDebug("ROpenGLWidget::loadGLTexture():: ERROR LOADING TEXTURES");
      exit(1);
    }


//    texture = textureVector.at(frameIndex);

    ///// --------- PURE OPENGL COMMANDS -------- /////////

//    //unsigned char* image = matImageListRGB.at(0).data;

//    float borderColor[] = { 0.5f, 0.5f, 0.5f, 1.0f };

//    //glEnable(GL_TEXTURE_2D);
//    glGenTextures(1, &texture);
//    //glActiveTexture(GL_TEXTURE0);
//    glBindTexture(GL_TEXTURE_2D, texture);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
//    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

////    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST_MIPMAP_NEAREST);
////    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST_MIPMAP_NEAREST);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

//    // Load the image in the texture.
//    cv::Mat matImageRGB16;
//    matImageListRGB.at(0).convertTo(matImageRGB16, CV_16UC3);

//    QElapsedTimer timer;
//    timer.start();

//    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, naxis1, naxis2, 0, GL_RGB, GL_FLOAT, matImageListRGB.at(0).data);
//   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16UI, naxis1, naxis2, 0, GL_RGB_INTEGER, GL_UNSIGNED_SHORT, matImageRGB16.data);

////   glGenerateMipmap(GL_TEXTURE_2D);

//    qDebug()<<"ROpenGLWidget::initializeGL():: Time for opengl initialize():" <<  timer.elapsed() << "ms";

//   glBindTexture(GL_TEXTURE_2D, 0);

}


void ROpenGLWidget::paintGL()
{

    texture = textureVector.at(frameIndex);

    glClear(GL_COLOR_BUFFER_BIT);

    m_shader.bind();
    // When I just move that window around, I'm also using those shaders...
    int alphaLocation = m_shader.uniformLocation("alpha");
    int betaLocation = m_shader.uniformLocation("beta");
    int gammaLocation = m_shader.uniformLocation("gamma");

    int wbRGBLocation = m_shader.uniformLocation("wbRGB");

    m_shader.setUniformValue(alphaLocation, alpha);
    m_shader.setUniformValue(betaLocation, beta);
    m_shader.setUniformValue(gammaLocation, gamma);

    QVector3D wbRGB(wbRed, wbGreen, wbBlue);
    m_shader.setUniformValue(wbRGBLocation, wbRGB);


    m_vao.bind();
    {
        // Profile here
        texture->bind();
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        texture->release(); // unbind?
        // End Profile here
    }
    m_vao.release();
    m_shader.release(); // unbind?


//    glBindVertexArray(vao);

//    glBindTexture(GL_TEXTURE_2D, texture);
//    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
//    glBindTexture(GL_TEXTURE_2D, 0);

//    glBindVertexArray(0);
//    m_shader.release();

}

void ROpenGLWidget::resizeGL(int w, int h)
{
//    const qreal retinaScale = devicePixelRatio();
//    glViewport(0, 0, width() * retinaScale, height() * retinaScale);
      glViewport( 0, 0, w, qMax( h, 1 ) );
      qDebug("resizeGL size = (%i, %i)", w, h);
}

bool ROpenGLWidget::prepareShaderProgram(const QString vertexShaderPath, const QString fragmentShaderPath)
{
    //qDebug("Loading shaders.");
    // First we load and compile the vertex shader…
     bool result = m_shader.addShaderFromSourceFile( QOpenGLShader::Vertex, vertexShaderPath );
     if ( !result )
     qWarning() << m_shader.log();

    // fragment shader
     result = m_shader.addShaderFromSourceFile( QOpenGLShader::Fragment, fragmentShaderPath );
     if ( !result )
     qWarning() << m_shader.log();

    // link the shaders
     result = m_shader.link();
     if ( !result )
     qWarning() << "Could not link shader program:" << m_shader.log();

    return result;
}


void ROpenGLWidget::mousePressEvent(QMouseEvent *event)
{
    setCursor(Qt::CrossCursor);
    //Coordinates in pixels in the local (resized) Widget
    // They are not in image coordinates.
    int cursorX = event->x();
    int cursorY = event->y();
    //qDebug("Mouse at (%i, %i)", cursorX, cursorY);

    imageCoordX = round((float) (cursorX) / (float) (resizeFac));
    imageCoordY = round((float) (cursorY) / (float) (resizeFac));
    qDebug("Image coordinates: (%i, %i)", imageCoordX, imageCoordY);

    updateSubQImage();

    emit mousePressed();
}

void ROpenGLWidget::mouseMoveEvent(QMouseEvent *event)
{
    setCursor(Qt::CrossCursor);
    //Coordinates in pixels in the local (resized) Widget
    // They are not in image coordinates.
    int cursorX = event->x();
    int cursorY = event->y();

    imageCoordX = cursorX / resizeFac;
    imageCoordY = cursorY / resizeFac;

    updateSubQImage();
    emit mousePressed();
}

void ROpenGLWidget::setupTableWidget(int value)
{
    if (rMatImageList.empty() || rMatImageList.at(value)->isBayer())
    {
        return;
    }

    rListImageManager->getTableWidgetList().at(value)->clearSelection();
    tableRSubWindow->setWidget(rListImageManager->getTableWidgetList().at(value));
    tableRSubWindow->setWindowTitle(windowTitleList.at(value));
    tableRSubWindow->resize(tableSize);
    rListImageManager->getTableWidgetList().at(value)->clearSelection();
    rListImageManager->getTableWidgetList().at(value)->setCurrentCell(0, 0);

    emit sendNewTitle(windowTitleList.at(value));

}

void ROpenGLWidget::changeFrame(int value)
{
    this->frameIndex = value;
    updateSubQImage();
    update();

    emit sendNewTitle(windowTitleList.at(value));
}


void ROpenGLWidget::wheelEvent(QWheelEvent *wheelEvent)
{
    if (wheelEvent->modifiers() & Qt::ShiftModifier)
    {
        qreal numDegrees = (qreal) (wheelEvent->angleDelta().y() / 8.0f);
        qDebug() << "numDegrees =" << -numDegrees;
        qreal zoomScale = (qreal) (-numDegrees/100.0f);
        qDebug() << "zoomScale =" << zoomScale;
        QSize oglSize = this->size();
        oglSize *= (1 + zoomScale);
        //qDebug() << "oglSize =" << oglSize;
        if (oglSize.width() > 10*naxis1 || oglSize.height() > 10*naxis2)
        {
            return;
        }

        if (oglSize.width() < naxis1/10 || oglSize.height() < naxis2/10)
        {
            return;
        }

        resize(oglSize);
        // Update the resizeFac
        resizeFac = ((float) (oglSize.width())/ (float) (naxis1));

        qDebug() << "wheelEvent::resizeFac =" << resizeFac;

    }
    else
    {
        QWidget::wheelEvent(wheelEvent);
    }

}

void ROpenGLWidget::initSubQImage()
{

    subNaxis = 60;
    int frameRows = naxis2 + subNaxis +1;
    int frameCols = naxis1 + subNaxis +1;

    matFrame = cv::Mat(frameRows, frameCols, CV_32FC3, cv::Scalar(32000, 32000, 32000));
    subQImage = new QImage(subNaxis, subNaxis, QImage::Format_ARGB32);
}

void ROpenGLWidget::updateSubQImage()
{
    matImageRGB = matImageListRGB.at(frameIndex);
    matImageRGB.copyTo(matFrame(cv::Rect(subNaxis/2-1, subNaxis/2-1, naxis1, naxis2)));


    float alpha2 = (float) (255.0f * alpha);
    float beta2 = (float) (255.0f * beta);

    // Convert the cursor coord. in the large frame coord. system
    int x = (imageCoordX - subNaxis/2) + subNaxis/2 ;
    int y = (naxis2 - imageCoordY) - (subNaxis/2 + 1) + subNaxis/2 ;

    // Define the Xmax, Ymax value that the coordinate of the top left rectangle vertex must not exceed.
    // And constrain that vertex coordinate  (x,y) to be within [0 ; Xmax] and [0 ; Ymax].
    int Xmax = matFrame.cols-1 - (subNaxis-1);
    int Ymax = matFrame.rows-1 - (subNaxis-1);
    FOV.x = std::max(std::min(x, Xmax), 0);
    FOV.y = std::max(std::min(y, Ymax), 0);
    FOV.width = subNaxis;
    FOV.height = subNaxis;

    cv::Mat tempMat = matFrame(FOV);
    tempMat.convertTo(subMatImage, CV_8UC3, alpha2, beta2);

     for (int i=0; i<subNaxis; i++)
     {
         for (int j=0; j < subNaxis; j++)
         {
             cv::Vec3b color = subMatImage.at<cv::Vec3b>(i, j);
             int red = (int) (color.val[0]);
             int green = (int) (color.val[1]);
             int blue = (int) (color.val[2]);
             QRgb pcolor = qRgb(red, green, blue);
             subQImage->setPixel(j, i, pcolor);
         }
     }

     // update info
     x = imageCoordX;
     y = naxis2 - (int) (imageCoordY -1);

     x = std::max(std::min(x, naxis1-1), 0);
     y = std::max(std::min(y, naxis2-1), 0);

     if (matImageRGB.type() == CV_32FC3)
     {
         cv::Vec3f color = matImageRGB.at<cv::Vec3f>(y, x);
         this->intensity = (float) ((color.val[0] + color.val[1] + color.val[2])/3.0f);
     }
     else
     {
         cv::Vec3w color = matImageRGB.at<cv::Vec3w>(y, x);
         this->intensity = (float) ((color.val[0] + color.val[1] + color.val[2])/3.0f);
     }

     emit sendSubQImage(subQImage, intensity, x, y);
}


void ROpenGLWidget::focusInEvent(QFocusEvent * event)
{
    qDebug("ROpenGLWidget:: emitting gotSelected(this) signal");
    emit gotSelected(this);
}

void ROpenGLWidget::focusOutEvent(QFocusEvent * event)
{

}

// Getters


int ROpenGLWidget::getFrameIndex()
{
    return frameIndex;
}

quint32 ROpenGLWidget::getImageCoordX()
{
    return imageCoordX;
}

quint32 ROpenGLWidget::getImageCoordY()
{
    return imageCoordY;
}

RListImageManager *ROpenGLWidget::getRListImageManager()
{
    return rListImageManager;
}

QList<RMat*> ROpenGLWidget::getRMatImageList()
{
    return rMatImageList;
}


QString ROpenGLWidget::getCalibrationType()
{
    return calibrationType;
}

QList<QString> ROpenGLWidget::getWindowTitleList()
{
    return windowTitleList;
}

QString ROpenGLWidget::getWindowTitle()
{
    return windowTitle;
}

RSubWindow *ROpenGLWidget::getTableRSubWindow()
{
    return tableRSubWindow;
}

float ROpenGLWidget::getGamma()
{
    return gamma;
}

float ROpenGLWidget::getNewMax()
{
    return newMax;
}

float ROpenGLWidget::getNewMin()
{
    return newMin;
}

float ROpenGLWidget::getWbRed()
{
    return wbRed;
}

float ROpenGLWidget::getWbGreen()
{
    return wbGreen;
}

float ROpenGLWidget::getWbBlue()
{
    return wbBlue;
}

void ROpenGLWidget::setNewMax(float newMax)
{
    this->newMax = newMax;
}

void ROpenGLWidget::setNewMin(float newMin)
{
    this->newMin = newMin;
}

// Setters


void ROpenGLWidget::setAlpha(float newAlpha)
{
    alpha = newAlpha;
}

void ROpenGLWidget::setBeta(float newBeta)
{
    beta = newBeta;
}

void ROpenGLWidget::setGamma(float newGamma)
{
    gamma = newGamma;
}


void ROpenGLWidget::setFrameIndex(int newFrameIndex)
{
    frameIndex = newFrameIndex;
}

void ROpenGLWidget::setCalibrationType(QString calibrationType)
{
    this->calibrationType = calibrationType;
}

void ROpenGLWidget::setWindowTitleList(QList<QString> newWindowTitleList)
{
    this->windowTitleList = newWindowTitleList;
}

void ROpenGLWidget::setWindowTitle(QString newWindowTitle)
{
    this->windowTitle = newWindowTitle;
}

void ROpenGLWidget::setResizeFac(float resizeFac)
{
    this->resizeFac = resizeFac;
}

void ROpenGLWidget::setWbRed(float wbRed)
{
    this->wbRed = wbRed;
}

void ROpenGLWidget::setWbGreen(float wbGreen)
{
    this->wbGreen = wbGreen;
}

void ROpenGLWidget::setWbBlue(float wbBlue)
{
    this->wbBlue = wbBlue;
}

void ROpenGLWidget::setTableSize(QSize size)
{
    this->tableSize = size;
}

