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

/// opencv
#include <opencv2/world.hpp>

ROpenGLWidget::ROpenGLWidget(QWidget *parent) :
    QOpenGLWidget(parent)
{

}

ROpenGLWidget::ROpenGLWidget(RListImageManager *rListImageManager, QWidget *parent)
    :QOpenGLWidget(parent)
    , m_shader(0)
    , m_vertexBuffer(QOpenGLBuffer::VertexBuffer)
    , m_indexBuffer(QOpenGLBuffer::IndexBuffer)
    , tableRSubWindow(NULL), tableWidget(NULL)
{
    this->rListImageManager = rListImageManager;
    this->rMatImageList = rListImageManager->getRMatImageList();
    this->nFrames = rMatImageList.size();

    initialize();
    //setupTableWidget();
    tableRSubWindow = new RSubWindow();
    tableSize = rListImageManager->getTableWidgetList().at(0)->size();

}

ROpenGLWidget::ROpenGLWidget(QList<RMat *> rMatImageList, QWidget *parent)
    :QOpenGLWidget(parent), rListImageManager(NULL)
    , m_shader(0)
    , m_vertexBuffer(QOpenGLBuffer::VertexBuffer)
    , m_indexBuffer(QOpenGLBuffer::IndexBuffer)
    , tableRSubWindow(NULL), tableWidget(NULL)
{
    this->rMatImageList = rMatImageList;
    this->nFrames = rMatImageList.size();

    initialize();
}

ROpenGLWidget::ROpenGLWidget(RMat *rMatImage, QWidget *parent)
    :QOpenGLWidget(parent), rListImageManager(NULL)
    , m_shader(0)
    , m_vertexBuffer(QOpenGLBuffer::VertexBuffer)
    , m_indexBuffer(QOpenGLBuffer::IndexBuffer)
    , tableRSubWindow(NULL), tableWidget(NULL)
{
    this->rMatImageList << rMatImage;
    //this->rMatImageList.push_back(rMatImage);
    this->nFrames = 1;

    initialize();
}

ROpenGLWidget::~ROpenGLWidget()
{
    qDebug("ROpenGLWidget::~ROpenGLWidget() destructor");
    m_indexBuffer.destroy();
    m_vertexBuffer.destroy();
    m_vao.destroy();

    delete rListImageManager;


    /// Deleting the rMatImageList below is causing a crash.
    ///
//    if (!rMatImageList.isEmpty())
//    {
//        qDeleteAll(rMatImageList);
//        rMatImageList.clear();
//    }

}


void ROpenGLWidget::initialize()
{

    customPlot = NULL;
    vertLineHigh = NULL;
    vertLineLow = NULL;

    setFocusPolicy(Qt::StrongFocus);
    frameIndex = 0;
    counter = 0;
    dataMax = rMatImageList.at(frameIndex)->getDataMax();
    dataMin = rMatImageList.at(frameIndex)->getDataMin();
    dataRange = (float) (dataMax - dataMin);
    newMin = dataMin;
    newMax = dataMax;

    calibrationType = QString("");
    naxis1 = rMatImageList.at(frameIndex)->matImage.cols;
    naxis2 = rMatImageList.at(frameIndex)->matImage.rows;

    intensity = 0;
    alpha = 1.0f/dataRange;
    beta = (float) (-dataMin / dataRange);
    gamma = 1.0;

    // White balance default
    wbRed = 1.0f;
    wbGreen = 1.0f;
    wbBlue = 1.0f;

    iMax = 5000.0;
    lambda = 100.0;
    mu = 100.0;
    applyToneMapping = false;

    imageCoordX = naxis1/2;
    imageCoordY = naxis2/2;

    prepImage();
    matImageRGB = matImageListRGB.at(frameIndex);
    initSubQImage();
    setupHistoPlots();

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
//    if (!matImageListRGB.isEmpty())
//    {
//        matImageListRGB.clear();
//    }

    for (int ii = 0; ii < nFrames; ii++)
    {
        cv::Mat tempMatRGB(naxis2, naxis1, rMatImageList.at(ii)->matImage.type());

        if (rMatImageList.at(ii)->isBayer())
        {
            cv::Mat tempMat16;
            rMatImageList.at(ii)->matImage.convertTo(tempMat16, CV_16U);
            cv::cvtColor(tempMat16, tempMatRGB, CV_BayerBG2RGB);
        }
        else
        {
            tempMatRGB = rMatImageList.at(ii)->matImage.clone();
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
    else if (matImageRGB.type() == CV_16UC3 || matImageRGB.type() == CV_16SC3)
    {
        if ( !prepareShaderProgram( ":/shaders/vertexShader.vert", ":/shaders/fragment16.frag") )
         return;
    }
    else if (matImageRGB.type() == CV_16U)
    {
        if ( !prepareShaderProgram( ":/shaders/vertexShader.vert", ":/shaders/fragment16uc1.frag") )
         return;
    }
    else if (matImageRGB.type() == CV_8UC3)
    {
        if ( !prepareShaderProgram( ":/shaders/vertexShader.vert", ":/shaders/fragment888.frag") )
         return;
    }
    else if (matImageRGB.type() == CV_8U)
    {
        if ( !prepareShaderProgram( ":/shaders/vertexShader.vert", ":/shaders/fragment8.frag") )
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
//    if (!textureVector.isEmpty())
//    {
//        qDeleteAll(textureVector);
//        textureVector.clear();
//    }

    for (int ii =0; ii < nFrames; ii++)
    {
        cv::Mat tempImageRGB = matImageListRGB.at(ii);

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
        else if (tempImageRGB.type() == CV_8UC3)
        {
            qDebug("ROpenGLWidget::loadGLTexture():: image is 8-bit unsigned integer, 3 channels");
            oglt->setFormat(QOpenGLTexture::RGB8U);
            oglt->allocateStorage(QOpenGLTexture::RGB_Integer, QOpenGLTexture::UInt8);
            oglt->setData(QOpenGLTexture::RGB_Integer, QOpenGLTexture::UInt8, tempImageRGB.data);
        }
        else if (tempImageRGB.type() == CV_8U)
        {
            qDebug("ROpenGLWidget::loadGLTexture():: image is 8-bit unsigned integer, 1 channel");
            oglt->setFormat(QOpenGLTexture::R8U);
            oglt->allocateStorage(QOpenGLTexture::Red_Integer, QOpenGLTexture::UInt8);
            oglt->setData(QOpenGLTexture::Red_Integer, QOpenGLTexture::UInt8, tempImageRGB.data);
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

}


void ROpenGLWidget::paintGL()
{

    texture = textureVector.at(frameIndex);

    glClear(GL_COLOR_BUFFER_BIT);

    m_shader.bind();
    /// When I just move that window around, I'm also using those shaders...
    int alphaLocation = m_shader.uniformLocation("alpha");
    int betaLocation = m_shader.uniformLocation("beta");
    int gammaLocation = m_shader.uniformLocation("gamma");
    int iMaxLocation = m_shader.uniformLocation("iMax");
    int lambdaLocation = m_shader.uniformLocation("lambda");
    int muLocation = m_shader.uniformLocation("mu");
    int applyToneMappingLocation = m_shader.uniformLocation("applyToneMapping");

    int wbRGBLocation = m_shader.uniformLocation("wbRGB");

    m_shader.setUniformValue(alphaLocation, alpha);
    m_shader.setUniformValue(betaLocation, beta);
    m_shader.setUniformValue(gammaLocation, gamma);
    /// Set ToneMapping parameters
    m_shader.setUniformValue(iMaxLocation, iMax);
    m_shader.setUniformValue(lambdaLocation, lambda);
    m_shader.setUniformValue(muLocation, mu);
    m_shader.setUniformValue(applyToneMappingLocation, applyToneMapping);



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
    int channels = matImageRGB.channels();

    if (channels == 3)
    {
        matFrame = cv::Mat(frameRows, frameCols, CV_32FC3, cv::Scalar(32000, 32000, 32000));
    }
    else if (channels == 1)
    {
        matFrame = cv::Mat(frameRows, frameCols, CV_32F, cv::Scalar(32000));
    }

    subQImage = new QImage(subNaxis, subNaxis, QImage::Format_ARGB32_Premultiplied);

}

void ROpenGLWidget::updateSubQImage()
{

    matImageRGB = matImageListRGB.at(frameIndex);
    matImageRGB.copyTo(matFrame(cv::Rect(subNaxis/2-1, subNaxis/2-1, naxis1, naxis2)));

    int channels = matFrame.channels();

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
    tempMat.convertTo(subMatImage, CV_8UC(channels), alpha2, beta2);

    int red, green, blue;
     for (uint i=0; i<subNaxis; i++)
     {
         for (uint j=0; j < subNaxis; j++)
         {

             if (channels == 3)
             {
                 cv::Vec3b color = subMatImage.at<cv::Vec3b>(i, j);
                 red = (int) (color.val[0]);
                 green = (int) (color.val[1]);
                 blue = (int) (color.val[2]);
             }
             else
             {
                 cv::Scalar color = subMatImage.at<uchar>(i, j);
                 red = (int) std::min(wbRed * (float) color.val[0], 255.0f);
                 green = (int) std::min(wbBlue * (float) color.val[0], 255.0f);
                 blue = (int) std::min(wbGreen * (float) color.val[0], 255.0f);
             }

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
     else if (matImageRGB.type() == CV_16UC3)
     {
         cv::Vec3w color = matImageRGB.at<cv::Vec3w>(y, x);
         this->intensity = (float) ((color.val[0] + color.val[1] + color.val[2])/3.0f);
     }
     else if (matImageRGB.type() == CV_16U)
     {
         cv::Scalar color = matImageRGB.at<ushort>(y, x);
         this->intensity = (float) color.val[0];
     }
     else if (matImageRGB.type() == CV_8UC3)
     {
         cv::Vec3b color = matImageRGB.at<cv::Vec3b>(y, x);
         this->intensity = (float) ((color.val[0] + color.val[1] + color.val[2])/3.0f);
     }
     else if (matImageRGB.type() == CV_8UC1)
     {
         cv::Scalar color = matImageRGB.at<uchar>(y, x);
         this->intensity = (float) color.val[0];
     }


     emit sendSubQImage(subQImage, intensity, x, y);
}

void ROpenGLWidget::setupHistoPlots()
{
    // Create list of QcustomPlot

    for (int i = 0 ; i < nFrames ; i++)
    {
        customPlotList << new QCustomPlot();
        vertLineHigh = new QCPItemLine(customPlotList.at(i));
        vertLineLow = new QCPItemLine(customPlotList.at(i));
        itemTextHigh = new QCPItemText(customPlotList.at(i));
        itemTextLow = new QCPItemText(customPlotList.at(i));

        vertLineHighList << vertLineHigh;
        vertLineLowList << vertLineLow;

        itemTextHighList << itemTextHigh;
        itemTextLowList << itemTextLow;

        cv::Mat tempHist;
        cv::Mat matHist = rMatImageList.at(i)->getMatHist();
        matHist.convertTo(tempHist, CV_64F);

        int nBins = matHist.rows;
        double minHist;
        double maxHist;
        cv::minMaxLoc(matHist, &minHist, &maxHist);
        double histWidth = rMatImageList.at(i)->getHistWidth();


        // This needs to be fully consistent with RMat::calcStats()
        double minXRange = rMatImageList.at(i)->getMinHistRange();
        double maxXRange = rMatImageList.at(i)->getMaxHistRange();

        QVector<double> x(nBins);

        for (int ii=0; ii< nBins; ++ii)
        {
            // This needs to be fully consistent with RMat::calcStats()
            // x goes from 0 to nBins-1
            // y must be the histogram values
            x[ii] = minXRange + histWidth*ii;
        }

        maxXRange = 1.01 * rMatImageList.at(i)->getMaxHistRange();

        QVector<double> xTicks(2);
        xTicks[0] = std::min(0.0f, rMatImageList.at(i)->getMinHistRange());
        xTicks[1] = rMatImageList.at(i)->getMaxHistRange();

        // Pointer to matHist data.
        const double* matHistPtr = tempHist.ptr<double>(0);
        std::vector<double> matHistStdVect(matHistPtr, matHistPtr + matHist.rows);
        QVector<double> y = QVector<double>::fromStdVector(matHistStdVect);

        customPlotList.at(i)->addGraph();
        customPlotList.at(i)->plottable(0)->setPen(QPen(QColor(125, 125, 125, 50))); // line color gray
        customPlotList.at(i)->plottable(0)->setBrush(QBrush(QColor(125, 125, 125, 50)));

        // setup axis
        //customPlotList.at(i)->xAxis->setTicks(false);
        //customPlotList.at(i)->xAxis->setTickLabels(false);
        customPlotList.at(i)->xAxis->setAutoTicks(false);
        customPlotList.at(i)->xAxis->setTickVector(xTicks);

        customPlotList.at(i)->yAxis->setTicks(false);

        //customPlotList.at(i)->yAxis->setTickLabels(false);
        // Here we have put two plot axes. So we need to make left and bottom axes always transfer their ranges to right and top axes:
        connect(customPlotList.at(i)->xAxis, SIGNAL(rangeChanged(QCPRange)), customPlotList.at(i)->xAxis2, SLOT(setRange(QCPRange)));
        connect(customPlotList.at(i)->yAxis, SIGNAL(rangeChanged(QCPRange)), customPlotList.at(i)->yAxis2, SLOT(setRange(QCPRange)));

        customPlotList.at(i)->yAxis->setScaleType(QCPAxis::stLogarithmic);
        customPlotList.at(i)->graph(0)->setData(x, y);
        customPlotList.at(i)->rescaleAxes();


        customPlotList.at(i)->xAxis->setRange(0, 1.05*maxXRange);
        customPlotList.at(i)->yAxis->setRange(1, maxHist);

        // Allow the user to zoom in / zoom out and scroll horizontally
        customPlotList.at(i)->setInteraction(QCP::iRangeDrag, true);
        customPlotList.at(i)->setInteraction(QCP::iRangeZoom, true);
        customPlotList.at(i)->axisRect(0)->setRangeDrag(Qt::Horizontal);
        customPlotList.at(i)->axisRect(0)->setRangeZoom(Qt::Horizontal);


        double nPixels = (double) rMatImageList.at(i)->getNPixels();
        vertLineHigh->start->setCoords(rMatImageList.at(i)->getDataMax(), 0);
        vertLineHigh->end->setCoords(rMatImageList.at(i)->getDataMax(), nPixels);
        vertLineLow->start->setCoords(rMatImageList.at(i)->getDataMin(), 0);
        vertLineLow->end->setCoords(rMatImageList.at(i)->getDataMin(), nPixels);

        // Text item
        // Max Threshold
        itemTextHigh->setClipToAxisRect(false);
        itemTextHigh->setRotation(90);
        // Min Threshold
        itemTextLow->setClipToAxisRect(false);
        itemTextLow->setRotation(90);

        customPlotList.at(i)->resize(200, 100);
    }


}

void ROpenGLWidget::updateCustomPlotLineItems()
{
    double xRangeUpper = customPlotList.at(frameIndex)->xAxis->range().upper;
    double yRangeUpper = customPlotList.at(frameIndex)->yAxis->range().upper;
    //double xRangeLower = customPlotList.at(frameIndex)->xAxis->range().lower;

    double nPixels = (double) rMatImageList.at(frameIndex)->getNPixels();
    vertLineHighList.at(frameIndex)->start->setCoords(newMax, 0);
    vertLineHighList.at(frameIndex)->end->setCoords(newMax, nPixels);
    vertLineLowList.at(frameIndex)->start->setCoords(newMin, 0);
    vertLineLowList.at(frameIndex)->end->setCoords(newMin, nPixels);

    // Max Threshold
    itemTextHighList.at(frameIndex)->position->setCoords(newMax + 0.05*xRangeUpper, 0.1*yRangeUpper);
    QString maxQString = QString::number(100.0f * newMax / rMatImageList.at(frameIndex)->getMaxHistRange(), 'f', 0) + QString("%");
    itemTextHighList.at(frameIndex)->setText(maxQString);
    // Min Threshold
    QString minQString = QString::number(100.0f * newMin / rMatImageList.at(frameIndex)->getMaxHistRange(), 'f', 0) + QString("%");
    itemTextLowList.at(frameIndex)->position->setCoords(newMin + 0.05*xRangeUpper, 0.1*yRangeUpper);
    itemTextLowList.at(frameIndex)->setText(minQString);

    customPlotList.at(frameIndex)->replot();
}

void ROpenGLWidget::focusInEvent(QFocusEvent * event)
{
    event->accept();
    qDebug("ROpenGLWidget:: emitting gotSelected(this) signal");
    emit gotSelected(this);
}

void ROpenGLWidget::focusOutEvent(QFocusEvent * event)
{
    event->accept();
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

float ROpenGLWidget::getResizeFac()
{
    return resizeFac;
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

QList<QCustomPlot*> ROpenGLWidget::getCustomPlotList()
{
    return customPlotList;
}

QCustomPlot *ROpenGLWidget::fetchCurrentCustomPlot()
{
    return customPlotList.at(frameIndex);
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

float ROpenGLWidget::getAlpha()
{
    return alpha;
}

float ROpenGLWidget::getBeta()
{
    return beta;
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

void ROpenGLWidget::setRMatImageList(QList<RMat *> rMatImageList)
{
    this->rMatImageList = rMatImageList;
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

void ROpenGLWidget::setImax(float iMax)
{
    this->iMax = iMax;
}

void ROpenGLWidget::setLambda(float lambda)
{
    this->lambda = lambda;
}

void ROpenGLWidget::setMu(float mu)
{
    this->mu = mu;
}

void ROpenGLWidget::setApplyToneMapping(bool state)
{
    applyToneMapping = state;
}

void ROpenGLWidget::setImageCoordX(quint32 x)
{
    imageCoordX = x;
}

void ROpenGLWidget::setImageCoordY(quint32 y)
{
   imageCoordY = y;
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

