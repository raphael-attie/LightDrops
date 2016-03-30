#include "winsockwrapper.h"
#include "rmainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    //format.setVersion(3, 3);
    format.setVersion(3, 3);
    QSurfaceFormat::setDefaultFormat(format);

    QApplication a(argc, argv);

    //a.setStyle("fusion");

    qDebug() << "opengGL version:" << format.version();

//    cv::Mat ellMat = cv::Mat::zeros(500,500, CV_8U);
//    QImage image2(ellMat.data, ellMat.cols, ellMat.rows, QImage::Format_Grayscale8);

//    QGraphicsPixmapItem item( QPixmap::fromImage(image2));
//    QGraphicsScene* scene = new QGraphicsScene;
//    scene->addItem(&item);

//    QGraphicsView* graphicsView= new QGraphicsView(scene);
//    graphicsView->show();


    RMainWindow w;
    w.show();

    return a.exec();
}
