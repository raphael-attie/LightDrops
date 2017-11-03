#include "winsockwrapper.h"
#include "rmainwindow.h"
#include <QApplication>

// Arrayfire
#include <arrayfire.h>

int main(int argc, char *argv[])
{

    qDebug() << af::infoString();


    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setVersion(4, 4);
    QSurfaceFormat::setDefaultFormat(format);


    QApplication a(argc, argv);

    //a.setStyle("fusion");

    RMainWindow w;
    w.show();

    return a.exec();
}
