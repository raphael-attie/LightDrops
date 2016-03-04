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

    RMainWindow w;
    w.show();

    return a.exec();
}
