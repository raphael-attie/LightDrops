#include "winsockwrapper.h"
#include "rmainwindow.h"
#include <QApplication>

// Arrayfire
#include <arrayfire.h>

int main(int argc, char *argv[])
{

//    string filePath("/Users/raphaela/Pictures/Astrophotography/Joris/Andromede_iso3200_20s/F36A7292.CR2");
//    ifstream inputF(filePath, ios::binary);
//    inputF >> noskipws;
//    qDebug() << "position = " << inputF.tellg();
    //istream_iterator<unsigned char> it(inputF);
    //inputF.seekg(0, ios::beg);

//    char buff;
//    qDebug() << "sizeof(buff) = " << sizeof(buff);
//    int i = 0;

//    while (inputF.tellg() < 60)
//    {
//        qDebug() << "inputF.tellg() = " << inputF.tellg();
//        inputF.read( reinterpret_cast<char *> (&buff), sizeof(buff));
//        i++;
//        qDebug() << "buff = " << buff;
//    }


//    cout << hex << setw(2) << setfill('0') << (int)*it << endl;
//    ++it;
//    qDebug() << "position = " << inputF.tellg();
//    cout << hex << setw(2) << setfill('0') << (int)*it << endl;
//    ++it;
//    qDebug() << "position = " << inputF.tellg();
//    cout << hex << setw(2) << setfill('0') << (int)*it << endl;

//    cout << *it;
//    qDebug() << "position = " << inputF.tellg();
//    ++it;
//    cout << *it;
//    qDebug() << "position = " << inputF.tellg();
//    ++it;
//    cout << *it << endl;

/////////////

//    char buff;
//    unsigned short number;

//    inputF.seekg(0, ios::beg);
//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( &buff, 1);
//    qDebug() << "buff = " << buff;

//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( &buff, 1);
//    qDebug() << "buff = " << buff;

//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( &buff, 1);
//    qDebug() << "buff = " << buff;



//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( reinterpret_cast<char *> (&number), 1);
//    qDebug() << "number = " << number;

//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( reinterpret_cast<char *> (&number), 1);
//    qDebug() << "number = " << number;
//    cout << "number (hex offset to first IFD) = " << hex << showbase << number << endl;

//    inputF.seekg(8, ios::beg);
//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( &buff, 1);
//    qDebug() << "buff = " << buff;

//    qDebug() << "position = " << inputF.tellg();
//    inputF.read( &buff, 1);
//    qDebug() << "buff = " << buff;

    qDebug() << af::infoString();


    QSurfaceFormat format;
    format.setDepthBufferSize(24);
    format.setStencilBufferSize(8);
    format.setProfile(QSurfaceFormat::CoreProfile);
    format.setVersion(3, 3);
    QSurfaceFormat::setDefaultFormat(format);


    QApplication a(argc, argv);

    //a.setStyle("fusion");

    RMainWindow w;
    w.show();

    return a.exec();
}
