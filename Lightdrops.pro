#-------------------------------------------------
#
# Project created by QtCreator 2016-02-16T21:10:09
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Lightdrops
TEMPLATE = app

win32 {
    DEFINES += WIN32
    INCLUDEPATH += C:/dev/cfitsio_64
    INCLUDEPATH += C:/dev/libraw/libraw
    INCLUDEPATH += C:/dev/opencv/build/include
    INCLUDEPATH += C:\dev\ArrayFire\v3\include

    LIBS += -LC:\dev\cfitsio_64 -lcfitsio
    LIBS += -LC:\dev\libraw\x64\Debug -llibraw
    LIBS += -LC:\dev\opencv\build\x64\vc14\lib -lopencv_world310
    LIBS += -LC:/dev/ArrayFire/v3/lib -lafopencl
}

macx {
    INCLUDEPATH += /opt/local/include /usr/local/include/libiomp
    INCLUDEPATH += /Users/raphaela/Dev/Libraw_source/LibRaw-0.17.0/install_noAddDemosaic/include/libraw
    INCLUDEPATH += /Users/raphaela/Dev/opencv3_tbb_opencl/include
    INCLUDEPATH += /usr/local/include

    LIBS += -L/opt/local/lib -lcfitsio
    LIBS += -L../Libraw_source/LibRaw-0.17.0/install_noAddDemosaic/lib -lraw
    LIBS += -L../opencv3_tbb_opencl/lib -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_video
    LIBS += -L/usr/local/lib -lafopencl
    QMAKE_CXXFLAGS += -mmacosx-version-min=10.7
    QMAKE_LFLAGS += -mmacosx-version-min=10.7
    QMAKE_CXXFLAGS += -stdlib=libc++
    # Setup Qt so Clang works with C++11
    LIBS += -stdlib=libc++
    QMAKE_CXXFLAGS += -std=c++11
}

SOURCES += main.cpp\
        rmainwindow.cpp \
    scrollarea.cpp \
    ropenglwidget.cpp \
    MyFitsImage.cpp \
    RawImage.cpp \
    rmat.cpp \
    rlistimagemanager.cpp \
    imagemanager.cpp \
    rtreewidget.cpp \
    rprocessing.cpp \
    parallelcalibration.cpp \
    rlineedit.cpp \
    RFrame.cpp \
    qcustomplot/qcustomplot.cpp \
    rtableworker.cpp

HEADERS  += rmainwindow.h \
    scrollarea.h \
    ropenglwidget.h \
    MyFitsImage.h \
    RawImage.h \
    rmat.h \
    rlistimagemanager.h \
    imagemanager.h \
    rtreewidget.h \
    rprocessing.h \
    parallelcalibration.h \
    rlineedit.h \
    RFrame.h \
    qcustomplot/qcustomplot.h \
    rtableworker.h \
    winsockwrapper.h \
    myunistd.h

FORMS    += rmainwindow.ui \
    scrollarea.ui

RESOURCES += \
    shaders.qrc

DISTFILES +=
