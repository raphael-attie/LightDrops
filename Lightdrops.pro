#-------------------------------------------------
#
# Project created by QtCreator 2016-02-16T21:10:09
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Lightdrops
TEMPLATE = app

INCLUDEPATH += /opt/local/include
INCLUDEPATH += /Usr/local/include/libraw

INCLUDEPATH += ~/Dev/opencv3_tbb_opencl/include
INCLUDEPATH += /usr/local/include

LIBS += -L/opt/local/lib -lcfitsio
LIBS += -L/usr/local/lib -lraw -lafopencl
LIBS += -L../opencv3_tbb_opencl/lib -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_video


# Setup Qt so Clang works with C++11
LIBS += -stdlib=libc++

QMAKE_CXXFLAGS += -stdlib=libc++
QMAKE_CXXFLAGS += -std=c++11
QMAKE_CXXFLAGS += -mmacosx-version-min=10.7
QMAKE_LFLAGS += -mmacosx-version-min=10.7



SOURCES += main.cpp\
        rmainwindow.cpp \
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
    rsubwindow.cpp

HEADERS  += rmainwindow.h \
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
    rsubwindow.h

FORMS    += rmainwindow.ui

RESOURCES += \
    ressources.qrc

DISTFILES +=
