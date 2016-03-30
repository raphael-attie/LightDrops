#-------------------------------------------------
#
# Project created by QtCreator 2016-02-16T21:10:09
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = Lightdrops
TEMPLATE = app

<<<<<<< Updated upstream
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

    INCLUDEPATH += /opt/local/include
    INCLUDEPATH += /Usr/local/include/libraw
    INCLUDEPATH += /usr/local/include/libiomp

#choose open cv setup
#    INCLUDEPATH += ~/Dev/opencv3_tbb_opencl/include
    INCLUDEPATH += /Usr/local/include

    LIBS += -L/opt/local/lib -lcfitsio
    LIBS += -L/usr/local/lib -lraw -lafopencl

#choose open cv setup
#    LIBS += -L../opencv3_tbb_opencl/lib -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_video
    LIBS += -L/usr/local/lib -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_video

    QMAKE_CXXFLAGS += -mmacosx-version-min=10.7
    QMAKE_LFLAGS += -mmacosx-version-min=10.7
    QMAKE_CXXFLAGS += -stdlib=libc++
    # Setup Qt so Clang works with C++11
    LIBS += -stdlib=libc++
    QMAKE_CXXFLAGS += -std=c++11
}
=======
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


>>>>>>> Stashed changes

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
    rsubwindow.cpp \
    data.cpp \
    circle.cpp \
    utilities.cpp

HEADERS  += winsockwrapper.h \
    rmainwindow.h \
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
<<<<<<< Updated upstream
    rsubwindow.h
    rtableworker.h \

=======
    rsubwindow.h \
    data.h \
    typedefs.h \
    circle.h \
    utilities.h
>>>>>>> Stashed changes

FORMS    += rmainwindow.ui

RESOURCES += \
    ressources.qrc

DISTFILES +=
