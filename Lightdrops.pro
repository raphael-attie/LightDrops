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
    INCLUDEPATH += C:/dev/libraw
    INCLUDEPATH += C:/dev/opencv/include
    INCLUDEPATH += C:\dev\ArrayFire\v3\include

    LIBS += -LC:\dev\cfitsio_64 -lcfitsio
    LIBS += -LC:/dev/libraw/lib -llibraw
    LIBS += -LC:/dev/ArrayFire/v3/lib -lafopencl
    LIBS += -LC:/dev/opencv/lib -lopencv_world310

}

macx {

    INCLUDEPATH += /usr/local/include
    INCLUDEPATH += /opt/local/include

# cfitsio
    LIBS += -L/usr/local/lib -lcfitsio

# libraw
    LIBS += -L/usr/local/lib -lraw
# exiv2
    LIBS += -L/opt/local/lib -lexiv2

# opencv "world" (it has all modules)
    LIBS += -L/usr/local/lib -lopencv_world

# Caveats regarding linking options:
# In Projects>run menu
# Original DYLD_FRAMEWORK_PATH = ~/Qt/5.7/clang_64/lib:/usr/local/lib
# must change to ~/Qt/5.7/clang_64/lib (remove e.g. /usr/local/lib if it's there)
# DYLD_LIBRARY_PATH must be removed (click UNSET)
# see here:
# http://stackoverflow.com/questions/17643509/conflict-between-dynamic-linking-priority-in-osx
# Although the answer at this link does not talk about DYLD_FRAMEWORK_PATH, it is necessary to
# change it here as advised above.

# arrayfire (not used at the moment)
    LIBS += -L/usr/local/lib -lafopencl

# Setup Qt so Clang works with C++11
    QMAKE_CXXFLAGS += -stdlib=libc++
    LIBS += -stdlib=libc++
    QMAKE_CXXFLAGS += -std=c++11
}

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
    utilities.cpp \
    rgraphicsscene.cpp \
    werner/limb.cpp \
    rscrollarea.cpp

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
    rsubwindow.h \
    utilities.h \
    rsubwindow.h \
    data.h \
    typedefs.h \
    circle.h \
    rgraphicsscene.h \
    werner/limb.h \
    werner/circle.h \
    werner/mystuff.h \
    rscrollarea.h


FORMS    += rmainwindow.ui \
    rscrollarea.ui

RESOURCES += \
    ressources.qrc

DISTFILES += \
    shaders/fragment16uc1.frag \
    shaders/fragment8.frag \
    Icons/icons2/AI/Outlined_Icons_AI.ai \
    Icons/icons3/168 Stroke Icons.ai \
    background/background.pxm \
    Icons/icons2/Icon Font/fonts/outlined-iconset.eot \
    Icons/icons2/Icon Font/fonts/outlined-iconset.woff \
    Icons/myIcons/extract_ROI.pxm \
    Icons/myIcons/extract_ROI_1.pxm \
    Icons/myIcons/extract_ROI_2.pxm \
    Icons/myIcons/selection_1.pxm \
    Icons/myIcons/selection_2.pxm \
    Windows/windown-x64-dependenceis-except ArrayFire.7z \
    Icons/icons2/Icon Font/fonts/outlined-iconset.ttf \
    background/background.jpg \
    screenshots/screenshot_Solar_M74.jpg \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (1).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (10).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (100).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (101).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (102).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (103).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (104).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (105).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (106).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (107).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (108).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (109).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (11).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (110).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (111).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (112).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (113).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (114).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (115).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (116).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (117).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (118).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (119).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (12).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (120).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (121).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (122).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (123).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (124).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (125).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (126).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (127).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (128).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (129).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (13).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (130).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (131).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (132).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (133).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (134).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (135).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (136).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (137).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (138).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (139).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (14).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (140).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (141).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (142).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (143).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (144).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (145).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (146).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (147).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (148).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (149).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (15).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (150).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (151).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (152).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (153).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (154).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (155).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (156).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (157).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (158).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (159).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (16).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (160).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (161).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (162).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (163).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (164).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (165).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (166).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (167).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (168).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (17).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (18).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (19).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (2).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (20).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (21).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (22).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (23).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (24).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (25).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (26).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (27).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (28).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (29).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (3).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (30).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (31).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (32).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (33).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (34).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (35).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (36).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (37).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (38).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (39).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (4).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (40).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (41).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (42).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (43).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (44).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (45).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (46).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (47).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (48).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (49).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (5).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (50).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (51).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (52).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (53).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (54).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (55).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (56).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (57).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (58).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (59).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (6).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (60).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (61).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (62).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (63).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (64).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (65).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (66).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (67).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (68).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (69).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (7).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (70).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (71).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (72).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (73).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (74).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (75).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (76).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (77).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (78).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (79).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (8).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (80).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (81).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (82).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (83).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (84).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (85).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (86).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (87).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (88).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (89).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (9).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (90).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (91).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (92).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (93).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (94).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (95).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (96).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (97).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (98).png \
    Icons/icons3/PNG/Grey/128px/Stroke icons by Dreamstale (99).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (1).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (10).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (100).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (101).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (102).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (103).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (104).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (105).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (106).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (107).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (108).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (109).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (11).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (110).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (111).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (112).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (113).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (114).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (115).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (116).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (117).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (118).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (119).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (12).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (120).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (121).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (122).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (123).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (124).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (125).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (126).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (127).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (128).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (129).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (13).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (130).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (131).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (132).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (133).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (134).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (135).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (136).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (137).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (138).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (139).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (14).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (140).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (141).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (142).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (143).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (144).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (145).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (146).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (147).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (148).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (149).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (15).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (150).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (151).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (152).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (153).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (154).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (155).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (156).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (157).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (158).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (159).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (16).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (160).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (161).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (162).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (163).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (164).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (165).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (166).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (167).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (168).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (17).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (18).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (19).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (2).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (20).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (21).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (22).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (23).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (24).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (25).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (26).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (27).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (28).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (29).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (3).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (30).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (31).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (32).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (33).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (34).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (35).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (36).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (37).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (38).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (39).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (4).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (40).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (41).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (42).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (43).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (44).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (45).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (46).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (47).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (48).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (49).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (5).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (50).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (51).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (52).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (53).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (54).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (55).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (56).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (57).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (58).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (59).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (6).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (60).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (61).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (62).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (63).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (64).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (65).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (66).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (67).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (68).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (69).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (7).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (70).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (71).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (72).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (73).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (74).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (75).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (76).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (77).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (78).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (79).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (8).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (80).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (81).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (82).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (83).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (84).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (85).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (86).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (87).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (88).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (89).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (9).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (90).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (91).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (92).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (93).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (94).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (95).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (96).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (97).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (98).png \
    Icons/icons3/PNG/Grey/64px/Stroke icons by Dreamstale (99).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (1).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (10).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (100).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (101).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (102).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (103).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (104).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (105).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (106).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (107).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (108).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (109).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (11).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (110).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (111).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (112).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (113).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (114).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (115).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (116).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (117).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (118).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (119).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (12).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (120).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (121).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (122).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (123).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (124).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (125).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (126).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (127).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (128).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (129).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (13).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (130).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (131).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (132).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (133).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (134).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (135).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (136).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (137).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (138).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (139).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (14).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (140).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (141).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (142).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (143).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (144).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (145).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (146).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (147).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (148).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (149).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (15).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (150).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (151).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (152).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (153).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (154).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (155).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (156).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (157).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (158).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (159).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (16).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (160).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (161).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (162).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (163).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (164).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (165).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (166).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (167).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (168).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (17).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (18).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (19).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (2).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (20).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (21).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (22).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (23).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (24).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (25).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (26).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (27).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (28).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (29).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (3).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (30).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (31).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (32).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (33).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (34).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (35).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (36).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (37).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (38).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (39).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (4).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (40).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (41).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (42).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (43).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (44).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (45).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (46).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (47).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (48).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (49).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (5).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (50).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (51).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (52).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (53).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (54).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (55).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (56).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (57).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (58).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (59).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (6).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (60).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (61).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (62).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (63).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (64).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (65).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (66).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (67).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (68).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (69).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (7).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (70).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (71).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (72).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (73).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (74).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (75).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (76).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (77).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (78).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (79).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (8).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (80).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (81).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (82).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (83).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (84).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (85).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (86).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (87).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (88).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (89).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (9).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (90).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (91).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (92).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (93).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (94).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (95).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (96).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (97).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (98).png \
    Icons/icons3/PNG/White/128px/Stroke icons by Dreamstale (99).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (1).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (10).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (100).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (101).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (102).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (103).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (104).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (105).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (106).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (107).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (108).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (109).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (11).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (110).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (111).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (112).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (113).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (114).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (115).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (116).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (117).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (118).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (119).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (12).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (120).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (121).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (122).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (123).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (124).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (125).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (126).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (127).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (128).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (129).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (13).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (130).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (131).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (132).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (133).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (134).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (135).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (136).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (137).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (138).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (139).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (14).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (140).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (141).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (142).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (143).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (144).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (145).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (146).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (147).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (148).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (149).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (15).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (150).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (151).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (152).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (153).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (154).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (155).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (156).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (157).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (158).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (159).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (16).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (160).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (161).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (162).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (163).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (164).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (165).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (166).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (167).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (168).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (17).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (18).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (19).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (2).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (20).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (21).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (22).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (23).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (24).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (25).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (26).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (27).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (28).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (29).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (3).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (30).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (31).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (32).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (33).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (34).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (35).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (36).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (37).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (38).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (39).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (4).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (40).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (41).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (42).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (43).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (44).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (45).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (46).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (47).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (48).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (49).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (5).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (50).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (51).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (52).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (53).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (54).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (55).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (56).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (57).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (58).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (59).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (6).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (60).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (61).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (62).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (63).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (64).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (65).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (66).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (67).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (68).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (69).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (7).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (70).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (71).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (72).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (73).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (74).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (75).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (76).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (77).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (78).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (79).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (8).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (80).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (81).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (82).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (83).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (84).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (85).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (86).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (87).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (88).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (89).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (9).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (90).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (91).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (92).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (93).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (94).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (95).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (96).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (97).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (98).png \
    Icons/icons3/PNG/White/48px/Stroke icons by Dreamstale (99).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (1).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (10).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (100).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (101).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (102).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (103).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (104).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (105).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (106).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (107).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (108).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (109).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (11).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (110).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (111).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (112).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (113).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (114).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (115).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (116).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (117).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (118).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (119).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (12).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (120).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (121).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (122).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (123).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (124).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (125).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (126).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (127).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (128).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (129).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (13).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (130).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (131).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (132).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (133).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (134).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (135).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (136).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (137).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (138).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (139).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (14).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (140).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (141).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (142).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (143).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (144).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (145).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (146).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (147).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (148).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (149).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (15).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (150).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (151).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (152).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (153).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (154).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (155).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (156).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (157).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (158).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (159).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (16).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (160).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (161).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (162).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (163).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (164).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (165).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (166).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (167).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (168).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (17).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (18).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (19).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (2).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (20).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (21).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (22).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (23).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (24).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (25).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (26).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (27).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (28).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (29).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (3).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (30).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (31).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (32).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (33).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (34).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (35).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (36).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (37).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (38).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (39).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (4).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (40).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (41).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (42).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (43).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (44).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (45).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (46).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (47).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (48).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (49).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (5).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (50).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (51).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (52).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (53).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (54).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (55).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (56).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (57).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (58).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (59).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (6).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (60).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (61).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (62).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (63).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (64).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (65).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (66).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (67).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (68).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (69).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (7).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (70).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (71).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (72).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (73).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (74).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (75).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (76).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (77).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (78).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (79).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (8).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (80).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (81).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (82).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (83).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (84).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (85).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (86).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (87).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (88).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (89).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (9).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (90).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (91).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (92).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (93).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (94).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (95).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (96).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (97).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (98).png \
    Icons/icons3/PNG/White/64px/Stroke icons by Dreamstale (99).png \
    Icons/icons4/PNG/png-32px/adjustments_32px.png \
    Icons/icons4/PNG/png-32px/alarmclock_32px.png \
    Icons/icons4/PNG/png-32px/anchor_32px.png \
    Icons/icons4/PNG/png-32px/aperture_32px.png \
    Icons/icons4/PNG/png-32px/bargraph_32px.png \
    Icons/icons4/PNG/png-32px/basket_32px.png \
    Icons/icons4/PNG/png-32px/beaker_32px.png \
    Icons/icons4/PNG/png-32px/bike_32px.png \
    Icons/icons4/PNG/png-32px/book-open_32px.png \
    Icons/icons4/PNG/png-32px/briefcase_32px.png \
    Icons/icons4/PNG/png-32px/browser_32px.png \
    Icons/icons4/PNG/png-32px/calendar_32px.png \
    Icons/icons4/PNG/png-32px/camera_32px.png \
    Icons/icons4/PNG/png-32px/caution_32px.png \
    Icons/icons4/PNG/png-32px/chat_32px.png \
    Icons/icons4/PNG/png-32px/circle-compass_32px.png \
    Icons/icons4/PNG/png-32px/clipboard_32px.png \
    Icons/icons4/PNG/png-32px/clock_32px.png \
    Icons/icons4/PNG/png-32px/cloud_32px.png \
    Icons/icons4/PNG/png-32px/compass_32px.png \
    Icons/icons4/PNG/png-32px/desktop_32px.png \
    Icons/icons4/PNG/png-32px/dial_32px.png \
    Icons/icons4/PNG/png-32px/document_32px.png \
    Icons/icons4/PNG/png-32px/documents_32px.png \
    Icons/icons4/PNG/png-32px/download_32px.png \
    Icons/icons4/PNG/png-32px/dribbble_32px.png \
    Icons/icons4/PNG/png-32px/edit_32px.png \
    Icons/icons4/PNG/png-32px/envelope_32px.png \
    Icons/icons4/PNG/png-32px/expand_32px.png \
    Icons/icons4/PNG/png-32px/facebook_32px.png \
    Icons/icons4/PNG/png-32px/flag_32px.png \
    Icons/icons4/PNG/png-32px/focus_32px.png \
    Icons/icons4/PNG/png-32px/gears_32px.png \
    Icons/icons4/PNG/png-32px/genius_32px.png \
    Icons/icons4/PNG/png-32px/gift_32px.png \
    Icons/icons4/PNG/png-32px/global_32px.png \
    Icons/icons4/PNG/png-32px/globe_32px.png \
    Icons/icons4/PNG/png-32px/googleplus_32px.png \
    Icons/icons4/PNG/png-32px/grid_32px.png \
    Icons/icons4/PNG/png-32px/happy_32px.png \
    Icons/icons4/PNG/png-32px/hazardous_32px.png \
    Icons/icons4/PNG/png-32px/heart_32px.png \
    Icons/icons4/PNG/png-32px/hotairballoon_32px.png \
    Icons/icons4/PNG/png-32px/hourglass_32px.png \
    Icons/icons4/PNG/png-32px/key_32px.png \
    Icons/icons4/PNG/png-32px/laptop_32px.png \
    Icons/icons4/PNG/png-32px/layers_32px.png \
    Icons/icons4/PNG/png-32px/lifesaver_32px.png \
    Icons/icons4/PNG/png-32px/lightbulb_32px.png \
    Icons/icons4/PNG/png-32px/linegraph_32px.png \
    Icons/icons4/PNG/png-32px/link_32px.png \
    Icons/icons4/PNG/png-32px/linkedin_32px.png \
    Icons/icons4/PNG/png-32px/lock_32px.png \
    Icons/icons4/PNG/png-32px/magnifying-glass_32px.png \
    Icons/icons4/PNG/png-32px/map-pin_32px.png \
    Icons/icons4/PNG/png-32px/map_32px.png \
    Icons/icons4/PNG/png-32px/megaphone_32px.png \
    Icons/icons4/PNG/png-32px/mic_32px.png \
    Icons/icons4/PNG/png-32px/mobile_32px.png \
    Icons/icons4/PNG/png-32px/newspaper_32px.png \
    Icons/icons4/PNG/png-32px/notebook_32px.png \
    Icons/icons4/PNG/png-32px/paintbrush_32px.png \
    Icons/icons4/PNG/png-32px/paperclip_32px.png \
    Icons/icons4/PNG/png-32px/pencil_32px.png \
    Icons/icons4/PNG/png-32px/phone_32px.png \
    Icons/icons4/PNG/png-32px/picture_32px.png \
    Icons/icons4/PNG/png-32px/pictures_32px.png \
    Icons/icons4/PNG/png-32px/piechart_32px.png \
    Icons/icons4/PNG/png-32px/presentation_32px.png \
    Icons/icons4/PNG/png-32px/pricetags_32px.png \
    Icons/icons4/PNG/png-32px/printer_32px.png \
    Icons/icons4/PNG/png-32px/profile-female_32px.png \
    Icons/icons4/PNG/png-32px/profile-male_32px.png \
    Icons/icons4/PNG/png-32px/puzzle_32px.png \
    Icons/icons4/PNG/png-32px/quote_32px.png \
    Icons/icons4/PNG/png-32px/recycle_32px.png \
    Icons/icons4/PNG/png-32px/refresh_32px.png \
    Icons/icons4/PNG/png-32px/ribbon_32px.png \
    Icons/icons4/PNG/png-32px/rss_32px.png \
    Icons/icons4/PNG/png-32px/sad_32px.png \
    Icons/icons4/PNG/png-32px/scissors_32px.png \
    Icons/icons4/PNG/png-32px/scope_32px.png \
    Icons/icons4/PNG/png-32px/search_32px.png \
    Icons/icons4/PNG/png-32px/shield_32px.png \
    Icons/icons4/PNG/png-32px/speedometer_32px.png \
    Icons/icons4/PNG/png-32px/strategy_32px.png \
    Icons/icons4/PNG/png-32px/streetsign_32px.png \
    Icons/icons4/PNG/png-32px/tablet_32px.png \
    Icons/icons4/PNG/png-32px/target_32px.png \
    Icons/icons4/PNG/png-32px/telescope_32px.png \
    Icons/icons4/PNG/png-32px/toolbox_32px.png \
    Icons/icons4/PNG/png-32px/tools-2_32px.png \
    Icons/icons4/PNG/png-32px/tools_32px.png \
    Icons/icons4/PNG/png-32px/trophy_32px.png \
    Icons/icons4/PNG/png-32px/tumblr_32px.png \
    Icons/icons4/PNG/png-32px/twitter_32px.png \
    Icons/icons4/PNG/png-32px/upload_32px.png \
    Icons/icons4/PNG/png-32px/video_32px.png \
    Icons/icons4/PNG/png-32px/wallet_32px.png \
    Icons/icons4/PNG/png-32px/wine_32px.png \
    Icons/icons4/PNG/png-64px/adjustments_64px.png \
    Icons/icons4/PNG/png-64px/alarmclock_64px.png \
    Icons/icons4/PNG/png-64px/anchor_64px.png \
    Icons/icons4/PNG/png-64px/aperture_64px.png \
    Icons/icons4/PNG/png-64px/bargraph_64px.png \
    Icons/icons4/PNG/png-64px/basket_64px.png \
    Icons/icons4/PNG/png-64px/beaker_64px.png \
    Icons/icons4/PNG/png-64px/bike_64px.png \
    Icons/icons4/PNG/png-64px/book-open_64px.png \
    Icons/icons4/PNG/png-64px/briefcase_64px.png \
    Icons/icons4/PNG/png-64px/browser_64px.png \
    Icons/icons4/PNG/png-64px/calendar_64px.png \
    Icons/icons4/PNG/png-64px/camera_64px.png \
    Icons/icons4/PNG/png-64px/caution_64px.png \
    Icons/icons4/PNG/png-64px/chat_64px.png \
    Icons/icons4/PNG/png-64px/circle-compass_64px.png \
    Icons/icons4/PNG/png-64px/clipboard_64px.png \
    Icons/icons4/PNG/png-64px/clock_64px.png \
    Icons/icons4/PNG/png-64px/cloud_64px.png \
    Icons/icons4/PNG/png-64px/compass_64px.png \
    Icons/icons4/PNG/png-64px/desktop_64px.png \
    Icons/icons4/PNG/png-64px/dial_64px.png \
    Icons/icons4/PNG/png-64px/document_64px.png \
    Icons/icons4/PNG/png-64px/documents_64px.png \
    Icons/icons4/PNG/png-64px/download_64px.png \
    Icons/icons4/PNG/png-64px/dribbble_64px.png \
    Icons/icons4/PNG/png-64px/edit_64px.png \
    Icons/icons4/PNG/png-64px/envelope_64px.png \
    Icons/icons4/PNG/png-64px/expand_64px.png \
    Icons/icons4/PNG/png-64px/facebook_64px.png \
    Icons/icons4/PNG/png-64px/flag_64px.png \
    Icons/icons4/PNG/png-64px/focus_64px.png \
    Icons/icons4/PNG/png-64px/gears_64px.png \
    Icons/icons4/PNG/png-64px/genius_64px.png \
    Icons/icons4/PNG/png-64px/gift_64px.png \
    Icons/icons4/PNG/png-64px/global_64px.png \
    Icons/icons4/PNG/png-64px/globe_64px.png \
    Icons/icons4/PNG/png-64px/googleplus_64px.png \
    Icons/icons4/PNG/png-64px/grid_64px.png \
    Icons/icons4/PNG/png-64px/happy_64px.png \
    Icons/icons4/PNG/png-64px/hazardous_64px.png \
    Icons/icons4/PNG/png-64px/heart_64px.png \
    Icons/icons4/PNG/png-64px/hotairballoon_64px.png \
    Icons/icons4/PNG/png-64px/hourglass_64px.png \
    Icons/icons4/PNG/png-64px/key_64px.png \
    Icons/icons4/PNG/png-64px/laptop_64px.png \
    Icons/icons4/PNG/png-64px/layers_64px.png \
    Icons/icons4/PNG/png-64px/lifesaver_64px.png \
    Icons/icons4/PNG/png-64px/lightbulb_64px.png \
    Icons/icons4/PNG/png-64px/linegraph_64px.png \
    Icons/icons4/PNG/png-64px/link_64px.png \
    Icons/icons4/PNG/png-64px/linkedin_64px.png \
    Icons/icons4/PNG/png-64px/lock_64px.png \
    Icons/icons4/PNG/png-64px/magnifying-glass_64px.png \
    Icons/icons4/PNG/png-64px/map-pin_64px.png \
    Icons/icons4/PNG/png-64px/map_64px.png \
    Icons/icons4/PNG/png-64px/megaphone_64px.png \
    Icons/icons4/PNG/png-64px/mic_64px.png \
    Icons/icons4/PNG/png-64px/mobile_64px.png \
    Icons/icons4/PNG/png-64px/newspaper_64px.png \
    Icons/icons4/PNG/png-64px/notebook_64px.png \
    Icons/icons4/PNG/png-64px/paintbrush_64px.png \
    Icons/icons4/PNG/png-64px/paperclip_64px.png \
    Icons/icons4/PNG/png-64px/pencil_64px.png \
    Icons/icons4/PNG/png-64px/phone_64px.png \
    Icons/icons4/PNG/png-64px/picture_64px.png \
    Icons/icons4/PNG/png-64px/pictures_64px.png \
    Icons/icons4/PNG/png-64px/piechart_64px.png \
    Icons/icons4/PNG/png-64px/presentation_64px.png \
    Icons/icons4/PNG/png-64px/pricetags_64px.png \
    Icons/icons4/PNG/png-64px/printer_64px.png \
    Icons/icons4/PNG/png-64px/profile-female_64px.png \
    Icons/icons4/PNG/png-64px/profile-male_64px.png \
    Icons/icons4/PNG/png-64px/puzzle_64px.png \
    Icons/icons4/PNG/png-64px/quote_64px.png \
    Icons/icons4/PNG/png-64px/recycle_64px.png \
    Icons/icons4/PNG/png-64px/refresh_64px.png \
    Icons/icons4/PNG/png-64px/ribbon_64px.png \
    Icons/icons4/PNG/png-64px/rss_64px.png \
    Icons/icons4/PNG/png-64px/sad_64px.png \
    Icons/icons4/PNG/png-64px/scissors_64px.png \
    Icons/icons4/PNG/png-64px/scope_64px.png \
    Icons/icons4/PNG/png-64px/search_64px.png \
    Icons/icons4/PNG/png-64px/shield_64px.png \
    Icons/icons4/PNG/png-64px/speedometer_64px.png \
    Icons/icons4/PNG/png-64px/strategy_64px.png \
    Icons/icons4/PNG/png-64px/streetsign_64px.png \
    Icons/icons4/PNG/png-64px/tablet_64px.png \
    Icons/icons4/PNG/png-64px/target_64px.png \
    Icons/icons4/PNG/png-64px/telescope_64px.png \
    Icons/icons4/PNG/png-64px/toolbox_64px.png \
    Icons/icons4/PNG/png-64px/tools-2_64px.png \
    Icons/icons4/PNG/png-64px/tools_64px.png \
    Icons/icons4/PNG/png-64px/trophy_64px.png \
    Icons/icons4/PNG/png-64px/tumblr_64px.png \
    Icons/icons4/PNG/png-64px/twitter_64px.png \
    Icons/icons4/PNG/png-64px/upload_64px.png \
    Icons/icons4/PNG/png-64px/video_64px.png \
    Icons/icons4/PNG/png-64px/wallet_64px.png \
    Icons/icons4/PNG/png-64px/wine_64px.png \
    Icons/myIcons/extract_ROI.png \
    screenshots/windows-screen1.png \
    Icons/icons2/Icon Font/fonts/outlined-iconset.svg \
    Icons/icons2/PSD/Outlined_Icons_PSD.psd \
    Icons/icons3/168 Stroke Icons - Grey.eps \
    Icons/icons3/168 Stroke Icons - White.eps \
    Icons/icons2/Icon Font/styles.css \
    Icons/icons2/Icon Font/icons-reference.html \
    Icons/icons3/License.txt \
    fragment32.frag \
    vertex32.vert \
    README.md
