#ifndef RTREEWIDGET_H
#define RTREEWIDGET_H

#include "winsockwrapper.h"
#include <QtCore>
#include <QTreeWidget>
#include <QDragEnterEvent>
#include <QDragLeaveEvent>
#include <QDragMoveEvent>
#include <QDropEvent>
#include <QUrl>

#include <rmat.h>

class RTreeWidget : public QTreeWidget
{
    Q_OBJECT
public:
    RTreeWidget(QWidget *parent = 0);

    QList<RMat*> rMatLightList, rMatBiasList, rMatDarkList, rMatFlatList;

    // getters
    QList<QUrl> getLightUrls() const;
    QList<QUrl> getBiasUrls() const;
    QList<QUrl> getDarkUrls() const;
    QList<QUrl> getFlatUrls() const;
    QDir getLightsDir() const;
    QDir getBiasDir() const;
    QDir getDarkDir() const;
    QDir getFlatDir() const;
    QTreeWidgetItem* getLightItem() const;
    QTreeWidgetItem* getBiasItem() const;
    QTreeWidgetItem* getDarkItem() const;
    QTreeWidgetItem* getFlatItem() const;

    void addItems(QTreeWidgetItem* parentItem, QList<RMat*> rMatImageList);

protected:
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);
    void dropEvent(QDropEvent *event);

signals:


public slots:

    void dispatchUrls(QTreeWidgetItem* parentItem, QList<QUrl> urls);
    void rMatFromLightRButton(QList<RMat*> rMatImageList);
    void rMatFromBiasRButton(QList<RMat*> rMatImageList);
    void rMatFromDarkRButton(QList<RMat*> rMatImageList);
    void rMatFromFlatRButton(QList<RMat*> rMatImageList);


private:

    QTreeWidgetItem *lightItem, *biasItem, *darkItem, *flatItem;
    QList<QUrl> lightUrls, biasUrls, darkUrls, flatUrls;
    QDir lightsDir, biasDir, darkDir, flatDir;

};

#endif // RTREEWIDGET_H
