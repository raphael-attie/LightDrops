#include "rtreewidget.h"
#include <QMimeData>

RTreeWidget::RTreeWidget(QWidget *parent) : QTreeWidget(parent)
{
    lightItem = new QTreeWidgetItem(this, QStringList(tr("Light")));
    biasItem = new QTreeWidgetItem(this, QStringList(tr("Bias")));
    darkItem = new QTreeWidgetItem(this, QStringList(tr("Dark")));
    flatItem = new QTreeWidgetItem(this, QStringList(tr("Flat")));

}

void RTreeWidget::dragEnterEvent(QDragEnterEvent *event)
{
    event->acceptProposedAction();
}

void RTreeWidget::dragMoveEvent(QDragMoveEvent *event)
{
    event->acceptProposedAction();
    QModelIndex index = indexAt(event->pos());
    QTreeWidgetItem *treeItem = itemFromIndex(index);
    setCurrentItem(treeItem);
}

void RTreeWidget::dropEvent(QDropEvent *event)
{
    event->acceptProposedAction();

    QModelIndex index = indexAt(event->pos());
    // just in case the selection or mouse hover is off limit
    if (!index.isValid())
    {
        event->setDropAction(Qt::IgnoreAction);
        return;
     }

    QList<QUrl> urls = event->mimeData()->urls();
    QTreeWidgetItem *parentItem = itemFromIndex(index);

    dispatchUrls(parentItem, urls);
}

void RTreeWidget::dispatchUrls(QTreeWidgetItem *parentItem, QList<QUrl> urls)
{

    for (int i = 0 ; i < urls.size() ; i++)
    {
        QTreeWidgetItem *treeItem = new QTreeWidgetItem(parentItem);
        treeItem->setText(0, urls.at(i).fileName());
    }

    QFileInfo filesInfo(urls.at(0).toLocalFile());

    if (parentItem == lightItem)
    {
        lightUrls = urls;
        lightsDir = filesInfo.absoluteDir();
    }
    else if (parentItem == biasItem)
    {
        biasUrls = urls;
        biasDir = filesInfo.absoluteDir();
    }
    else if (parentItem == darkItem)
    {
        darkUrls = urls;
        darkDir = filesInfo.absoluteDir();
    }
    else if (parentItem == flatItem)
    {
        flatUrls = urls;
        flatDir = filesInfo.absoluteDir();
    }
}

void RTreeWidget::rMatFromLightRButton(QList<RMat> *rMatImageList)
{
    addItems(lightItem, rMatImageList);
}

void RTreeWidget::rMatFromBiasRButton(QList<RMat> *rMatImageList)
{
    addItems(biasItem, rMatImageList);
}

void RTreeWidget::rMatFromDarkRButton(QList<RMat> *rMatImageList)
{
    addItems(darkItem, rMatImageList);
}

void RTreeWidget::rMatFromFlatRButton(QList<RMat> *rMatImageList)
{
    addItems(flatItem, rMatImageList);
}

void RTreeWidget::addItems(QTreeWidgetItem *parentItem, QList<RMat> *rMatImageList)
{
    cleanUpItems(lightItem, rMatImageList);
    cleanUpItems(biasItem, rMatImageList);
    cleanUpItems(darkItem, rMatImageList);
    cleanUpItems(flatItem, rMatImageList);

    for (int i = 0 ; i < rMatImageList->size() ; i++)
    {
        QTreeWidgetItem *treeItem = new QTreeWidgetItem(parentItem);
        treeItem->setText(0, rMatImageList->at(i).getImageTitle());
        rMatImageList->operator [](i).setItem(treeItem);
//        qDebug("RTreeWidget::addItems:: setting items in rMatImageList");
//        qDebug() << "RTreeWidget::addItems:: rMatImageList->value(i)->getItem() =" << rMatImageList->at(i).getItem();
    }
}

void RTreeWidget::cleanUpItems(QTreeWidgetItem *parentItem, QList<RMat> *rMatImageList)
{
    for (int i = 0 ; i < rMatImageList->size() ; i++)
    {
        //qDebug() << "RTreeWidget::cleanUpItems:: rMatImageList.[i].getItem() =" << rMatImageList[i].getItem();
        if (rMatImageList->at(i).getItem() != NULL)
        {
            parentItem->removeChild(rMatImageList->at(i).getItem());
            //qDebug("RTreeWidget::cleanUpItems:: Cleaned up items.");
        }
    }

}

QList<QUrl> RTreeWidget::getLightUrls() const
{
    return lightUrls;
}

QList<QUrl> RTreeWidget::getBiasUrls() const
{
    return biasUrls;
}

QList<QUrl> RTreeWidget::getDarkUrls() const
{
    return darkUrls;
}

QList<QUrl> RTreeWidget::getFlatUrls() const
{
    return flatUrls;
}

QDir RTreeWidget::getLightsDir() const
{
    return lightsDir;
}

QDir RTreeWidget::getBiasDir() const
{
    return biasDir;
}

QDir RTreeWidget::getDarkDir() const
{
    return darkDir;
}

QDir RTreeWidget::getFlatDir() const
{
    return flatDir;
}

QTreeWidgetItem *RTreeWidget::getLightItem() const
{
    return lightItem;
}

QTreeWidgetItem *RTreeWidget::getBiasItem() const
{
    return biasItem;
}

QTreeWidgetItem *RTreeWidget::getDarkItem() const
{
    return darkItem;
}

QTreeWidgetItem *RTreeWidget::getFlatItem() const
{
    return flatItem;
}



