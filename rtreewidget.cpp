#include "rtreewidget.h"
#include <QMimeData>

RTreeWidget::RTreeWidget(QWidget *parent) : QTreeWidget(parent)
{
    lightItem = new QTreeWidgetItem(this, QStringList(tr("Light")));
    lightItem->setFlags(Qt::ItemIsEnabled);
    biasItem = new QTreeWidgetItem(this, QStringList(tr("Bias")));
    biasItem->setFlags(Qt::ItemIsEnabled);
    darkItem = new QTreeWidgetItem(this, QStringList(tr("Dark")));
    darkItem->setFlags(Qt::ItemIsEnabled);
    flatItem = new QTreeWidgetItem(this, QStringList(tr("Flat")));
    flatItem->setFlags(Qt::ItemIsEnabled);


    setFocusPolicy(Qt::ClickFocus);
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

void RTreeWidget::keyPressEvent(QKeyEvent *key)
{
    if ( (key->key() == Qt::Key_Backspace) )
    {
        qDebug("Backspace was pressed");
        if (currentItem()->parent() == biasItem)
        {
            biasUrls.removeAt(0);
            removeItem(currentItem());
        }
        else if (currentItem()->parent() == darkItem)
        {
            darkUrls.removeAt(0);
            removeItem(currentItem());

        }
        else if (currentItem()->parent() == flatItem)
        {
            flatUrls.removeAt(0);
            removeItem(currentItem());
        }
        else
        {
            QTreeWidget::keyPressEvent(key);
            return;
        }

    }
    else
    {
        QTreeWidget::keyPressEvent(key);
    }
}

void RTreeWidget::dispatchUrls(QTreeWidgetItem *parentItem, QList<QUrl> urls)
{

    for (int i = 0 ; i < urls.size() ; i++)
    {
        QTreeWidgetItem *treeItem = new QTreeWidgetItem(parentItem);
        treeItem->setText(0, urls.at(i).fileName());
        qDebug() << "url =" << urls.at(i).toLocalFile();
    }

    QFileInfo filesInfo(urls.at(0).toLocalFile());

    if (parentItem == lightItem)
    {
        qDebug("RTreeWidget:: append to lights");
        lightUrls.append(urls);
        lightsDir = filesInfo.absoluteDir();
    }
    else if (parentItem == biasItem)
    {
        qDebug("RTreeWidget:: append to bias");
        biasUrls.append(urls);
        biasDir = filesInfo.absoluteDir();

    }
    else if (parentItem == darkItem)
    {
        qDebug("RTreeWidget:: append to dark");
        darkUrls.append(urls);
        darkDir = filesInfo.absoluteDir();
    }
    else if (parentItem == flatItem)
    {
        qDebug("RTreeWidget:: append to flat");
        flatUrls.append(urls);
        flatDir = filesInfo.absoluteDir();
        qDebug() << "RTreeWidget:: flatDir (path):" << flatDir.path();
        qDebug() << "RTreeWidget:: flatPath:" << flatDir.filePath(QString("masterFlat.fits"));
    }
}

void RTreeWidget::rMatFromLightRButton(QList<RMat*> rMatImageList)
{
    addItems(lightItem, rMatImageList);
    rMatLightList = rMatImageList;
    lightsDir = rMatImageList.at(0)->getFileInfo().absoluteDir();
}

void RTreeWidget::rMatFromBiasRButton(QList<RMat*> rMatImageList)
{
    addItems(biasItem, rMatImageList);
    rMatBiasList = rMatImageList;
    biasDir = rMatImageList.at(0)->getFileInfo().absoluteDir();

    double min = 0;
    double max = 0;
    cv::minMaxLoc(rMatBiasList.at(0)->matImage, &min, &max);
    qDebug("RTreeWidget::rMatFromBiasRButton:: min =%f , max =%f", min, max );
}

void RTreeWidget::rMatFromDarkRButton(QList<RMat*> rMatImageList)
{
    addItems(darkItem, rMatImageList);
    rMatDarkList = rMatImageList;
    darkDir = rMatImageList.at(0)->getFileInfo().absoluteDir();
}

void RTreeWidget::rMatFromFlatRButton(QList<RMat*> rMatImageList)
{
    addItems(flatItem, rMatImageList);
    rMatFlatList = rMatImageList;
    flatDir = rMatImageList.at(0)->getFileInfo().absoluteDir();
}

void RTreeWidget::removeItem(QTreeWidgetItem *treeItem)
{
    treeItem->parent()->removeChild(treeItem);
}

void RTreeWidget::cleanup(QList<RMat *> rMatImageList)
{
        for (int i=0; i<rMatImageList.size(); i++)
        {
            qDebug("Removing child from parent item");
            QTreeWidgetItem *item = rMatImageList.at(i)->getItem();
            QTreeWidgetItem *parent = item->parent();
            parent->removeChild(rMatImageList.at(i)->getItem());
            if (parent == lightItem)
            {
                // The URL should be unique in that list (if not, that's a user mistake and I need to prevent that.
                // So let's see where the URL of the RMatImageList is (if any) and remove it from the lightUrls.
                QUrl url = rMatImageList.at(i)->getUrl();
                for (int j=0; j<lightUrls.size(); j++)
                {
                    if (url == lightUrls.at(j))
                    {
                        lightUrls.removeAt(j);
                    }
                }

            }
            qDebug("lightUrls.size() = %i", lightUrls.size());
        }
}

void RTreeWidget::addItems(QTreeWidgetItem *parentItem, QList<RMat*> rMatImageList)
{
    QList<QUrl> urls;

    for (int i = 0 ; i < rMatImageList.size() ; i++)
    {
        if (rMatImageList.at(i)->getItem() != NULL)
        {
            lightItem->removeChild(rMatImageList.at(i)->getItem());
            biasItem->removeChild(rMatImageList.at(i)->getItem());
            darkItem->removeChild(rMatImageList.at(i)->getItem());
            flatItem->removeChild(rMatImageList.at(i)->getItem());
        }

        QTreeWidgetItem *treeItem = new QTreeWidgetItem(parentItem);
        treeItem->setText(0, rMatImageList.at(i)->getImageTitle());
        rMatImageList[i]->setItem(treeItem);

        /// Caveat: There may not be a URL. This only make sense if there's one.
        urls.append(rMatImageList.at(i)->getUrl());

    }

    QFileInfo filesInfo(urls.at(0).toLocalFile());

    if (parentItem == lightItem)
    {
        qDebug("RTreeWidget:: append to lights");
        lightUrls.append(urls);
        lightsDir = filesInfo.absoluteDir();
    }
    else if (parentItem == biasItem)
    {
        qDebug("RTreeWidget:: append to bias");
        biasUrls.append(urls);
        biasDir = filesInfo.absoluteDir();

    }
    else if (parentItem == darkItem)
    {
        qDebug("RTreeWidget:: append to dark");
        darkUrls.append(urls);
        darkDir = filesInfo.absoluteDir();
    }
    else if (parentItem == flatItem)
    {
        qDebug("RTreeWidget:: append to flat");
        flatUrls.append(urls);
        flatDir = filesInfo.absoluteDir();
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



