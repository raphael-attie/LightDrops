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

    QString fileName;
    QList<QUrl> urls;
    urls = event->mimeData()->urls();
    QTreeWidgetItem *parentItem = itemFromIndex(index);
    for (int i = 0 ; i < urls.size() ; i++)
    {
        fileName = urls.at(i).fileName();
        QTreeWidgetItem *treeItem = new QTreeWidgetItem(parentItem);
        treeItem->setText(0, fileName);
    }

    if (parentItem == lightItem)
    {
        lightUrls = urls;
    }
    else if (parentItem == biasItem)
    {
        biasUrls = urls;
    }
    else if (parentItem == darkItem)
    {
        darkUrls = urls;
    }
    else if (parentItem == flatItem)
    {
        flatUrls = urls;
    }

}


QList<QUrl> RTreeWidget::getLightUrls()
{
    return lightUrls;
}

QList<QUrl> RTreeWidget::getBiasUrls()
{
    return biasUrls;
}

QList<QUrl> RTreeWidget::getDarkUrls()
{
    return darkUrls;
}

QList<QUrl> RTreeWidget::getFlatUrls()
{
    return flatUrls;
}
