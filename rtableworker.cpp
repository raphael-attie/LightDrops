#include <QApplication>
#include <QHeaderView>

#include "rtableworker.h"


RTableWorker::RTableWorker(RListImageManager *rListImageManager)
{
    frameIndex = 0;
    this->rListImageManager = rListImageManager;
    newFitsImage = rListImageManager->fetchImageManager(frameIndex)->getNewFitsImage();
    nKeys = newFitsImage->getNKeys();
    keyNames = newFitsImage->getKeyNames();
    keyValues = newFitsImage->getKeyValues();
    keyComments = newFitsImage->getKeyComments();

    newTableWidget = new QTableWidget(newFitsImage->getNKeys(), 3);
    newTableWidget->verticalHeader()->setVisible(false);
    newTableWidget->horizontalHeader()->setVisible(true);
    newTableWidget->setMinimumWidth(400);
    newTableWidget->setEditTriggers(QAbstractItemView::NoEditTriggers);

    newTableWidget->setHorizontalHeaderItem(0, new QTableWidgetItem("Keywords"));
    newTableWidget->setHorizontalHeaderItem(1, new QTableWidgetItem("Values"));
    newTableWidget->setHorizontalHeaderItem(2, new QTableWidgetItem("Comments"));

    newTableWidget->horizontalHeaderItem(0)->setTextAlignment(Qt::AlignLeft);
    newTableWidget->horizontalHeaderItem(1)->setTextAlignment(Qt::AlignLeft);
    newTableWidget->horizontalHeaderItem(2)->setTextAlignment(Qt::AlignLeft);

    newTableWidget->setColumnWidth(1, 150);
    newTableWidget->setColumnWidth(2, 400);

    process();

}

RTableWorker::~RTableWorker()
{

}

void RTableWorker::setFrameIndex(int frameIndex)
{
    QMutexLocker locker(&mutex);
    this->frameIndex = frameIndex;
}

void RTableWorker::process()
{

    int rowHeight = 20;
    QFont tFont;
    tFont.setPointSize(12);
    tFont.setFamily("Arial");

    QMutexLocker locker(&mutex);

    newFitsImage = rListImageManager->fetchImageManager(frameIndex)->getNewFitsImage();

    for (int ii=0; ii<nKeys; ii++)
    {
       newTableWidget->setItem(ii, 0,  new QTableWidgetItem(newFitsImage->getKeyNames().at(ii)));
       newTableWidget->setItem(ii, 1,  new QTableWidgetItem(newFitsImage->getKeyValues().at(ii)));
       newTableWidget->setItem(ii, 2,  new QTableWidgetItem(newFitsImage->getKeyComments().at(ii)));
       newTableWidget->item(ii, 0)->setFont(tFont);
       newTableWidget->item(ii, 1)->setFont(tFont);
       newTableWidget->item(ii, 2)->setFont(tFont);
       newTableWidget->setRowHeight(ii, rowHeight);
    }

    newTableWidget->viewport()->update();

    emit finished();
}

QTableWidget* RTableWorker::getNewTableWidget()
{
    return newTableWidget;
}

