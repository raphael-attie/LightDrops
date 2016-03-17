#ifndef RTABLEWORKER
#define RTABLEWORKER

#include "winsockwrapper.h"
#include <QtCore>
#include <QTableWidget>
#include <QDockWidget>

#include "MyFitsImage.h"
#include "rlistimagemanager.h"
#include "imagemanager.h"

class RTableWorker : public QObject
{
    Q_OBJECT

public:

    RTableWorker(RListImageManager *rListImageManager);
    ~RTableWorker();

    void setFrameIndex(int frameIndex);
    QTableWidget *getNewTableWidget();

public slots:
    void process();

signals:
    void finished();

private:
    QMutex mutex;

    MyFitsImage* newFitsImage;
    QTableWidget *newTableWidget;
    RListImageManager *rListImageManager;
    int frameIndex;
    int nKeys;
    QVector<QString> keyNames;
    QVector<QString> keyValues;
    QVector<QString> keyComments;

};

#endif // RTABLEWORKER

