#ifndef RTREEWIDGET_H
#define RTREEWIDGET_H

#include <QtCore>
#include <QTreeWidget>
#include <QDragEnterEvent>
#include <QDragLeaveEvent>
#include <QDragMoveEvent>
#include <QDropEvent>
#include <QUrl>

class RTreeWidget : public QTreeWidget
{
    Q_OBJECT
public:
    RTreeWidget(QWidget *parent = 0);

    // getters
    QList<QUrl> getLightUrls();
    QList<QUrl> getBiasUrls();
    QList<QUrl> getDarkUrls();
    QList<QUrl> getFlatUrls();



protected:
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);
    void dropEvent(QDropEvent *event);

signals:

public slots:

private:

    QTreeWidgetItem *lightItem, *biasItem, *darkItem, *flatItem;
    QList<QUrl> lightUrls, biasUrls, darkUrls, flatUrls;
};

#endif // RTREEWIDGET_H
