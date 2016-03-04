#ifndef RLINEEDIT_H
#define RLINEEDIT_H

#include <QtCore>
#include <QLineEdit>
#include <QDragEnterEvent>
#include <QDragLeaveEvent>
#include <QDragMoveEvent>
#include <QDropEvent>
#include <QUrl>

class RLineEdit : public QLineEdit
{
    Q_OBJECT
public:
    explicit RLineEdit(QWidget *parent = 0);

protected:
    void dragEnterEvent(QDragEnterEvent *event);
    void dragMoveEvent(QDragMoveEvent *event);
    void dropEvent(QDropEvent *event);

signals:

public slots:

private:

    QUrl url;
};

#endif // RLINEEDIT_H
