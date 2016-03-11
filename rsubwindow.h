#ifndef RSUBWINDOW_H
#define RSUBWINDOW_H

#include <QMdiSubWindow>
#include <QCloseEvent>

class RSubWindow : public QMdiSubWindow
{
    Q_OBJECT
public:
    RSubWindow(QWidget* parent = 0);

signals:
    void gotHidden();

protected:
    void closeEvent(QCloseEvent *closeEvent);

};

#endif // RSUBWINDOW_H
