#include "rsubwindow.h"

RSubWindow::RSubWindow(QWidget* parent) :
    QMdiSubWindow(parent)
{

}

void RSubWindow::closeEvent(QCloseEvent *closeEvent)
{
    this->hide();
    emit gotHidden();
    closeEvent->accept();
}

