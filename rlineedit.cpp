#include "rlineedit.h"

RLineEdit::RLineEdit(QWidget *parent) : QLineEdit(parent)
{



}

void RLineEdit::dragEnterEvent(QDragEnterEvent *event)
{
    event->acceptProposedAction();
}


void RLineEdit::dragMoveEvent(QDragMoveEvent *event)
{
    event->acceptProposedAction();
}

void RLineEdit::dropEvent(QDropEvent *event)
{
    event->acceptProposedAction();
    url = event->mimeData()->urls().at(0);
    this->setText(url.toLocalFile());
}
