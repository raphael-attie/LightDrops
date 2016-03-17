#ifndef SCROLLAREA_H
#define SCROLLAREA_H

#include "winsockwrapper.h"
#include <QScrollArea>

namespace Ui {
class ScrollArea;
}

class ScrollArea : public QScrollArea
{
    Q_OBJECT

public:
    explicit ScrollArea(QWidget *parent = 0);
    ~ScrollArea();

private:
    Ui::ScrollArea *ui;
};

#endif // SCROLLAREA_H
