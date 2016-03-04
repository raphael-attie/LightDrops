#include "scrollarea.h"
#include "ui_scrollarea.h"

ScrollArea::ScrollArea(QWidget *parent) :
    QScrollArea(parent),
    ui(new Ui::ScrollArea)
{
    ui->setupUi(this);
}

ScrollArea::~ScrollArea()
{
    delete ui;
}
