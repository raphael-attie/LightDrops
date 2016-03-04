#include "rscrollarea.h"
#include "ui_rscrollarea.h"

RScrollArea::RScrollArea(QWidget *parent) :
    QScrollArea(parent),
    ui(new Ui::RScrollArea)
{
    ui->setupUi(this);
}

RScrollArea::~RScrollArea()
{
    delete ui;
}
