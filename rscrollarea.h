#ifndef RSCROLLAREA_H
#define RSCROLLAREA_H

#include <QScrollArea>

namespace Ui {
class RScrollArea;
}

class RScrollArea : public QScrollArea
{
    Q_OBJECT

public:
    explicit RScrollArea(QWidget *parent = 0);
    ~RScrollArea();

private:
    Ui::RScrollArea *ui;
};

#endif // RSCROLLAREA_H
