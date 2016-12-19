#ifndef RGRAPHICSSCENE_H
#define RGRAPHICSSCENE_H

#include <QKeyEvent>
#include <QGraphicsScene>
#include <QGraphicsItem>
#include <QGraphicsRectItem>

//opencv
#include <opencv2/world.hpp>

class RGraphicsScene : public QGraphicsScene
{

    Q_OBJECT
public:
    RGraphicsScene();
    RGraphicsScene(QObject* parent);
    bool squareMode;
    int blkSize;

signals:
    void ROIsignal(QRect rect);


protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void mouseMoveEvent(QGraphicsSceneMouseEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void keyPressEvent(QKeyEvent *event);
    QRect convertRect(QRect &rect);

private:


//    enum Mode {NoMode, SelectObject, DrawLine};
//    Mode sceneMode;
    QPointF origPoint;
    QGraphicsRectItem* rectItem;

};

#endif // RGRAPHICSSCENE_H
