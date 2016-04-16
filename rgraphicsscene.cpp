#include "rgraphicsscene.h"
#include <QDebug>
#include <QGraphicsSceneMouseEvent>
#include <QRectF>

RGraphicsScene::RGraphicsScene() : QGraphicsScene()
{

    rectItem = 0;
}

RGraphicsScene::RGraphicsScene(QObject* parent) : QGraphicsScene(parent)
{
    rectItem = 0;
}

void RGraphicsScene::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
//    qDebug("origin cursor at (%f , %f)", event->scenePos().x(), event->scenePos().y());
    origPoint = event->scenePos();
    QGraphicsScene::mousePressEvent(event);
}

void RGraphicsScene::mouseMoveEvent(QGraphicsSceneMouseEvent *event)
{
    int mouseX = event->scenePos().x();
    int mouseY = event->scenePos().y();

    int width  = abs(mouseX - origPoint.x());
    int height = abs(mouseY -origPoint.y());

    if ( (mouseX == origPoint.x()) ||  (mouseY == origPoint.y()) )
    {
        return;
    }

        /// If the rectangle has odd size, we have to make it even
        /// by reducing the size by 1 px to make it even.

    int moduloWidth = width % 2;
    int moduloHeight = height % 2;

        if ( moduloWidth != 0 && mouseX < origPoint.x())
        {
            // The cursor is slided left from origin
             mouseX += 1;
        }
        else if ( moduloWidth != 0 && mouseX > origPoint.x())
        {
            // The cursor is slided right from origin
            mouseX -= 1;
        }

        if ( moduloHeight != 0 && mouseY < origPoint.y())
        {
            /// The cursor is slided downward from origin,
            /// move it back upward.
             mouseY += 1;
        }
        else if ( moduloHeight != 0 && mouseY > origPoint.y())
        {
            /// The cursor is slided upward from origin,
            /// move it back downward;
            mouseY -= 1;
        }

        QRect rect;

        if (mouseX < origPoint.x())
        {
            rect.setLeft(mouseX);
            rect.setWidth(abs(mouseX - origPoint.x()));
            //rect.setRight(origPoint.x());
        }
        else
        {
            rect.setLeft(origPoint.x());
            rect.setWidth(abs(mouseX - origPoint.x()));
            //rect.setRight(mouseX);
        }

        if (mouseY < origPoint.y())
        {
            rect.setTop(mouseY);
            rect.setHeight(abs(mouseY - origPoint.y()));
            //rect.setBottom(origPoint.y());
        }
        else
        {
            rect.setTop(origPoint.y());
            rect.setHeight(abs(mouseY - origPoint.y()));
            //rect.setBottom(mouseY);
        }

        if(!rectItem)
        {
            rectItem = new QGraphicsRectItem();
            rectItem->setPen(QPen(Qt::green, 2, Qt::SolidLine));
            this->addItem(rectItem);
        }
        rectItem->setRect(rect);

        QRect correctedRect = rect;
        correctedRect.setTop(this->height() - rect.y() - rect.height() - 1);
        correctedRect.setHeight(rect.height());

        emit ROIsignal(correctedRect);
}

void RGraphicsScene::mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
{
    QGraphicsScene::mouseReleaseEvent(event);
}

void RGraphicsScene::keyPressEvent(QKeyEvent *event)
{
    qDebug("key press event");

    if(event->key() == Qt::Key_Backspace || event->key() == Qt::Key_Delete)
    {
        removeItem(rectItem);
        delete rectItem;
        rectItem = 0;
    }
    else
        QGraphicsScene::keyPressEvent(event);
}
