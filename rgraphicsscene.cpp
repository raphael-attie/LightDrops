#include "rgraphicsscene.h"
#include <QDebug>
#include <QGraphicsSceneMouseEvent>
#include <QRectF>

RGraphicsScene::RGraphicsScene() : QGraphicsScene()
{
    rectItem = NULL;
    squareMode = false;
    blkSize = 0;
}

RGraphicsScene::RGraphicsScene(QObject* parent) : QGraphicsScene(parent)
{
    rectItem = NULL;
    squareMode = false;
    blkSize = 0;
}

void RGraphicsScene::mousePressEvent(QGraphicsSceneMouseEvent *event)
{
//    qDebug("origin cursor at (%f , %f)", event->scenePos().x(), event->scenePos().y());
    origPoint = event->scenePos();

    if (squareMode)
    {
        QRect rect;
        rect.setLeft(origPoint.x() - blkSize/2);
        rect.setTop(origPoint.y() - blkSize/2);
        rect.setWidth(blkSize);
        rect.setHeight(blkSize);

        if(!rectItem)
        {
            rectItem = new QGraphicsRectItem();
            rectItem->setPen(QPen(Qt::green, 2, Qt::SolidLine));
            this->addItem(rectItem);
        }
        // Draw the rectangle
        rectItem->setRect(rect);
        // Send it to the caller for extracting the square ROI.
        QRect correctedRect = convertRect(rect);
        emit ROIsignal(correctedRect);
    }

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
        /// If mouse is just at the origin point, do not draw anything
        /// as a rectangle cannot be a point.
        return;
    }

    QRect rect;

    /// If the rectangle has odd size, make it even by reducing the size by 1 px.
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

    /// Below the "top" coordinate of the rectangle is counter-intuitive
    /// This is due to the inverted orientation of the y-axis.
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
    // Draw the rectangle
    rectItem->setRect(rect);
    // Send it to the caller for cropping the data later if needed.
    QRect correctedRect = convertRect(rect);
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

QRect RGraphicsScene::convertRect(QRect &rect)
{
    /// For ROI selection, convert to the proper coordinates system of the data
    /// reference frame
    QRect correctedRect = rect;
    correctedRect.setTop(this->height() - rect.y() - rect.height() - 1);
    correctedRect.setHeight(rect.height());
    return correctedRect;
}
