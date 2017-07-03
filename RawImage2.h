#ifndef RAWIMAGE2_H
#define RAWIMAGE2_H

// OpenCV
#include <opencv2/world.hpp>
// libraw
#include <libraw.h>
//include <libraw/libraw.h>


class RawImage2
{



public:
    RawImage2();
    cv::Mat getMatCFA() const;

private:
    cv::Mat matCFA;
    LibRaw rawProcess;

};

#endif // RAWIMAGE2_H
