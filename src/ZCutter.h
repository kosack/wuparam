#ifndef ZCUTTER_H
#define ZCUTTER_H

#include "SuperCutter.h"
#include "ImageAnalyzer.h"
class Camera;

/**
 * Originial zenith cutter, from ZCuts99.  Uses SuperCuts-like cuts,
 * but with zenith angle corrections based on cos(theta).
 */
class ZCutter : public SuperCutter {

 public:

    ZCutter();

    void applyCorrections( HillasParameterization &p );
    void setCamera(int npix, int utdate);
    
 protected:

    

 private:

    Camera *_cam;
					   

};

#endif
