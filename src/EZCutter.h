#ifndef EZCUTTER_H
#define EZCUTTER_H

#include "SuperCutter.h"
#include "ImageAnalyzer.h"
class Camera;

/**
 * Extended zenith cutterm Uses SuperCuts-like cuts, but with zenith
 * angle corrections based on cos(theta) and 3rd order poly in
 * log(SIZE).  Good for spectral cuts.
 */
class EZCutter : public SuperCutter {

 public:

    EZCutter();

    void applyCorrections( HillasParameterization &p );
    void setCamera(int npix,int utdate);
    
 protected:

    
 private:

    Camera *_cam;
					   

};

#endif
