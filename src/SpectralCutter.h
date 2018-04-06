#ifndef SPECTRALCUTTER_H
#define SPECTRALCUTTER_H

#include "SuperCutter.h"


/**
 * Implements cuts which scale with zenith angle and SIZE, otherwise
 * very similar to SuperCuts. Re-implements the length and width cuts.
 */
class SpectralCutter : public SuperCutter {

 public:

    SpectralCutter() : SuperCutter() {;}

    bool passLength( HillasParameterization &param );
    bool passWidth( HillasParameterization &param );

    double expectedLength( double size, double zenith );
    double modelLength( double lnsize, const double *param );

    double expectedWidth( double size, double zenith );
    double modelWidth( double lnsize, const double *param );

    enum ParameterEnum { A,B,C,D };

 private:



};



#endif
