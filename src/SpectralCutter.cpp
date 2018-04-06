

#include <iostream>
#include <stdexcept>
#include "SpectralCutter.h"

using namespace std;


#define DEBUG_SPECTRALCUTS 0

bool 
SpectralCutter::
passLength( HillasParameterization &param ) {

   double len = expectedLength( param.size, param.zenith );
   double lower = len - _cuts.dlength;
   double upper = len + _cuts.dlength;

   if (DEBUG_SPECTRALCUTS) {
       cout << "DEBUG: SpectralCuts: length = "<<param.length
	    << " model="<<len
	    << " low="<<lower
	    << " up="<<upper
	    << endl;
   }

   if (lower <= param.length && param.length <= upper ) {
       return true;
   }
   else{
       return false;
   }

}


bool
SpectralCutter:: 
passWidth( HillasParameterization &param ) {

   double wid = expectedWidth( param.size, param.zenith );
   double lower = wid - _cuts.dwidth;
   double upper = wid + _cuts.dwidth;

   if (lower <= param.width && param.width <= upper )
       return true;
   else
       return false;

    
}


/**
 * Returns expected length (from simulations) for a given zenith angle
 * and size
 **/
double
SpectralCutter::
expectedLength( double size, double zenith ) {
    
    const int NZBINS = 7;
    const int NPARAMS = 4;
    
    // ==================================================================
    // Set up the data from fits:
    
    // here are the simulated zenith angles: (note: these must be
    // consecutive integer bins, or you need to change the bracketing
    // calculation below!)
    
    const double zbin[NZBINS] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
    
    // here are the fit parameters that should be passed to modelLength
    
    const double fitparam[NZBINS][NPARAMS] = {
        //  A         B          C           D
        {-1.90459,0.784166,-0.0979834,0.00413616},     //  0
        {-1.49174,0.626015,-0.078051,0.00330978},      // 10
        {-1.36292,0.571065,-0.0703576,0.00294821},     // 20
        {-1.64037,0.638111,-0.0733976,0.00283247},     // 30
        {-1.4705,0.612931,-0.0767032,0.00328061},      // 40
        {-1.15645,0.47459,-0.0573313,0.00238501},      // 50
        {-0.809668,0.404381,-0.0598228,0.00309438 }    // 60
    };
    
    
    // ==================================================================
    // First, check bounds and figure out which two zenith angle bins
    // bracket the given zenith angle:
    
    //     if (zenith >  80*M_PI/180.0) {
    //      throw out_of_range("given zenith angle out of interpolation range");
    //     }
    
    if (size <= 1e-50) {
        throw out_of_range("given size is < 0");
    }
    
    
    // HACK: the zenith angle bins are integers right now, so the
    // fastest way to get the bracketing is to truncate (won't work
    // for non-integer bins though!)
    
    int lowbin = static_cast<int>(floor(zenith/10.0*180/M_PI));
    int highbin = lowbin+1;
    
    if (highbin >= NZBINS ) {
        // just do linear extrapolation if out of range...
        highbin = NZBINS-1;
        lowbin = NZBINS-2;
    }
    
    // ==================================================================
    // now, evaluate modelLength(SIZE) using the bracketing zenith
    // angles and do a linear interpolation between them to get the
    // actual value for the given zenith angle:
    
    double lnsize = log(size);
    double l0 = modelLength( lnsize, fitparam[lowbin] );
    double l1 = modelLength( lnsize, fitparam[highbin] );
    double z0 = zbin[lowbin] *M_PI/180.0;
    double z1 = zbin[highbin]*M_PI/180.0;
    
    double slope = (l1-l0)/(z1-z0);
    double l = slope*(zenith-z0) + l0;
    
    return l;

}



/**
 * Returns expected length (from simulations) for a given zenith angle
 * and size
 **/
double
SpectralCutter::
expectedWidth( double size, double zenith ) {
    
    const int NZBINS = 7;
    const int NPARAMS = 4;
    
    // ==================================================================
    // Set up the data from fits:
    
    // here are the simulated zenith angles: (note: these must be
    // consecutive integer bins, or you need to change the bracketing
    // calculation below!)
    
    const double zbin[NZBINS] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0};
    
    // here are the fit parameters that should be passed to modelWidth
    
    const double fitparam[NZBINS][NPARAMS] = {
        //  A         B          C           D
        {-0.407304,0.178678,-0.0231434,0.00113588},    //  0
        {-0.0869331,0.0363275,-0.0023061,0.000130217}, // 10
        {-0.485086,0.20325,-0.025185,0.00114646} ,     // 20
        {-0.274123,0.115825,-0.0133312,0.000616812},   // 30
        {-0.570203,0.243562,-0.0311984,0.0014157},     // 40
        {-0.189972,0.0792605,-0.00755413,0.000276128}, // 50
        {-0.408642,0.193252,-0.0265448,0.00128923 }    // 60
    };
    
    
    // ==================================================================
    // First, check bounds and figure out which two zenith angle bins
    // bracket the given zenith angle:
    
    //     if (zenith >  80*M_PI/180.0) {
    //      throw out_of_range("given zenith angle out of interpolation range");
    //     }
    
    if (size <= 1e-50) {
        throw out_of_range("given size is < 0");
    }
    
    
    // HACK: the zenith angle bins are integers right now, so the
    // fastest way to get the bracketing is to truncate (won't work
    // for non-integer bins though!)
    
    int lowbin = static_cast<int>(floor(zenith/10.0*180/M_PI));
    int highbin = lowbin+1;
    
    if (highbin >= NZBINS ) {
        // just do linear extrapolation if out of range...
        highbin = NZBINS-1;
        lowbin = NZBINS-2;
    }
    
    // ==================================================================
    // now, evaluate modelLength(SIZE) using the bracketing zenith
    // angles and do a linear interpolation between them to get the
    // actual value for the given zenith angle:
    
    double lnsize = log(size);
    double l0 = modelWidth( lnsize, fitparam[lowbin] );
    double l1 = modelWidth( lnsize, fitparam[highbin] );
    double z0 = zbin[lowbin] *M_PI/180.0;
    double z1 = zbin[highbin]*M_PI/180.0;
    
    double slope = (l1-l0)/(z1-z0);
    double l = slope*(zenith-z0) + l0;
    
    return l;

}


/**
 * Returns value of function:
 * A + B*x + C*x^2 + D*x^3, where x = ln(SIZE)
 */
double
SpectralCutter::
modelLength( double lnsize, const double *fitparam ) {
    
    return (fitparam[A] + fitparam[B]*lnsize +
            fitparam[C]*pow(lnsize,2) + fitparam[D]*pow(lnsize,3));
    
}

/**
 * Returns value of function:
 * A + B*x + C*x^2 + D*x^3, where x = ln(SIZE)
 */
double
SpectralCutter::
modelWidth( double lnsize, const double *fitparam ) {
    
    return (fitparam[A] + fitparam[B]*lnsize +
            fitparam[C]*pow(lnsize,2) + fitparam[D]*pow(lnsize,3));
    
}





