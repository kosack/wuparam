#ifndef SUPERCUTTER_H
#define SUPERCUTTER_H

// void
// example( HillasParameterization &param ) {
//
//     SuperCutter cutter;
//
//     if (cutter.pass( param, ASYMM|ALPHA )) {
// 	    // passed orientation cuts
//
//     }
//
//     if (cutter.pass( param, ALL_CUTS | ~ASYMM ) {
//          // passed all, ignoring asymmetry
//     }
// 
// }

#include <iostream>
#include "Config.h"
#include "ImageAnalyzer.h"

using namespace std;

/**
 * Contains enumerations for SuperCutter and related cut objects
 */
namespace Cuts {
    enum CutMask { SIZE     = 1 << 0, 
		   MAX      = 1 << 1, 
		   FRAC     = 1 << 2,  
		   LENSIZE  = 1 << 3, 
		   DISTANCE = 1 << 4,  
		   LENGTH   = 1 << 5,  
		   WIDTH    = 1 << 6, 
		   ALPHA    = 1 << 7,  
		   ASYMM    = 1 << 8, 
		   ALL_CUTS = 0xFFFFFFFF
    };
};

/**
 * Applies SuperCuts to a given HillasParameterization.  To use, first
 * call setCuts() to supply the cuts to use, then you can check if an
 * event passed by calling the pass() function with appropriate mask
 * (ALL_CUTS will apply all cuts, or you can choose a list of cuts
 * such as SIZE|FRAC|MAX|DISTANCE.
 */
class SuperCutter {

 public:

    void setCuts( CutInfo &cuts );
    virtual void applyCorrections(HillasParameterization &param);
    bool pass( HillasParameterization &param, int cutmask );

    virtual bool passSize( HillasParameterization &param );
    virtual bool passMax( HillasParameterization &param );
    virtual bool passFrac( HillasParameterization &param );
    virtual bool passLensize( HillasParameterization &param );
    virtual bool passLength( HillasParameterization &param );
    virtual bool passWidth( HillasParameterization &param );
    virtual bool passDistance( HillasParameterization &param );
    virtual bool passAlpha( HillasParameterization &param );
    virtual bool passAsymm( HillasParameterization &param );
    
 protected:

    CutInfo _cuts;

};

inline void SuperCutter::setCuts( CutInfo &cuts) {
    _cuts = cuts;
}

inline bool SuperCutter::passSize( HillasParameterization &param ) {
    return (param.size >= _cuts.size.lower && 
	    param.size <= _cuts.size.upper) ;
}

inline bool SuperCutter::passMax( HillasParameterization &param ){
    if (param.max[0] >= _cuts.max1.lower && 
	param.max[0] <= _cuts.max1.upper) {
	if (param.max[1] >= _cuts.max2.lower && 
	    param.max[1] <= _cuts.max2.upper) {
	    if (param.max[2] >= _cuts.max3.lower && 
		param.max[2] <= _cuts.max3.upper) {
		return true;
	    }
	}
    }
    return false;
}

inline bool SuperCutter::passFrac( HillasParameterization &param ){
    return (param.frac[2] >= _cuts.frac3.lower &&
	    param.frac[2] <= _cuts.frac3.upper);
}

inline bool SuperCutter::passLensize( HillasParameterization &param ) {
    return (param.length_over_size >= _cuts.lensize.lower&&
	    param.length_over_size <= _cuts.lensize.upper);
}

inline bool SuperCutter::passLength( HillasParameterization &param ){
    return (param.length >= _cuts.length.lower && 
	    param.length <= _cuts.length.upper);

}

inline bool SuperCutter::passWidth( HillasParameterization &param ){
    return(param.width >= _cuts.width.lower && 
	   param.width <= _cuts.width.upper);
}

inline bool SuperCutter::passDistance( HillasParameterization &param ) {
    return (param.distance >= _cuts.distance.lower && 
	    param.distance <= _cuts.distance.upper);
}

inline bool SuperCutter::passAlpha( HillasParameterization &param ){
    return(param.alpha*180.0/M_PI >= _cuts.alpha.lower && 
	   param.alpha*180/M_PI <= _cuts.alpha.upper);
}

inline bool SuperCutter::passAsymm( HillasParameterization &param ){
    return(param.distance >= _cuts.asymmdist.lower && 
	   param.distance <= _cuts.asymmdist.upper);   
}



/**
 * Responsible for locating files and creating the correct data reader
 * object for the detected file type.  Implemented as a Singleton.
 */
class CutFactory {

 public:
    
    static CutFactory* instance();
    SuperCutter* newCutter( const RunInfo &ri );

 protected:

    CutFactory();
    

 private:

    static CutFactory *pinstance;


};



#endif
