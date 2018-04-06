#include "ImageAnalyzer.h"
#include "SuperCutter.h"
#include "ZCutter.h"
#include "EZCutter.h"
#include "SpectralCutter.h"
#include "Exceptions.h"


/**
 * Checks the parameterization of an event to see if it passes cuts.
 * Subclasses of this can override the protected methods to
 * re-implement different cut types.  Also, the prePass() and
 * postPass() functions can be re-implemented to allow for special
 * manipulations.  As a rule, no SuperCutter shuold modify the passed
 * in parameters (even though they are passed by reference, which is
 * done for speed).
 *
 * \param param set of hillas parameters to check
 * \param cutmask bitmask of Cuts::CutMask values
 * \returns true if param passwd all cuts specified in cutmast.
 *
 * \todo: right now, prePass occurs every time a cut is
 * specified. That means it gets redone several times for each event!
 * This shuold be done more efficiently without making it more
 * confusing.
 */
bool
SuperCutter::pass( HillasParameterization &param, int cutmask ) {

    bool passed=true;

    // check for invalid event

    if (param.invalid > 0)
	return false;

    // process cutmask:

    if ( cutmask & Cuts::SIZE ) {
	passed &= passSize(param);
    } 
    
    if (cutmask & Cuts::MAX ) {
	passed &= passMax(param);
    }
    
    if (cutmask & Cuts::FRAC ) {
	passed &= passFrac(param);
    }

    if (cutmask & Cuts::LENSIZE ) {
	passed &= passLensize(param);
    }

    if (cutmask & Cuts::DISTANCE ) {
	passed &= passDistance(param);
    }

    if (cutmask & Cuts::LENGTH ) {
	passed &= passLength(param);
    }

    if (cutmask & Cuts::WIDTH ) {
	passed &= passWidth(param);
    }

    if (cutmask & Cuts::ALPHA ) {
	passed &= passAlpha(param);
    }

    if (cutmask & Cuts::ASYMM ) {
	passed &= passAsymm(param);
    }

    return passed;

}



/**
 * Apply a filter to the parameters.  Should always be run before the
 * actual cutting occurs.  Modifies param.
 */
void
SuperCutter::
applyCorrections( HillasParameterization &param ) {

    Coordinate_t oa, ob; // two points of origin
    double elong;   // for elongation correction later...

    elong = _cuts.elongation; // use elong from conf file 

    // should be same as used in ImageAnalyzer:    
    const double NOMINAL_ELONGATION = 1.68; 

    // re-correct the point of origin for the specified elongation
    // factor.  Why? so optimization of cut parameters is faster! It's
    // much harder to optimize elongation when you have to
    // re-parameterize every time - much quicker to do it here.
    
    oa.x = (param.point_of_origin_a.x - param.centroid.x)
	* elong/NOMINAL_ELONGATION + param.centroid.x;
    oa.y = (param.point_of_origin_a.y - param.centroid.y)
	* elong/NOMINAL_ELONGATION + param.centroid.y;
    
    ob.x = (param.point_of_origin_b.x - param.centroid.x)
	* elong/NOMINAL_ELONGATION + param.centroid.x;
    ob.y = (param.point_of_origin_b.y - param.centroid.y)
	* elong/NOMINAL_ELONGATION + param.centroid.y;
    
    param.point_of_origin_a = oa;
    param.point_of_origin_b = ob;
    
}


CutFactory* CutFactory::pinstance = NULL;

/**
 * Returns pointer to the global CutFactory instance
 */
CutFactory* 
CutFactory::
instance() {

    if (pinstance == NULL)
	pinstance = new CutFactory();
    
    return pinstance;

}

CutFactory::CutFactory(){ }

SuperCutter*
CutFactory::newCutter( const RunInfo &ri ) {

    SuperCutter *cutter;
    
    switch (ri.cuttype) {
    case Config::SUPERCUTTER:
	cutter = new SuperCutter();
	break;
    case Config::ZCUTTER:
	cutter = new ZCutter();
	break;
    case Config::SPECTRALCUTTER:
	cutter = new SpectralCutter();
	break;
    case Config::EZCUTTER:
	cutter = new EZCutter();
	break;
    default:
	throw CriticalAnalysisException("Unknown Cutter type specified");
	break;
    }

    return cutter;

}


