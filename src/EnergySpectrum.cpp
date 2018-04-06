#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include "EnergySpectrum.h"
#include "Exceptions.h"

using namespace std;

EnergyEstimatorFactory* EnergyEstimatorFactory::pinstance = NULL;

/**
 * Energy estimator for LZA (Z ~= 30 degrees) data
 * 
 * This function, which was fit to simulation data, returns the
 * log10 of the estimated energy (log10(EstimatedEnergy))
 */
double 
LZAEnergyEstimator::getEstimate(const double size, const double dist, 
				const double zenith ) {
    

    double lnsize = log(size);
    double e1,e2;


    // from fit of log(size) vs log(true energy)
    // third order poly in log(SIZE)
    const double A = 0.34821;
    const double B = -0.32983;
    const double C = 0.084692;
    
    // from fit of distance vs (e1 - true energy)
    // two linear fits of DISTANCE (not log)
    const double DA1 = 0.39281;
    const double DB1 = -0.79686;
    const double DA2 = -3.0901;
    const double DB2 = 3.8822;
    
    // distance at which we should use the second fit line instead of
    // the first
    const double distbreak = 0.75;


    if ( (size>0) && (dist>0)) {

	// The lowest order energy estimator is just a function of log(SIZE)
	e1 = A + B*lnsize + C*lnsize*lnsize;
	
	// The first order correction is a function of DISTANCE
	if (dist <= distbreak) {
	    e2 = DA1 + DB1*dist;
	}
	else {
	    e2 = DA2 + DB2*dist;
	}

	// Finally, the log(energyEstimate) is e1+e2.  To be
	// equivalent to Henric Krawczinski's original function, this
	// needs to return log10(energyEstimate) not the natural log,
	// so we need to multiply by a conversion factor. log10(x) =
	// ln(x)/ln(10)
	
	return (e1+e2)/log(10.0);

    }
    else
	return -100;


};


/**
 * Energy estimator for SZA (Z ~= 60 degrees) data
 * 
 * This function, which was fit to simulation data, returns the
 * log10 of the estimated energy (log10(EstimatedEnergy))
 */
double 
SZAEnergyEstimator::getEstimate(const double size, const double dist,
				const double zenith) {
    
    double lnsize;
    double e1,e2;

    // scale SIZE by zenith angle (fit by eye shows scaling goes like cos^3)
    //double size_correction = pow( cos(21*M_PI/180.0)/cos(zenith), 3);
    double size_correction = 1.0;
    lnsize = log(size*size_correction);

    
    // from fit of log(size) vs log(true energy)
    // third order poly in log(SIZE)
    // logE' = A + Bx + Cx^2

     const double A = -7.0527;
     const double B = 1.2952;
     const double C = -0.034273;
    
    //my z=21deg fits to the MC data with scaled electronics factors
    //const double A = -5.386;
    //const double B = 0.73754;
    //const double C = 5.6625e-3;

    // from fit of distance vs (e1 - log(true energy))
    // two linear fits of DISTANCE (not log)

     const double DA1 = 0.057208;
     const double DB1 = -0.19915;
     const double DA2 = -1.9642;
     const double DB2 = 2.4432;
    
    //my z=21deg fits to the MC data with scaled electronics factors
    //const double DA1 = 0.4433;
    //const double DB1 = -5.846e-2;
    //const double DA2 = -3.6807;
    //const double DB2 = 5.1725;
    
    // distance at which we should use the second fit line instead of
    // the first

     const double distbreak = 0.75;
   
    //my z=21deg fits to the MC data with scaled electronics factors
    //const double distbreak = 0.81;

    if ( (size>0) && (dist>0)) {

	// The lowest order energy estimator is just a function of log(SIZE)
	e1 = A + B*lnsize + C*lnsize*lnsize;
	
	// The first order correction is a function of DISTANCE
	if (dist <= distbreak) {
	    e2 = DA1 + DB1*dist;
	}
	else {
	    e2 = DA2 + DB2*dist;
	}

	// Finally, the log(energyEstimate) is e1+e2.  To be
	// equivalent to Henric Krawczinski's original function, this
	// needs to return log10(energyEstimate) not the natural log,
	// so we need to multiply by a conversion factor. log10(x) =
	// ln(x)/ln(10)
	
	return (e1+e2)/log(10.0);

    }
    else
	return -100;

};


/**
 * Energy estimator for (Z ~= 40 degrees) data
 * 
 * This function, which was fit to simulation data, returns the
 * log10 of the estimated energy (log10(EstimatedEnergy))
 */
double
Z40EnergyEstimator::getEstimate(const double size, const double dist,
				const double zenith) {
    
    double logsize;
    double e1,e2;

    // scale SIZE by zenith angle (fit by eye shows scaling goes like cos^3)
    //double size_correction = pow( cos(21*M_PI/180.0)/cos(zenith), 3);
    double size_correction = 1.0;
    logsize = log(size*size_correction);
    
    // from fit of log(size) vs log(true energy)
    // third order poly in log(SIZE)
    // logE' = A + Bx + Cx^2

    const double A = -4.2;
    const double B = 0.59432;
    const double C = 0.01855;

    
    // from fit of distance vs (e1 - log(true energy))
    // two linear fits of DISTANCE (not log)

    const double DA1 = 0.1238;
    const double DB1 = -0.3477;
    const double DA2 = -4.924;
    const double DB2 = 5.9146;
    
    // distance at which we should use the second fit line instead of
    // the first
    const double distbreak = 0.8;

    if ( (size>0) && (dist>0)) {

	// The lowest order energy estimator is just a function of log(SIZE)
	e1 = A + B*logsize + C*logsize*logsize;
	
	// The first order correction is a function of DISTANCE
	if (dist <= distbreak) {
	    e2 = DA1 + DB1*dist;
	}
	else {
	    e2 = DA2 + DB2*dist;
	}

	// Finally, the log(energyEstimate) is e1+e2.  To be
	// equivalent to Henric Krawczinski's original function, this
	// needs to return log10(energyEstimate) not the natural log:
	
	return (e1+e2);

    }
    else
	return -100;
    
}


/**
 * Energy estimator for SZA (Z ~= 50 degrees) data
 * 
 * This function, which was fit to simulation data, returns the
 * log10 of the estimated energy (log10(EstimatedEnergy))
 */
double 
Z50EnergyEstimator::getEstimate(const double size, const double dist,
				const double zenith) {
    
    double logsize;
    double e1,e2;
    
    // scale SIZE by zenith angle (fit by eye shows scaling goes like cos^3)
    //double size_correction = pow( cos(21*M_PI/180.0)/cos(zenith), 3);
    double size_correction = 1.0;
    logsize = log(size*size_correction);
    
    // from fit of log(size) vs log(true energy)
    // third order poly in log(SIZE)
    // logE' = A + Bx + Cx^2

    const double A = -4.00526;
    const double B = 0.701737;
    const double C = 0.0103804;

    
    // from fit of distance vs (e1 - log(true energy))
    // two linear fits of DISTANCE (not log)

    const double DA1 = 0.263418;
    const double DB1 = -0.537669;
    const double DA2 = -4.64593;
    const double DB2 = 5.66858;
    
    // distance at which we should use the second fit line instead of
    // the first
    const double distbreak = 0.78;

    if ( (size>0) && (dist>0)) {

	// The lowest order energy estimator is just a function of log(SIZE)
	e1 = A + B*logsize + C*logsize*logsize;
	
	// The first order correction is a function of DISTANCE
	if (dist <= distbreak) {
	    e2 = DA1 + DB1*dist;
	}
	else {
	    e2 = DA2 + DB2*dist;
	}

	// Finally, the log(energyEstimate) is e1+e2.  To be
	// equivalent to Henric Krawczinski's original function, this
	// needs to return log10(energyEstimate) not the natural log:
	
	return (e1+e2)/log(10.0);

    }
    else
	return -100;

}




/**
 * Returns pointer to the global CutFactory instance
 */
EnergyEstimatorFactory* 
EnergyEstimatorFactory::
instance() {

    if (pinstance == NULL)
	pinstance = new EnergyEstimatorFactory();
    
    return pinstance;

}

EnergyEstimatorFactory::EnergyEstimatorFactory() { ;}

EnergyEstimator*
EnergyEstimatorFactory::getEstimator( string type ) {
    if (type == "LZA" ) {
	return new LZAEnergyEstimator();
    }
    else if (type == "SZA" ) {
	return new SZAEnergyEstimator();
    }
    else if (type == "Z40") {
	return new Z40EnergyEstimator();
    }
    else if (type == "Z50") {
	return new Z50EnergyEstimator();
    }
    else
	throw CriticalAnalysisException("Unknown EnergyEstimator type: '"
					+type+"'");
}
