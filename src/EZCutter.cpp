#include <iostream>
#include <stdexcept>
#include <iomanip>

#include "Camera.h"
#include "EZCutter.h"

using namespace std;


EZCutter::EZCutter() : SuperCutter() {

    _cam = new Camera(490);

}


/**
 * EZCutter defaults to loading the 490 camera, but if you want
 * something different call this function to load another.
 */
void
EZCutter::setCamera( int npix, int utdate ) {

    delete _cam;
    _cam = new Camera(npix,utdate);

}


/**
 * Perform zenith correction. Overrides the SuperCutter::prePass()
 * function which is executed before cuts are applied by
 * SuperCutter::pass()
 *
 * \todo: make sigma_pix and the "0.023" free optimization parameters!
 */
void
EZCutter::
applyCorrections(HillasParameterization &param) {

    double dx,dy, cfact;
    
    // constants from fits to simulated data. The functional form for
    // width and length is as follows:
    //
    // f(f',x) = [( sqrt(f'^2-A^2) + B(x-C) + E(x-C)^2 + F(x-C)^3)^2 *
    // (cos(60)^COSPOW)/(cos(60)^\gamma) * cos(zenith)^\gamma +
    // A]^(1/2) where x is the log(size) and f' is the zwidth or
    // zlength value.  The function has to be inverted to convert from
    // width -> zwidth.  D is the fitted mean f' value from simulations. 
   
    const double WA = 0.003;
    const double WB = 0.04679;
    const double WC = 9.866;
    const double WD = 0.1832;
    const double WE = 0.01534;
    const double WF = 0.00248;
    const double WGAMMA = 0.949;
    const double WCOSPOW = 1.5;

    const double LA = 0.003;
    const double LB = 0.0792;
    const double LC = 9.867;
    const double LD = 0.3023;
    const double LE = 0.02700;
    const double LF = 0.00374;
    const double LGAMMA = 0.714;
    const double LCOSPOW = 1.2;

    double zen; 

    HillasParameterization oldparam = param; //for debugging
    zen = param.zenith;
    
    //============================================================
    // First, apply camera gain correction (PE/DC)
    //============================================================

    // l/s must be corrected after PE to DC correction, but before
    // zenith angle

    param.max[0] *= _cam->getPEToDC();
    param.max[1] *= _cam->getPEToDC();
    param.max[2] *= _cam->getPEToDC();
    param.size   *= _cam->getPEToDC();
    param.length_over_size = (param.length/(double)(param.size)); 

    //============================================================
    // Calculate ezwidth , ezlength
    //============================================================

    const double cos60 = cos(60.0*M_PI/180.0);
    double ezwidth2, ezlength2;
    double ezwidth, ezlength, wzenithfactor,lzenithfactor;
    double x,wshift2,lshift2,lterm1,wterm1;
    
    wzenithfactor = ( (pow(cos60, WGAMMA)/pow(cos60,WCOSPOW)) 
			    /pow( cos(zen), WGAMMA ));

    lzenithfactor = ( (pow(cos60, LGAMMA)/pow(cos60,LCOSPOW))
			     /pow( cos(zen), LGAMMA ));
    x=log(param.size/0.4489); // need to divide out .4489 since fit
			      // values assume 490 camera with no
			      // corrections


    wshift2 = param.width*param.width - WA;
    wterm1 = (wshift2*wzenithfactor>1e-20)? 
	sqrt(wshift2*wzenithfactor): 1e-20;

    ezwidth2 = WA + pow( wterm1 - WB*(x-WC) - WE*pow(x-WC,2) 
			    - WF*pow(x-WC,3) ,2 );
    ezwidth = (ezwidth2>0.0)? sqrt(ezwidth2) : 0.0;

    // ---- zlength

    lshift2 = param.length*param.length - LA;
    lterm1 = (lshift2*lzenithfactor>1e-20)? 
	sqrt(lshift2*lzenithfactor) :1e-20;
    
    ezlength2 = LA + pow( lterm1 - LB*(x-LC) - LE*pow(x-LC,2) 
			     - LF*pow(x-LC,3) ,2 );

    ezlength = (ezlength2>0.0)? sqrt(ezlength2) : 0.0;
    
    cout.precision(5);
//     cout << "DEBUG: width: "<<setw(8)<<param.width
// 	 << " size: "<< setw(8) << param.size 
// 	 << " x: "<< setw(8) << x
// 	 << " ezwidth: "<<setw(8)<<ezwidth
//  	 << " lzfact: "<<setw(8)<<lzenithfactor
//  	 << " length: "<<setw(8)<<param.length
//  	 << " ezlength2: "<<setw(8)<<ezlength2
//  	 << " ezlength: "<<setw(8)<<ezlength
//  	 << endl;


    //============================================================    
    // Correct the 2D points 
    //============================================================

    double newelongation;

    // First calculate the distance between the point of origin and
    // centroid

    dx = param.point_of_origin_a.x - param.centroid.x;
    dy = param.point_of_origin_a.y - param.centroid.y;

    // Now calculate the correction factor including two terms: the
    // first ratio corrects the elongation factor, the second replaces
    // the width/length (elongation) dependence with the intrinsic
    // elongation given by zwidth/zlength

    // KPK: the elongation factor is now a function of x=ln(size), fit
    // from simulations.  e(x) = 1.0894 + 0.092611*x.  Here x
    // shouldn't be divided by 0.4489 because the correction factor
    // was already applied in the simulations!
    
    x = log(param.size);
    newelongation = 1.0894 + 0.092611*x;

    // correct the old point of origin based on the new elongation factor:
    // the Parameterizer always assumes a elongation factor of 1.68. 

    cfact = (newelongation / 1.68) *  // TODO: 1.68 shouldn't be hardcoded!
	((1.0-ezwidth/ezlength)/(1.0-param.width/param.length));
    
    // Modified by JB 030201 should just add the rescaled distance
    // from centroid to point of origin to the orignal centroid
    // position for point a, subtract fo point b

    param.point_of_origin_a.x = param.centroid.x + cfact*dx;
    param.point_of_origin_a.y = param.centroid.y + cfact*dy;

    param.point_of_origin_b.x = param.centroid.x - cfact*dx;
    param.point_of_origin_b.y = param.centroid.y - cfact*dy;

    //============================================================
    // set all the values
    //============================================================
    
    // note miss is also scaled by cos(theta) in case somebody wants
    // to re-calculate alpha from miss and distance later.
    
    param.width = ezwidth;
    param.length = ezlength;
//  param.distance /= cos(zen); // from now on, don't correct distance
                                //  (it's not needed)
    param.miss /= cos(zen);


}



