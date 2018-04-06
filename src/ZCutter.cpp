#include <iostream>
#include <stdexcept>

#include "Camera.h"
#include "ZCutter.h"


using namespace std;


ZCutter::ZCutter() : SuperCutter() {

    _cam = new Camera(490);

}


/**
 * ZCutter defaults to loading the 490 camera, but if you want
 * something different call this function to load another.
 */
void
ZCutter::setCamera( int npix, int utdate ) {

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
ZCutter::
applyCorrections(HillasParameterization &param) {

    const int NOMINAL_SIZE = 1100;
    double zwidth2, zlength2;
    double sig2_w0, sig2_l0;
    double sig_w0, sig_l0 ;
    double zwidth, zlength;
    double dx,dy, cfact;

    double zen; 
    double cos2_zen; 
    double sig2_psf;
    double sig2_pix;
    double dlog;

    HillasParameterization oldparam = param; //for debugging
    
    zen = param.zenith;
    
//     if(param.event_number < 10) 
// 	cout << "DEBUG: zenith angle="<< zen*180.0/M_PI << endl;

    cos2_zen = cos(zen) ; cos2_zen *= cos2_zen;
    sig2_psf = _cam->getSigmaPSF( zen ) ; sig2_psf *= sig2_psf;
    sig2_pix = _cam->getSigmaPixel() ; sig2_pix *= sig2_pix;

//     cout << "DEBUG: zen = "<< zen << endl;
//     cout << "DEBUG: sigmapsf2 = "<< sig2_psf << endl;
//     cout << "DEBUG: sigmapix2 = "<< sig2_pix << endl;

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
    // Calculate zwidth 
    //============================================================

    sig2_w0 = (param.width * param.width) - (sig2_pix + sig2_psf);    
    dlog = log(param.size/0.449) - 8.0; // the 8.0 is hand fit

    if (sig2_w0 > 0.0) {
	sig2_w0  /= pow(cos(zen),1.5) ;
	
	sig_w0 = sqrt(sig2_w0) - 0.023 * dlog ;
	sig2_w0 = sig_w0 * sig_w0 ;
// 	cout << "DEBUG: SATURATED" << endl;
    }
    else sig2_w0 = 0.0;
 
    zwidth2 = sig2_w0 + sig2_pix + sig2_psf;
    zwidth = sqrt(zwidth2);
    

    //============================================================
    // Calculate zlength 
    //============================================================

//     sig2_l0 = (param.length*param.length - (sig2_pix + sig2_psf));    
//     if (sig2_l0 > 0.0) sig2_l0  /= cos2_zen ;
//     else sig2_l0 = 0.0;

//     sig_l0 = sqrt(sig2_l0) - 0.020 * dlog ;
//     sig2_l0 = sig_l0 * sig_l0 ;

//     zlength2 = sig2_l0 + sig2_pix + sig2_psf;
//     zlength = sqrt(zlength2);

    sig2_l0 = (param.length*param.length - (sig2_pix + sig2_psf));    
    if (sig2_l0 > 0.0) {
	sig2_l0  /= pow(cos(zen),1.2) ;
	sig_l0 = sqrt(sig2_l0) - 0.020 * dlog ;
	sig2_l0 = sig_l0 * sig_l0 ;
    }
    else sig2_l0 = 0.0;


    zlength2 = sig2_l0 + sig2_pix + sig2_psf;
    zlength = sqrt(zlength2);



    //============================================================    
    // Correct the 2D points 
    //============================================================

    // First calculate the distance between the point of origin and
    // centroid

    dx = param.point_of_origin_a.x - param.centroid.x;
    dy = param.point_of_origin_a.y - param.centroid.y;

    // Now calculate the correction factor including two terms: the
    // first ratio corrects the elongation factor, the second replaces
    // the width/length (elongation) dependence with the intrinsic
    // elongation given by zwidth/zlength

    cfact = (_cuts.elongation / 1.68) *  // TODO: shouldn't be hardcoded!
	((1.0-zwidth/zlength)/(1.0-param.width/param.length));
    
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
    
    param.width = zwidth;
    param.length = zlength;
    param.distance /= cos(zen);
    param.miss /= cos(zen);


    // DEBUG: try undoing the transformation and comparing...

//  double unzwidth2 =  pow( ( sqrt(zwidth*zwidth - (sig2_pix+sig2_psf)) 
// 			       + 0.023*(log(param.size/0.449)-8.0)),2 ) 
// 	* pow(cos(zen),1.5) + (sig2_pix+sig2_psf);
    
//     cout << "DEBUG: width="<<oldparam.width
// 	 <<"\t unzwidth="<<sqrt(unzwidth2)
// 	 <<"\t zwidth=" <<zwidth
// 	 << "\t diff="<<oldparam.width-sqrt(unzwidth2)<< endl;


}



