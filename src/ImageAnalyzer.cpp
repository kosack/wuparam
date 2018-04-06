// ImageAnalyzer.cpp
// Karl Kosack <kosack@hbar.wustl.edu>


#include <iostream>
#include <iomanip>
#include <valarray>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_math.h>

#include "Camera.h"
#include "Types.h"
#include "ImageAnalyzer.h"
#include "DataWriter.h"

using namespace std;
void ImageAnalyzer::init( Array_t &xcoords, Array_t &ycoords, int n ) {

    _nadc = n;

    int npix = xcoords.size();
    
    _x.resize( npix );
    _y.resize( npix  );
    _x2.resize( npix );
    _y2.resize( npix );
    _xy.resize( npix );

    //=========================================================================
    // Precalculate stuff for speed

    _x = xcoords;
    _y = ycoords;
    _x2  = _x * _x;			// these are valarrays, so we 
    _y2  = _y * _y;			// can just multiply them 
    _xy  = _x * _y;			// directly (no loop needed!)


}


/**
 * Construct an image analyzer
 */
ImageAnalyzer::ImageAnalyzer( Camera *cam )  {

    // KPK: used last pixel instead of getNumPixels to ignore outer tubes

    init(cam->xCoords(), cam->yCoords(),cam->getLastPixel() );
   
    

}

/**
 *  Create a new HillasImageAnalyzer with the camera geometry specified by
 *  cam.  
 *
 * @param cam The Camera object that the HillasImageAnalyzer should use to
 * determine the camera geometry
 */
HillasImageAnalyzer::
HillasImageAnalyzer( Camera *cam, double elongation ) 
    :  _elongation(elongation), ImageAnalyzer(cam)
{
    

    _2d_is_enabled = true;

}


/**
 * Parameterize the image. This is the most important (and probably
 * the most called) function in the analysis, so it should be as
 * optimized as possible!
 *
 * @param image the pixel data for the image you want to parameterize
 * @param cleanpixels is a list of the indices of the picture/boundary tubes
 * parameters should be stored.
 */
void
HillasImageAnalyzer::
parameterize(const Array_t &image, vector<int> &cleanpixels){

    register int i;
    vector<int>::iterator iter;
    double totalsignal=0;
    double momx=0,momy=0, momx2=0, momy2=0, momxy=0, momx3=0, momy3=0;
    double momx2y=0, momxy2=0;
    double vx2, vy2, vxy, vx3, vy3, vx2y, vxy2;
    double d,z, tanpsi_numerator, tanpsi_denominator, sinalpha, azwidth2,z2;
    double cospsi, sinpsi;
    double disp;
    Coordinate_t pa,pb,tempcoord;
    double tmp;
    double u,v;
    double apsi;
    Array_t qx, qy;

    _param.invalid = 0;

    //=========================================================================
    // Get total signal (sum over picture and boundary)

    _param.pixels_in_picture = 0;

    for (iter=cleanpixels.begin(); iter != cleanpixels.end(); iter++) {
	totalsignal += image[*iter];
	_param.pixels_in_picture++;
    }

    _param.size = totalsignal;

    if (totalsignal < EPSILON) {
	clearParam(_param);
	_param.invalid=1;
	return;
    }

    //=========================================================================
    // Calculate moments <x>,<y>,<x^2>,<y^2>,<xy> and use values
    // precalculated in the constructor to minimize multiplies.  Third
    // order moments are calculated later, when the the asymmetry
    // calculation is done (since they myst be centered on the point
    // of origin, not the camera origin)

    for (iter=cleanpixels.begin(); iter != cleanpixels.end(); iter++) {
	momx    += image[*iter] * _x[*iter]  ;
	momy    += image[*iter] * _y[*iter];  
	momx2   += image[*iter] * _x2[*iter];
	momy2   += image[*iter] * _y2[*iter];
	momxy   += image[*iter] * _xy[*iter];
    }

    momx   /= totalsignal;
    momy   /= totalsignal;
    momx2  /= totalsignal;
    momy2  /= totalsignal;
    momxy  /= totalsignal;
    
    //=========================================================================
    // Calculate variances (sigma = sqrt(variance))

    vx2  = (momx2 - momx*momx);	     // <x^2> - <x>^2
    vy2  = (momy2 - momy*momy);	     // <y^2> - <y>^2
    vxy  = (momxy - momx*momy);	     // <xy> - <x><y>

    //=========================================================================
    // Calculate some common factors...

    d = vy2 - vx2; 
    
    z2 = d*d + 4.0*vxy*vxy;
    if (z2>0) 
	z = sqrt(z2) ;
    else {
	z = 0.0;
	_param.invalid++;
    }

    //======================================================================== 
    // Calculate Centroid

    _param.centroid.x = momx;
    _param.centroid.y = momy;

    //=========================================================================
    // Calculate width and length

    tmp  = (vx2 + vy2 - z )/2.0;
    if (tmp>0) _param.width = sqrt(tmp);
    else {_param.width = 0.0; _param.invalid++;}

    tmp = ( vx2 + vy2 + z )/2.0;
    if (tmp>0) _param.length = sqrt(tmp);
    else {_param.length=0.0; _param.invalid++;}

    //=========================================================================
    // Calculate length over size
    
    _param.length_over_size = _param.length / _param.size;

    //=========================================================================
    // Calculate distance

    _param.distance = gsl_hypot(momx, momy);

    //=========================================================================
    // Calculate miss

    if (z > 0.0) {

	u = 1.0+d/z;
	v = 2.0-u;
	
	tmp = (u*momx*momx + v*momy*momy)/2.0 - momx*momy*2.0*vxy/z;
 	if(tmp>0.0) _param.miss = sqrt(tmp);
 	else { _param.miss=0; } // don't want small miss to be invalid here
    }
    else {
	_param.miss = _param.distance;
    }

    //=========================================================================
    // Calculate azwidth

    azwidth2 = (momx2 + momy2 - z);
    if (azwidth2>0.0) _param.azwidth = sqrt(azwidth2);
    else {_param.azwidth=0.0; _param.invalid++;}

    //=========================================================================
    // Calculate psi angle

    tanpsi_numerator = (d+z)*momy + 2.0*vxy*momx;
    tanpsi_denominator = (2*vxy*momy) - (d-z)*momx;

    if (tanpsi_numerator > 1.0) 
	tanpsi_numerator = 1.0;
    else if (tanpsi_numerator<-1.0) 
	tanpsi_numerator = -1.0;
    if (tanpsi_denominator > 1.0) 
	tanpsi_denominator = 1.0;
    else if (tanpsi_denominator<-1.0) 
	tanpsi_denominator = -1.0;
    
    if (fabs(tanpsi_denominator) > EPSILON ) 
	_param.psi = atan2( tanpsi_numerator, tanpsi_denominator );
    else {
	_param.psi = M_PI_2;
	_param.invalid++;
    }

    //    _param.psi = gsl_sf_angle_restrict_pos(_param.psi);

    //=========================================================================
    // Calculate alpha angle

    sinalpha = _param.miss/_param.distance;
    if (-1.0 <= sinalpha && sinalpha <= 1.0)
	_param.alpha = asin( sinalpha );
    else {
	_param.alpha = (sinalpha >=1)? LARGE : -LARGE;
	_param.invalid++;
    }

    //=========================================================================
    // Calculate phi angle (centroid polar angle)

    _param.phi = atan2( _param.centroid.y, _param.centroid.x );
    _param.phi = gsl_sf_angle_restrict_pos( _param.phi );

    //=========================================================================
    // Obtain the max adc values... 

    for (i=0; i<3; i++) {
	_param.max[i] = 0.0;
	_param.index_of_max[i] = 0;
    }
    
    for (iter=cleanpixels.begin(); iter != cleanpixels.end(); iter++) {

	if ( image[*iter] > _param.max[2] ){
	    if ( image[*iter] > _param.max[1] ){
		if ( image[*iter] > _param.max[0] ){
		    _param.max[2] = _param.max[1];
		    _param.max[1] = _param.max[0];
		    _param.max[0] = image[*iter];
		    _param.index_of_max[2] = _param.index_of_max[1];
		    _param.index_of_max[1] = _param.index_of_max[0];
		    _param.index_of_max[0] = *iter;
		}
		else {
		    _param.max[2] = _param.max[1];
		    _param.max[1] = image[*iter];
		    _param.index_of_max[2] = _param.index_of_max[1];
		    _param.index_of_max[1] = *iter;
		}
	    }
	    else {
		    _param.max[2] = image[*iter];
		    _param.index_of_max[2] = *iter;
	    }
	}
    }    

    //=========================================================================
    // Obtain frac[3] values
    
    _param.frac[0] = _param.max[0]/totalsignal;
    _param.frac[1] = _param.frac[0] + _param.max[1]/totalsignal;
    _param.frac[2] = _param.frac[1] + _param.max[2]/totalsignal;

    //=========================================================================
    // 2-D analysis

    _param.asymmetry = 0;

    if ( _2d_is_enabled ) {

	// ---------------------------------------------------------
	// Calculate 2 possible Points of Origin

	cospsi = cos( -_param.psi );
	sinpsi = sin( -_param.psi );

	// First calculate the displacement along the major-axis of the
	// point of origin (don't know whether it is plus or minus)


	if (_param.length > EPSILON) {
	    disp = fabs( _elongation - _elongation *
			_param.width/_param.length );
	}
	else {
	    disp = -LARGE;
	}
	
	// the points of origin in this frame are given by 
	// x' = x +/= disp, y' = y  

	pa.x = -disp;
	pb.x = disp;
	pa.y = 0.0;
	pb.y = 0.0;
						     
	// now put into correct orientation by rotating by psi
	// and then translating to the centroid 

	if (fabs(_param.psi) > EPSILON ) {
	    _param.point_of_origin_a.x = pa.x * cospsi + pa.y*sinpsi;
	    _param.point_of_origin_a.y = -pa.x * sinpsi + pa.y*cospsi;
	    _param.point_of_origin_b.x = pb.x * cospsi + pb.y*sinpsi;
	    _param.point_of_origin_b.y = -pb.x * sinpsi + pb.y*cospsi;
	}	
	else {
	    _param.point_of_origin_a.x = pa.x;
	    _param.point_of_origin_a.y = pa.y;
	    _param.point_of_origin_b.x = pb.x;
	    _param.point_of_origin_b.y = pb.y;
	}

	_param.point_of_origin_a.x += _param.centroid.x;
	_param.point_of_origin_a.y += _param.centroid.y;
	_param.point_of_origin_b.x += _param.centroid.x;
	_param.point_of_origin_b.y += _param.centroid.y;

	// ---------------------------------------------------------
	// Now, we need to know which of the two points is the correct
	// one, so calculate the Asymmetry about one of the points.
	// Put the favored point in "a" and the reject in "b"

	// For the asymmetry calculation, we need to calculate vx3,
	// vx2y, vxy2, vy3 centered around one of the possible points
	// of origin (not about the camera center). That means
	// recalculating all the moments about the new origin.  I
	// don't use the precalculated arrays (_x, _y, etc.) since
	// it's pointless to recalculate ALL the pixels, just the ones
	// in the image. The precalculated arrays are still useful for
	// the next event, so they shouldn't be modified.
	
	double x,y;
	int sign;
	momx=0; momy=0; momx2=0; momy2=0; momxy=0;
	momx3=0; momy3=0; momx2y=0; momxy2=0;

	for (iter=cleanpixels.begin(); iter != cleanpixels.end(); iter++) {
	    x = _x[*iter]-_param.point_of_origin_a.x;
	    y = _y[*iter]-_param.point_of_origin_a.y;
	    momx    += image[*iter] * x;
	    momy    += image[*iter] * y;
	    momx2   += image[*iter] * x*x;
	    momy2   += image[*iter] * y*y;
	    momxy   += image[*iter] * x*y;
	    momx3   += image[*iter] * x*x*x;
	    momy3   += image[*iter] * y*y*y;
	    momx2y  += image[*iter] * x*x*y;
	    momxy2  += image[*iter] * x*y*y;
	}

	momx   /= totalsignal;	
	momy   /= totalsignal;
	momx2  /= totalsignal;
	momy2  /= totalsignal;
	momxy  /= totalsignal;
	momx3  /= totalsignal;
	momy3  /= totalsignal;
	momx2y /= totalsignal;
	momxy2 /= totalsignal;

	vx3  = (momx3  - 3*momx2*momx + 2*momx*momx*momx);
	vy3  = (momy3  - 3*momy2*momy + 2*momy*momy*momy);
	vx2y = (momx2y - 2*momxy*momx + 2*momx*momx*momy - momx2*momy);
	vxy2 = (momxy2 - 2*momxy*momy + 2*momx*momy*momy - momx*momy2);
	
	// Recalculate psi angle from the new point:
	// Really, the angle should be the same as for the origin, but 
	// the sign may be incorrect (since the quandrants are different)

	tanpsi_numerator = (d+z)*momy + 2.0*vxy*momx;
	tanpsi_denominator = (2*vxy*momy) - (d-z)*momx;
	apsi = atan2( tanpsi_numerator, tanpsi_denominator );
	cospsi = cos(apsi);
	sinpsi = sin(apsi);
	
	// Now calculate the asymmetry

	tmp = ( vx3*pow(cospsi,3) + 3*vx2y*cospsi*cospsi*sinpsi 
		+ 3*vxy2*cospsi*sinpsi*sinpsi + vy3*pow(sinpsi,3) );

	// NOTE: unlike quicklook, assymmetry is an absolute value for
	// the correct point of origin (which is always "a")
	// A negative sign just means a swap occured

	tmp = tmp/pow(_param.length,3);
	sign = GSL_SIGN(tmp);
	_param.asymmetry = sign*pow(fabs(tmp), 0.3333333);

	// Swap points of origin if asymmetry is negative, so "a"
	// always contains the correct point of origin

	if (sign == -1 ) {
	    tempcoord = _param.point_of_origin_a;
	    _param.point_of_origin_a = _param.point_of_origin_b;
	    _param.point_of_origin_b = tempcoord;
	}

    }

}


void
clearParam(HillasParameterization &p ) {

    p.event_number =0;
    p.gpstime = 0;
    p.zenith=0;
    p.pixels_in_picture = 0;
    p.centroid.x = p.centroid.y = 0.0;
    p.point_of_origin_a.x = p.point_of_origin_a.y = 0.0;
    p.point_of_origin_b.x = p.point_of_origin_b.y = 0.0;
    p.length = p.width = 0.0;
    p.length_over_size = 0.0;
    p.miss = p.distance = p.azwidth = 0.0;
    p.size = 0.0;
    p.alpha=p.psi=p.phi=0.0;
    p.max[0]=p.max[1]=p.max[2]=0.0;
    p.index_of_max[0]= p.index_of_max[1]= p.index_of_max[2]=0;
    p.frac[0]=p.frac[1]=p.frac[2]=0.0;
    p.asymmetry=0.0;
    p.energy_estimate = 0.0;

}

/**
 * Inserter for outputting an HillasParameterization in text format
 */
ostream& operator<<(ostream &stream,const HillasParameterization &p){

    stream << p.sim 
	   << setw(8) << p.event_number << " " 
	   << setw(8) << (int)p.invalid << " "
	   << setprecision(14) << p.gpstime << " " 
	   << setw(8) << p.centroid.x << " " 
	   << setw(8) << p.centroid.y << " " 
	   << setw(8) << p.point_of_origin_a.x << " " 
	   << setw(8) << p.point_of_origin_a.y << " " 
	   << setw(8) << p.point_of_origin_b.x << " " 
	   << setw(8) << p.point_of_origin_b.y << " " 
	   << setw(8) << p.length << " " 
	   << setw(8) << p.width << " " 
	   << setw(8) << p.miss << " " 
	   << setw(8) << p.distance << " " 
	   << setw(8) << p.azwidth << " " 
	   << setw(8) << p.alpha << " " 
	   << setw(8) << p.psi << " "
	   << setw(8) << p.phi << " "
	   << setw(8) << p.size << " ";

    for (int i=0; i<3; i++) 
	stream  << setw(8) << p.max[i] << " " ;
    
    for (int i=0; i<3; i++) 
	stream << setw(3) <<p.index_of_max[i] << " " ;
    
    for (int i=0; i<3; i++) 
	stream << setw(8) <<p.frac[i] << " " ;
    
    stream  << setw(8) << p.asymmetry<< " " ;
    stream  << setw(8) << p.length_over_size << " ";
    stream  << setw(8) << p.zenith<< " ";
    stream  << setw(8) << p.energy_estimate<< " ";
    stream  << setw(5) << (int)p.telescope_id;
    
    return stream;

}


/**
 * Extractor for HillasParameterizations:
 */
// istream& operator>>(istream &stream, HillasParameterization &p){
    
//     stream  >> p.event_number  
// // 	    >> p.on 
// 	    >> p.invalid 
// 	    >> p.gpstime  
// 	    >> p.centroid.x  
// 	    >> p.centroid.y  
// 	    >> p.point_of_origin_a.x  
// 	    >> p.point_of_origin_a.y  
// 	    >> p.point_of_origin_b.x  
// 	    >> p.point_of_origin_b.y  
// 	    >> p.length  
// 	    >> p.width  
// 	    >> p.miss  
// 	    >> p.distance  
// 	    >> p.azwidth  
// 	    >> p.alpha  
// 	    >> p.psi 
// 	    >> p.phi 
// 	    >> p.size ;
    
//     for (int i=0; i<3; i++) 
// 	stream >> p.max[i]  ;
    
//     for (int i=0; i<3; i++) 
// 	stream >> p.index_of_max[i]  ;
    
//     for (int i=0; i<3; i++) 
// 	stream >> p.frac[i]  ;
    
//     stream >> p.asymmetry ;
//     stream >> p.length_over_size;
//     stream >> p.zenith;

//     return stream;
// }
