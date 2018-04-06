
// ImageAnalyzer.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef IMAGEANALYZER_H
#define IMAGEANALYZER_H

#include <vector>
#include "Camera.h"
#include "Types.h"
#include "DataRecords.h"

const double LARGE   =  1e20;
const double EPSILON =  1e-50;

/**
 * Information common to all parameterizations
 */
struct ImageParameterization {

    short telescope_id;         //!< telescope id number
    double osctime;             //!< oscillator time
    double gpstime;             //!< calibrated gps time
    double livetime;            //!< livetime of run
    short  invalid;             //!< was this event properly parameterized?
    int    event_number;        //!< Event number for debugging
    
    SimShowerRecord sim;        //!< simulation info

};

/**
 * Parameterization of an image, produced by the ImageAnalyzer object.
 * I kept this a pure struct, not a class so that the gsl N-tuple
 * routines can process it.
 */
struct HillasParameterization : public ImageParameterization {

    Coordinate_t centroid;    	//!< Center position of shower 
    Coordinate_t point_of_origin_a;	//!< First point of origin (x,y)
    Coordinate_t point_of_origin_b;	//!< second point of origin (x,y)
    double length;		//!< Length of ellipse
    double width;		//!< width of ellipse
    double size;                //!< total signal
    double miss;		//!< miss parameter
    double distance; 		//!< distance to centroid
    double azwidth;             //!< used by certain analysis techniques
    double alpha;		//!< Alpha angle (-Pi/2..Pi/2)
    double length_over_size;    //!< length over size
    double psi;			//!< Angle betw major axis and x-axis (-PI..PI)
    double phi;                 //!< Angle that centroid makes with x-axis
    double max[3];		//!< The three largest adc values
    int    index_of_max[3]; 	//!< The tube indices of the max[3] values
    double frac[3]; 		//!< fraction of maximum digital counts
    int    pixels_in_picture;   //!< number of tubes in the picture
    double asymmetry; 		//!< measure of shower skew
    int    on;                  //!< is this from an on run? (added by henric)
    double zenith;              //!< zenith angle of event
    double energy_estimate;     //!< estimated energy - computed at time of cut

    friend std::ostream &operator<<( std::ostream &stream,
				     const HillasParameterization &p);
/*     friend std::istream &operator>>( std::istream &stream, */
/* 				     HillasParameterization &p); */
    

};


void clearParam( HillasParameterization &p );


/**
 * Base class for ImageAnalyzers - components that look at an
 * individual image and derive a set of parameters.
 */
class ImageAnalyzer {

 public:

    ImageAnalyzer( Camera *cam );
    ~ImageAnalyzer() {;}

    int    getNumPixels() { return _nadc; }
    virtual void parameterize( const Array_t &image,
			       std::vector<int> &cleanpixels)=0;

    
    
protected:

    Array_t _x, _y;             // the camera coorinates with the origin moved
    Array_t _x2, _y2, _xy;	// precalculated values
    
private:

    void init( Array_t &xcoords, Array_t &ycoords, int n );
    double _xorigin;		// x origin of coord system
    double _yorigin;		// y origin of coord system 
    int _nadc;


};


/**
 * Calculates the Hillas parameters and 2D point of origin of a
 * cleaned image. 
 */
class HillasImageAnalyzer :public ImageAnalyzer {
    
 public:
    
    HillasImageAnalyzer( Camera *cam, double elongation=1.68 );
    void parameterize( const Array_t &image, std::vector<int> &cleanpixels );
    void setElongationFactor( double e ){ _elongation = e; }

    HillasParameterization getHillasParameters() {
	return _param;
    }
    void enable2D( bool val=true )      { _2d_is_enabled = val;}
    bool is2DEnabled()                  { return _2d_is_enabled; }


 private:

    double _elongation;         // elongation factor
    int _nadc;			// number of adc's (from Camera)
    bool _2d_is_enabled;	// enable the 2d analysis
    HillasParameterization _param;// The current parameterization


};



#endif
