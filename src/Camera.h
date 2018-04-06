//
// camera.h
//
// Karl Kosack <kosack@hbar.wustl.edu>
//

#ifndef _CAMERA_H
#define _CAMERA_H

#include "Types.h"
#include "PedestalFinder.h"
#include <cstdio>
#include <vector>
#include <list>
#include <string>

class HillasParameterization;
class MuonParameterization;

/**
 * Contains all operations pertaining to camera data. This includes
 * loading the camera pixel positions from a file and building the
 * pixel neighborlist.
 *
 *
 *  \todo: make separate classes for each camera 
 */
class Camera  {

 public:

    Camera( int npixels, Array_t &xcoords, Array_t &ycoords );
    Camera(char *);
    Camera( int npixels,int utdate=-1 );
    ~Camera();

    void init( char * );

    void   shiftCameraCoordinates( double dx, double dy );
    double getSigmaPSF(double zen);           //!< Point spread dispersion
    double getSigmaPixel(){ return _sigpix; } //!< Finite pix size dispersion
    double getPEToDC(){ return _pe2dc; }      //!< pe/dc ratio

    int  getNumPixels()   { return _npix; }   //!< Number of pixels in camera
    int  getFirstPixel() { return _firstpix;};//!< first pixel # to use
    int  getLastPixel() { return _lastpix;};  //!< last pixel # to use
    int  getCameraID() { return _camera_id; } //!< Camera ID number

    Array_t& xCoords() { return _x; }         //!< Array of x coordinates
    Array_t& yCoords() { return _y; }         //!< Array of y coordinates
    Array_t& radii()   { return _rad; }       //!< Array of tube radii

    double getMinX(){ return _minx; }
    double getMinY(){ return _miny; }
    double getMaxX(){ return _maxx; }
    double getMaxY(){ return _maxy; }

    std::list<int>* getNeighborListOfPixel(int n);

    static const int AUTOSCALE = -1;
    enum CameraID { ID_490, ID_331, ID_151, ID_109, ID_V500, ID_UNKNOWN };

 private:

    Array_t _x;
    Array_t _y;
    Array_t _rad;

    std::vector< std::list<int> > _neighborlist;
				// vector of tube neighborlists;

    void    loadCamera(char *filename);

    int     _camera_id;      	        // camera id number
    int     _npix;		        // number of camera pixels
    int     _firstpix;                  // first usable pixel
    int     _lastpix;                   // last usable pixel 

    double _sigpix;                     // error due to finite pixel size
    double _pe2dc;                      // p.e. to d.c. ratio
    
    double _maxx, _maxy, _minx, _miny;	// camera display size
    FILE *_fp;

};



/**
 * Stores information for an array of telescopes, including Camera objects
 */
class TelescopeArray {

 public:

    struct TelInfo {
	double x;
	double y;
	double z;
	double r;
	Camera *cam;
    };

    TelescopeArray( std::string filename, std::vector<int> nadc, int utdate );
    ~TelescopeArray();

    int     getNumTelescopes() { return _telescope.size(); }
    Camera* getCamera( int telescope_id ) {
	return _telescope[telescope_id].cam; 
    }
    
    double getLocationX( int id ) { return _telescope[id].x; }
    double getLocationY( int id ) { return _telescope[id].y; }
    double getLocationZ( int id ) { return _telescope[id].z; }

    double getMinX(){ return _minx ;}
    double getMaxX(){ return _maxx ;}
    double getMinY(){ return _miny ;}
    double getMaxY(){ return _maxy ;}

 private:

    std::vector<TelInfo> _telescope;

    double _minx, _miny, _maxx, _maxy;

};


#endif
