//
// Camera.cpp
// Karl Kosack <kosack@hbar.wustl.edu>
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <plot.h>
#include <string>

#include "Exceptions.h"
#include "ImageAnalyzer.h"
#include "MuonImageAnalyzer.h"
#include "Camera.h"
#include "PedestalFinder.h"
#include "Log.h"
#include "Config.h"

using namespace std;

/**
 * Create a Camera object.
 * @param filename filename of the camera definition file.
 */ 
Camera::Camera(char *filename) 
{
    
    init(filename);
    
}


Camera::Camera( int npixels, int utdate ) {


    // VERITAS CAMERA:
    if ( 499 <= npixels && npixels <= 501 ) {
	
	init( "veritas_499.camera" );
	_camera_id = ID_V500;
	_lastpix = 499;
	_firstpix = 0;
	_pe2dc = 1.0;  // This need to be updated

	cout << "DEBUG: Camera: this is a VERITAS camera"<<endl;
	
    }
    
    // Whipple cameras:
    else if ( 450 < npixels && npixels < 498 ) {
	init( "490.camera" );
	_camera_id = ID_490;
	_lastpix = 379;  // turn off outer tubes!


	/**
	 *  PFR: change year switch to year+month switch.  telescope seasons
	 *  are fiscal- they start in june of one year, and go until june of the 
	 *  next!  This should only use one pe2dc ratio, but previously it used
	 *  2.
	 */


	int year_month = (int)(utdate/100);
	
	if(year_month >= 109 && year_month <=206){
	  _pe2dc = 0.4467;
	}
	else if(year_month >= 209 && year_month <=306){
	    //_pe2dc = 0.5095*0.6251;
	    _pe2dc = 0.5095;
	}
	else if(year_month >= 309 && year_month <=406){
	  _pe2dc = 0.5352;
	}
	else if(year_month >= 409 && year_month <=506){
	  _pe2dc = 0.4776;
	}
	else{
	  _pe2dc = 0.4489;
	}

	/*
	int year = (int)(utdate/10000);

	switch (year) {
	case 01:
	    _pe2dc = 0.4467;
	    break;
	case 02:
	    _pe2dc = 0.5095;
	    break;
	case 03:
	    _pe2dc = 0.5352;
	    break;
	case 04:
	    _pe2dc = 0.4776;
	    break;
	default:
	    _pe2dc = 0.4489;
	    break;
	}
	*/

    }
    else if ( 300 < npixels && npixels < 340 ) {
	init( "331.camera" );
	_camera_id = ID_331;
	_pe2dc = 1.0;
	_lastpix = 330; 
    }
    else if ( 140 < npixels && npixels < 160 ) {
	init("151.camera");
	_camera_id = ID_151;
	_pe2dc = 1.0;
	_lastpix = 150; 
    }
    else if ( 100 < npixels && npixels < 139 ) {
	init("109.camera");
	_camera_id = ID_109;
 	_pe2dc = 1.0;
// 	_pe2dc = 0.75; // KPK: hand optimized using max2 distribution
	_lastpix = 108; 
    }
    else {
	_camera_id = ID_UNKNOWN;
	_pe2dc = 1.0;
	cout << "DEBUG: Camera: npixels = "<<npixels<<endl;
	throw CriticalAnalysisException("No camera definition exists with specified number of pixels!");
    }

}


/**
 * Instantiate a generic camera, not loaded from a file.  
 */
Camera::Camera( int npixels, Array_t &xcoords, Array_t &ycoords ) {
    _camera_id = ID_UNKNOWN;
    _npix = npixels;
    _firstpix = 0;
    _lastpix = npixels-1;
    _x.resize(npixels);
    _y.resize(npixels);
    _rad.resize(npixels);
    if (xcoords.size() < npixels || ycoords.size()<npixels) {
	throw CriticalAnalysisException("Camera: npixels larger than array size");
    }

    _maxx=-10000.0;
    _minx=10000.0;
    _maxy=-100000.0;
    _miny=10000.0;

    double drad = sqrt(pow(xcoords[2]-xcoords[1],2)
		       + pow(ycoords[2]-ycoords[1],2))/2.0;

    for (int i=0; i<npixels; i++) {
	_x[i] = xcoords[i];
	_y[i] = ycoords[i];
	_rad[i] = drad;
	if (_x[i]>_maxx) _maxx= _x[i];
	if (_y[i]>_maxy) _maxy= _y[i];
	if (_x[i]<_minx) _minx= _x[i];
	if (_y[i]<_miny) _miny= _y[i];
    }

    cout << "DEBUG: rad.size = "<<_rad.size()<<", npix="<<npixels<<endl;
    cout << "DEBUG: camera extent=("<<_minx<<","<<_miny
	 <<")-("<<_maxx<<","<<_maxy<<")"<<endl;

    _sigpix =  sqrt(pow(_x[1]-_x[0],2) + pow(_y[1]-_y[0],2))/8.0;
}


void
Camera::init( char *filename ) {

    _camera_id = ID_UNKNOWN;
    _npix = 0;
    _firstpix = 0;
    _sigpix = 0;
    _pe2dc = 0;
    _maxx = 0;
    _maxy = 0; 
    _minx = 0; 
    _miny = 0;
	
    double dia2;
	
    loadCamera(filename);
    dia2 = pow(_x[1]-_x[0],2) + pow(_y[1]-_y[0],2);
    _sigpix = sqrt(dia2/8.0);
    _lastpix = _npix;


} 


Camera::~Camera() {
   

}


/**
 * Load in the positions and sizes of each camera pixel.  This data is
 * stored in a text file in the proper format and is loaded into
 * position arrays _x[] and _y[]. 
 * @param filename name of camera datafile to load 
  */
void
Camera::loadCamera(char *filename) {

    char id[256];
    string header;
    float  version;
    int i;
    double xx,yy,rr;
    int nneigh,neighid;

    string sfile(filename);
    string sdir;

    sdir = getSupportDir();

    string sfilename(sdir+sfile);
    ifstream camerafile( sfilename.c_str(), ios::in );


    // Check if opened correctly...
    if (camerafile.fail()) {
	cerr << "** Couldn't open camera datafile '" <<filename<<"'"<< endl;
	cerr << "** Set the environment variable WUPARAMDIR to the directory"
	     << endl
	     << "** where it is found (the default is /usr/share/wuparam/)"
	     <<endl;
	exit(1);
    }

    // Read header and version number and see if they're OK
    
    camerafile >> header >> version;

    if (header != "WHIPPLE-CAMERA-DEFINITION") {
	cout << "The camera file you specified, '"<<filename<<"' "
	     << "doesn't appear to be in the correct format!"
	     << endl;
	throw CriticalAnalysisException("Bad Camera Definition File Format!");
    }

    if (version < 1.0) {
	cout << "The version of the camera file '"<<filename<<"' "
	     << "is " << version << ".  Make sure it is in the right format!"
	     << endl;
	throw CriticalAnalysisException( "Wrong Camera Definition version" );
  
    }
    

    // Read number of data points (should be 1st line)
    // and allocate memory:
    
    camerafile >> _npix;

    // resize arrays and neighborlist

    _x.resize( _npix );
    _y.resize( _npix );
    _rad.resize( _npix );
    _neighborlist.resize( _npix );

    // Read in the datafiles... 

    rr = 1.0; 
    i = 0;
    while ( !camerafile.eof() && i<_npix ){
	
	camerafile >> id;
	
	
	if (id[0] == 'R') {
	    // set radius
	    camerafile >> rr;
	}
	else {

	    if (i<_npix) {

		camerafile >> xx >> yy >> nneigh;

		// store the values
		_x[i] = xx;
		_y[i] = yy;
		_rad[i] = rr;

		// add the indices to the neighborlist

		for (int j=0; j<nneigh; j++) {
		    camerafile >> neighid;
		    _neighborlist[i].push_back( neighid );
		}

		// calculate camera viewport size
		if(_maxx < xx) _maxx = xx;
		if(_maxy < yy) _maxy = yy;
		if(_minx > xx) _minx = xx;
		if(_miny > yy) _miny = yy;
    
		i++;

	    }
	    else {
		cout << "** Camera: Pixel # "<< i<< " out of range: " 
		     << " ID = '" << id << "'   (SKIPPED)" << endl;
		continue;
	    }
	}
    }

    _minx = _minx - 2*rr;
    _miny = _miny - 2*rr; 
    _maxx = _maxx + 2*rr;
    _maxy = _maxy + 2*rr; 

    camerafile.close();
    
    cout << "Camera: loaded configuration from '"<<filename<<"'"<<endl;

}


/**
 * Returns a pointer to an STL list<int> of the indices of the neighbors
 * of the specified tube.
 */
list<int>* 
Camera::getNeighborListOfPixel(int n) {

    return &(_neighborlist[n]);

}


/**
 * Returns the estimated dispersion due to the point spread function
 * of the camera at the given zenith angle.
 */
double
Camera::getSigmaPSF( double zen ) {

    //return ((((0.12 - 0.18)/(TWOPI/4.0 - 0.401425728))
    //       *(zen - 0.401425728) + 0.18)/2.354);

    // KPK: changed 0.12 -> 0.1
//     return ((((0.13 - 0.13)/(M_PI_2 - 0.401425728))
//               *(zen - 0.401425728) + 0.13)/2.354);

    return 0.13/2.354;


}



/**
 *  translate the camera coordinates by the specified amount.
 *  This is used to fix offsets in the camera due to telescope
 *  flexing at low elevation, for example. 
 */
void
Camera::
shiftCameraCoordinates( double dx, double dy ) {
    
    _x += dx;    // note: this is a valarray operation
    _y += dy;

    cout << "DEBUG: Camera: camera coordinates shifted by "
	 <<dx<<","<<dy<<endl;

}


/**
 * Loads an array definition from the specified file and creates a
 * telescope array.
 *
 * The array datafile is in the following format:
 *
 * \code
 * WUPARAM-ARRAY-DEFINITION <version number>
 * <x location> <y location> <z location> <dish-size> <camera filename> 
 * ...
 * \endcode
 *
 * for example, 
 * \code
 * WUPARAM-ARRAY-DEFINITION 1.0
 * 0.0   0.0  0.0  10.0 490.camera
 * 80.0  0.0  0.0  12.0 veritas-499.camera 
 * \endcode
 */
TelescopeArray::TelescopeArray(string filename, vector<int> nadc, int utdate) {
    
    string fname = getSupportDir()+filename+".array";
    ifstream infile( fname.c_str() );
    if (infile.fail() ) 
	if (infile.fail()) {
	    cerr << "** Couldn't open array datafile '" <<filename<<"'"<< endl;
	    cerr << "** Set the environment variable WUPARAMDIR to the "
		 << "directory"	 << endl << "** where it is found"
		 <<endl;
	    throw AnalysisException("Couldn't open array definition: "
				    +filename);
    }
    
    string header;
    double version;

    infile >> header >> version;

    // could check if version is wrong here - should do that if you
    // add any colums later which break the format.

    if (header != "WUPARAM-ARRAY-DEFINITION") {
	cerr << "** The array definition '"<<filename<<"' appears to be"<<endl
	     << "** in the wrong format. Please check it"<<endl;
	throw AnalysisException("Bad format for: "+filename);
    }

    TelInfo tel;
    vector<string> tokens;
    string line, camfile;
    int i,linenum=1;

    while (infile &&  !infile.eof() ){
	
	getline( infile, line, '\n' );

	// remove comments
	i=line.find("#");
	line = line.substr(0,i);

	// tokenize
	
	tokenize( line, tokens, " \t," );
	
	// check 

	if (tokens.size() == 0) continue;
	if (tokens.size() != 5) {
	    cerr <<"Line "<<linenum<<": expected 5 columns, got "
		 <<tokens.size()<<endl;
	    throw AnalysisException("Syntax error in "+filename );
	}

	// read data:
	
	tel.x = atof(tokens[0].c_str());
	tel.y = atof(tokens[1].c_str());
	tel.z = atof(tokens[2].c_str());
	tel.r = atof(tokens[3].c_str())/2.0;
	camfile = tokens[4];
	
	int k = _telescope.size();

	if (camfile == "AUTO") {
	    if (nadc.size() < k+1) {
		Logger *logger = Logger::instance();
		logger->printf("Data reports %d telescopes, while array"
			       " specifies %d. Using default 499 pix camera"
			       " for telescope %d", 
			       nadc.size(), _telescope.size(), k);
		tel.cam = new Camera( 499, utdate );
	    }
	    else {
		try{
		    tel.cam = new Camera( nadc.at(k), utdate );
		}
		catch (out_of_range &e) {
		    cout << "DEBUG: nadc wrong dimension: " << e.what() << endl;
		    cout << "DEBUG: k="<<k<<", nadc.size()="<<nadc.size()<<endl;
		    throw e;
		}

	    }
	}
	else {
	    tel.cam = new Camera( (char*) camfile.c_str() );
	}

	// add the camera to the array:
	
	_telescope.push_back( tel );
	
	linenum++;

    }


    // compute array bounding box:

    _minx =10000;
    _miny =10000;
    _maxx =-10000;
    _maxy =-10000;

    double x0,y0,x1,y1;

    for (int i=0; i<getNumTelescopes(); i++) {

	cout << "TEL "<<i<<": "
	     <<getLocationX(i) << ","
	     <<getLocationY(i) <<endl;

	x0 = _telescope[i].x - _telescope[i].r;
	x1 = _telescope[i].x + _telescope[i].r;
	y0 = _telescope[i].y - _telescope[i].r;
	y1 = _telescope[i].y + _telescope[i].r;

	_minx = fmin( _minx, x0 );
	_miny = fmin( _miny, y0 );
	_maxx = fmax( _maxx, x1 );
	_maxy = fmax( _maxy, y1 );
    }

    // make it square!
//     _minx = fmin( _minx, _miny );
//     _miny = fmin( _minx, _miny );
//     _maxx = fmax( _maxx, _maxy );
//     _maxy = fmax( _maxx, _maxy );

    cout << "TelescopeArray: loaded "<<getNumTelescopes()<<" telescopes"<<endl;

}


TelescopeArray::~TelescopeArray() {

    // clean up all the cameras!
    for (int i=0; i<_telescope.size(); i++) {
	if (_telescope[i].cam) delete _telescope[i].cam;
    }

}
