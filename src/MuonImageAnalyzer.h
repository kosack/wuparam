#ifndef MUONIMAGEANALYZER_H
#define MUONIMAGEANALYZER_H

#include <vector>
#include "Camera.h"
#include "Types.h"
#include "ImageAnalyzer.h"
#include "PlotMaker.h"

static const int MAX_PLOT_CENTERS = 10000;

// Useful constants:

static const short NX = 48;
static const short NY = 48;
static const int NXD2 = NX/2;
static const int NYD2 = NY/2;
static const double RSQMAX = 4.0;
static const double SATURATION_RADIUS = 1.1459; // max cherenkov radius
static const double DRAD_COARSE = 0.24; // anulus radius width
static const double DRAD_FINE = 0.24; // finer anulus radius width
//  static const double DRAD = 0.15; // anulus radius width
static const double SQRT12 = 3.4641;
static const double DX = 0.125; // this should be comparible to the 
static const double DY = 0.125; // pixel size



/**
 * Describes the muon-like characteristics of an image
 */
struct MuonParameterization : public ImageParameterization {

    double  radius ;   		//!< Average radius 
    Coordinate_t ringcenter ;	//!< Average position
    double  arcstrength ;	//!< fract of centers in peak bin
    double  gain ;		//!< Total picture light per all pixels in arc
    double  mugain ;		//!< Light in annulus corrected for spillover
    double  ringfrac;         	//!< fraction of ring falling on camera
    double  soal ;		//!< Signal over arclength
    double  muskew ;		//!< centroid->center length over ring radius
    double  arclen ;		//!< Length of the muon arc
    double  philo;              //!< starting angle of arc
    double  phihi;              //!< ending angle of arc
    double  phimid;             //!< midpoint angle
    Coordinate_t plot[MAX_PLOT_CENTERS]; //!< Array of ring center positions 
    int	nplotpoints ;		//!< number of ring positions
    double muonness; 		//!< measure of muon-like characteristics
    double smoothness;          //!< measure of angular smoothness of ring
    double smoothness_var;      //!< measure of angular smoothness variation of ring
    double xcs;
    double ycs;
    double rspread;

    void clear() {
	
	radius = 0;
	ringcenter.x  = ringcenter.y = 0;
	arcstrength = 0;
	gain = 0;
	mugain = 0;
	ringfrac = 0;
	soal = 0;
	muskew = 0;
	arclen = 0;
	nplotpoints = 0;
	muonness =0;
	smoothness=0;
	smoothness_var=0;
	philo=0;
	phihi=0;
	phimid=0;
	rspread=0;
	xcs=0;
	ycs=0;

    };

};

std::ostream &operator<<( std::ostream &stream,const MuonParameterization &p);


/**
 * Analyzes the muon-like properties of an image. In particular, it
 * detects muon arcs and can calculates factors useful for absolute
 * gain calibration.
 *
 *  The algorithm goes as follows:
 *
 *  - For each unique triplet of PMTs (N*(N-1)*(N-2)/6 in total)
 *  calculate the center and radius of the uniquely defined circle.
 *  
 *  - OR the bits for PMTs which are elements of this triplet into the
 *  appropriate element (determined by the ring center) of a 2-d array
 *  of bitmasks, and increment the corresponding element of the array
 *  ntriplets[x,y]. Also increment the cumulative ring radius and
 *  cumulative x and y-center position arrays rad[x,y] xsum[x,y],
 *  ysum[x,y].  
 * 
 *  - Divide xsum, ysum, radsum by ntriplets in the end.  Then do a
 *  boxcar average over neighboring elements in the array ntriplets,
 *  and x,y,rad arrays.  
 *
 *  - Find the maximum of the boxcar averaged array ntripbox(x,y) and
 *  calculate the values xave,yave,rave from these.  
 * 
 *  - Calculate a bitmask (for which tubes are in the arc) by ORing
 *  together the bitmasks of the peak element and its neighbors.  From
 *  the resulting bitmask calculate the signal sum for these pmts and
 *  divide by the number of pmts.  
 *
 *  This gives a value proportional to the pe/dc ratio for this arc.
 *  A comparison of this value to the value derived for big rings
 *  (from the hadronicity) can be combined with limits on the ring
 *  radius to derive the likelihood that this event is in fact a muon
 *  arc.
 *
 */
class MuonImageAnalyzer : public ImageAnalyzer {
    
 public:
    
    MuonImageAnalyzer( Camera *cam );
    ~MuonImageAnalyzer(){;}
    
    void parameterize( const Array_t &image, std::vector<int> &cleanpixels );
    void enablePlotPoints( bool val=true ) { _plot_points_is_enabled =val;}
    void enableCalibration( bool val=true ) { _plot_points_is_enabled =val;}

    /**
     * Get the parameters calculated in the parameterize() routine.
     * \returns a MuonParameterization containing all the parameters
     */
    MuonParameterization getMuonParameters() { return _mparam; }



 private:
    
    // Some structures for the internal triplet lookup table
    // Since these won't ever be used outside this class, they are private

    struct Triplet {
	float	x0 ;
	float	y0 ;
	short     ix ;
	short	iy ;
	float	r ; 
    };
    
    struct WordBitMask {
	int		word ;
	unsigned int	bit ;
    };

    // useful typedefs for the lookup-tables
    
    typedef std::vector<std::vector< std::vector<double> > > Grid_3d;
    typedef std::vector<std::vector< std::vector<Triplet> > > Grid_3d_trip;
    typedef std::vector<std::vector<std::vector<unsigned int> > > Grid_3d_uint;
    typedef std::vector<std::vector<double> > Grid_2d ;
    typedef std::vector<std::vector<int> > Grid_2d_int ;
    


    // helper functions

    void getCenter( int,int,int,Triplet & t );
    
    // Parameters

    MuonParameterization _mparam; // precalculated values
    //    Grid_3d_trip	_centers ;
    Grid_3d_uint	_pmtmask ;
    Grid_2d		_xsum ;
    Grid_2d		_ysum ;
    Grid_2d		_rsum ;
    Grid_2d_int		_ntrip ;
    vector<WordBitMask>	_single_pmt_mask ;     
    int			_nwords ;
    double		_apix ;		// Approx. area of a pixel
    bool _plot_points_is_enabled;
    bool _calibration_is_enabled; 

    PlotMaker *pm;

};


#endif
