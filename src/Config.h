// 
// Config.h
//
// Karl Kosack <kosack@hbar.wustl.edu>
//


#ifndef _CONFIG_H
#define _CONFIG_H

#include <fstream>
#include <string>
#include <queue>
#include <vector>
#include <map>
#include "Types.h"

using std::string;
using std::vector;

// Helper functions

void removeWhiteSpace( string &str );
string getSupportDir();
void tokenize(const string& str, vector<string>& tokens,
	      const string& delimiters=" ");


/**
 * Describes range of values of a parameter to retain after cuts.
 */ 
struct Cut {
    double lower;
    double upper;

    Cut():lower(0),upper(1e9){;}


};

/**
 * Describes which cuts should be applied to the data.
 * Fields with upper and lower values are arrays with 2 elements.
 */ 
struct CutInfo {

    double tracking_ratio;
    double elongation;          //!< elongation factor to use for 2D
    Cut alpha;
    Cut distance;
    Cut length;
    Cut width;
    Cut size;
    Cut lensize;
    Cut max1;
    Cut max2;
    Cut max3;
    Cut frac3;
    Cut asymmdist;
    double dlength; // for zcutter
    double dwidth;  // for zcutter

    double smoothing_radius;    //!< radial smoothing radius for 2D
				//!< and offset analysis
    bool  radial_analysis;     //!< offset instead of alpha-cut analysis
    Coordinate_t radial_offset;      //!< offset for radial cut

    bool alignment;             //!< Do 2d alignment when stacking runs
    Coordinate_t align_offset;       //!< alignment offset in x direction

    CutInfo() {
	// set defaults
	elongation = 1.68;
	tracking_ratio = 3.0;
	alpha.lower = 0; alpha.upper = 1e9; 
	distance.lower = 0; distance.upper = 1e9; 
	length.lower = 0; length.upper = 1e9; 
	width.lower = 0; width.upper = 1e9; 
	size.lower = 0; size.upper = 1e9; 
	lensize.lower = 0; lensize.upper = 1e9; 
	max1.lower = 0; max1.upper = 1e9; 
	max2.lower = 0; max2.upper = 1e9; 
	max3.lower = 0; max3.upper = 1e9; 
	frac3.lower = 0; frac3.upper = 1e9; 
	asymmdist.lower = 0; asymmdist.upper = 1e9; 
	dlength=1e9;
	dwidth=1e9;
	smoothing_radius = 0.25;
	radial_analysis = false;
	radial_offset.x = 0.0;
	radial_offset.y = 0.0;
	alignment = false;
	align_offset.x = 0.0;
	align_offset.y = 0.0;
    }

    friend std::ostream &operator<<( std::ostream &stream,const CutInfo &c);

};


/**
 * Basic structure returned by a Config object.  Describes all
 * information about the run that is collected from the configuration
 * file - such as override values and data locations 
 *
 * \todo: Make RunInfo serializable so it can be sent as an MPI message?
 *
 * \todo: Make RunInfo saveable (i.e. write out a config file
 * containing all the settings.  This can be written to each run
 * output dir so you can check what settings were used at the last
 * parameterization.  You can then use a Config object to parse the
 * saved settings and check if the RunInfo has changed as a robust
 * check for reparameterizing.  Save() load() compare()
 *
 */
class RunInfo {
    
 public:

    char type;			//!< type of run (on/off vs tracking)
    string datadir;		//!< directory for data files
    string cachedir;	        //!< directory to write ped/gain data
    string sourcename;		//!< name of source
    string onid;		//!< id of on/tracking run
    string offid;		//!< id of off run for padding 
    string n2id;		//!< id of nitrogen gain run
    string utdate;		//!< UT date
    string arrayconfig;         //!< name of array config file

    double ra;			//!< override value source RA
    double dec;			//!< override value source declination
    double zenithoverride;      //!< override value for zenith angle
    int    ntubes;		//!< override for number of tubes
    double picthresh;		//!< picture threshold value
    double bndthresh;		//!< boundary threshold value
    double padlevel;            //!< pedvar to pad to if there is no off run

    double utbase;              //!< Reference time for lightcurve binning

    bool   padding;		//!< flag to enable padding
    bool   pad_track_with_run;  //!< flag to enable real padding for track runs
    bool   derotation;          //!< flag to enable derotation
    bool   muoncalibrate;   	//!< flag to enable muon calibration
    int    outtype;             //!< (see Config::OutType)
    CutInfo cuts;		//!< How this file should be cut
    double muon_threshold;      //!< Muonness threshold for muon selection
    int    cuttype;		//!< type of cutter to use

    string mcdatabase; 		//!< monte-carlo database for Fitter
    string energyestimator;     //!< energy estimator type string
    double emin_tev;            //!< minimum energy for energy spectrum 
    double emax_tev;            //!< maximum energy for energy spectrum 
    int    num_e_bins;          //!< number of energy bins for spectrum
    bool   tubemask[1024];      //!< bitmask of tubes to turn off manually
    
    bool camera_offset_analysis;//!< perform parameterization with offset
    Coordinate_t camera_offset;

    int max_ped_events;         //!< Maxum number of events to use for peds

    bool image_filter;          //!< Smooth images during
				// parameterization instead of using
				// standard cleaner
    double filter_amount;       //!< afiltr
    int filter_subdivide;       //!< msub
    int filter_iter;            //!< iterations for derivative calculation
    double filter_smoothness;	//!< spas 
    double filter_percent;
    double filter_threshold;   //!< Filter picture threshold (in sigma)
    
};


class TrackingRunInfo : public RunInfo {




};

class PairRunInfo : public RunInfo {


};


/**
 * Loads Analysis configuration information from a text file. 
 */
class Config {
    
 public:
    
    Config(char *filename);
    ~Config();

    enum RunType { TRACK,ONOFF };
    enum OutType { NTUPLE, TEXT };
    enum CutType { SUPERCUTTER, ZCUTTER, EZCUTTER,SPECTRALCUTTER };

    void    writeSample(char *file);
    RunInfo getNextRun();
    std::deque<RunInfo>& getRunQueue(){ return _runqueue; }
    bool    isDone() { return _runqueue.empty(); }
    int     getRunCount(){ return (int)_runqueue.size();}


 private:


    void openFailed(char *filename);
    void parse(char *filename, int level=0 ) ;
    void addRun();
    void loadAlignmentMap( std::string file );

    std::deque<RunInfo> _runqueue;
    RunInfo _curinfo;
    std::map<std::string,Coordinate_t> _alignmap;
    
};




#endif
