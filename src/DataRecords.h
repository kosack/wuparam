// DataRecords.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef DATARECORDS_H
#define DATARECORDS_H

// HERE ARE ALL THE DATA STRUCTURES NEEDED TO READ/WRITE DATA:

#include "Image2D.h"
#include <vector>

/**
 * Records are readable/writable (i.e. serializable) representations
 * of data. For example, a DataReader should read data from a file and
 * return a Record for each event read; the corresponding DataWriter
 * would write a Record to disk.
 */ 
struct Record {

};


/**
 * Run information (header) from raw datafiles
 *
 * \todo: add fields for picture/boundary threshold so they can be
 * checked by wucut - don't want to accidently use cuts which were
 * designed for one threshold after you parameterize with another!
 */ 
struct HeaderRecord : public Record {

    int num_telescopes;		//!< number of telescopes in the array
    std::vector<int> nadc;	//!< number of pixels in each camera 
    double ra;			//!< ra of run
    double dec;			//!< dec of run
    double starttime;		//!< starting time (MJD) 
    double endtime;		//!< ending time (MJD)
    double average_elevation;	//!< Average elevation (duh!)
    int windowsize;		//!< size of integration window (if FADC data)
    std::string sourcename;	//!< name of the source 
   
};


struct SimShowerRecord  {
    int event_number;
    int primary_type;
    double primary_energy;
    Coordinate_t impact_parameter;
    Coordinate_t direction_cos;
};



/**
 * Event information
 */
struct EventRecord : public Record {


};



/**
 *  One event from a raw data file
 */
struct RawEventRecord : public EventRecord {

    int     type;
    double  gpstime;
    double  osctime;
    double  livetime;
    short   elevation;
    short   azimuth;
    int     telescope_id;
    Array_t adc;   


};


struct RawHeaderRecord : public HeaderRecord {

};

/**
 * Output from Cutter for each run...
 */ 
struct CutRecord : public Record {

    std::string runid;
    int total; 		     	// total number of events
    int valid;			// valid events
    int trigger;		// events passing trigger cut
    int shape;			// events passing shape cut
    int orientation;		// events passing orientation cut
    int trackoff;		// events that are tracking "off" events
    double duration;            // Run duration in minutes
    double average_elevation;   // average zenith angle

    Image2D im2d;               // 2-D image
    Image2D skybright;          // Sky brightness
    Image2D tubeoff;            // Tubeoffness

    std::vector<Coordinate_t> poolist;// point-of-origin list for later binning

    CutRecord() {
	duration=0;
	total=0;
	valid=0;
	trigger=0;
	shape=0;
	orientation=0;
	trackoff =0;
	runid="unknown";
	im2d.setCoordinateBox( -im2d.getXDim()*0.1/2.0, 
			       -im2d.getYDim()*0.1/2.0,
			       im2d.getXDim()*0.1/2.0,
			       im2d.getYDim()*0.1/2.0 );
	skybright.setCoordinateBox( -im2d.getXDim()*0.1/2.0, 
			       -im2d.getYDim()*0.1/2.0,
			       im2d.getXDim()*0.1/2.0,
			       im2d.getYDim()*0.1/2.0 );
	tubeoff.setCoordinateBox( -im2d.getXDim()*0.1/2.0, 
			       -im2d.getYDim()*0.1/2.0,
			       im2d.getXDim()*0.1/2.0,
			       im2d.getYDim()*0.1/2.0 );
	poolist.reserve(2000);
			       
    }

    ~CutRecord() {

    }

    void addValuesFrom( CutRecord &cr ) {
	total += cr.total;
	valid += cr.valid;
	trigger += cr.trigger;
	shape += cr.shape;
	orientation += cr.orientation;
	trackoff += cr.trackoff;
	duration += cr.duration;
	average_elevation = (average_elevation + cr.average_elevation)/2.0;
	im2d.addImage( cr.im2d );
	skybright.addImage( cr.skybright );
	tubeoff.addImage( cr.tubeoff );

	std::vector<Coordinate_t>::iterator coord;
	for(coord=cr.poolist.begin(); coord != cr.poolist.end(); coord++){
	    poolist.push_back( *coord );
	}
    }

};

/**
 * Header Information read from the raw data
 */
class RawRunInfo {

 public:
    std::string sourcename;
    std::string runid;
    int nadc;
    double ra;
    double dec;
    double starttime;

};

/**
 * Adds info calculated in the parameterize routine
 */
class ParamRunInfo : public RawRunInfo {

 public:

    double endtime;
    double average_elevation;
   
};

/**
 * ParamRunInfo + cut information
 */
class CutRunInfo : public ParamRunInfo {
    
 public:
    
    int total; 		     	// total number of events
    int valid;			// valid events
    int trigger;		// events passing trigger cut
    int shape;			// events passing shape cut
    int orientation;		// events passing orientation cut
    int trackoff;		// events that are tracking "off" events
    double duration;            // Run duration in minutes

};


#endif
