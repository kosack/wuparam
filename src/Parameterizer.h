// Parameterizer.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef _PARAMETERIZER_H
#define _PARAMETERIZER_H

#include "Types.h"
#include "Config.h"
#include "ImageAnalyzer.h"
#include "MuonImageAnalyzer.h"
#include "Histogram.h"
#include "ImageCleaner.h"
#include "DataReader.h"

#include <string>
#include <vector>

struct RawEventRecord;
class GaussianDeviate;
class RawHeaderRecord;
class RawDataReader;
class ParamDataWriter;
struct Pedestal;


/**
 * Contains all functions pretaining to the parameterization of raw data.
 */
class Parameterizer {

 public:

    Parameterizer();
    ~Parameterizer();

    void process( RunInfo & );
    void setLatLonInRadians( double lat, double lon);
    void enableDisplay( bool val=true ){ _is_display_enabled=val; }  
    void enableVerbose( bool val=true ){ _is_verbose_enabled=val; }
    void enableInteractive( bool val=true ){ _is_interactive=val; }
    void enableZenithCorrection(bool val=true){_is_zenithcorrect_enabled=val;}
    void enableCleanup( bool val=true ) { _is_cleanup_enabled=val; }
    void enableOverwrite(bool val=true){_is_overwrite_enabled=val;}
    void printSummary( std::ostream &stream, RunInfo &ri, std::string &id, 
		       RawHeaderRecord &header, ParamDataWriter *writer, 
		       RawDataReader *data, TelescopeArray &array);

 private:


    /**
     *  All per-telescope information for the array
     */
    struct TelescopeData {
	HillasImageAnalyzer *analyzer;
	MuonImageAnalyzer *muonanalyzer;
	ImageCleaner  *cleaner;
	std::vector<Pedestal> combinedpeds;
	std::vector<double> sigmaped;
	Array_t gains;
	std::vector<Pedestal> peds;        // Pedestal array 
	std::vector<Pedestal> padpeds;        // Pedestal array 
	
	Histogram *lensize;
	Histogram *pictubes;
	Histogram *tubehits;
	Histogram *phihist;
	Histogram *alphahist;
	Histogram *lengthhist;
	Histogram *sizehist;
	Histogram *disthist;
	Histogram *widthhist;
	Histogram *rawrate;
	Histogram *deltat;
	Histogram *max2hist;
	Image2D  *centroids;
       
        Histogram *musoalhist;
        Histogram *musizehist;
        Histogram *mugainhist;
        Histogram *mulensizehist;	
        Histogram *muonnesshist;	

	

	TelescopeData() {

	    analyzer =NULL;
	    muonanalyzer =NULL;
	    cleaner =NULL;
	    lensize =NULL;
	    pictubes =NULL;
	    phihist =NULL;
	    alphahist =NULL;
	    lengthhist =NULL;
	    sizehist =NULL;
	    disthist =NULL;
	    widthhist =NULL;
	    rawrate =NULL;
	    max2hist =NULL;
	    centroids =NULL;
	    musoalhist =NULL;
	    musizehist =NULL;
	    mugainhist =NULL;
	    mulensizehist =NULL;
	    muonnesshist =NULL;
	}

	~TelescopeData() {
	    if (analyzer) delete analyzer;
	    if (muonanalyzer) delete muonanalyzer;
	    if (cleaner) delete cleaner;
	    if (lensize) delete lensize;
	    if (pictubes) delete pictubes;
	    if (phihist) delete phihist;
	    if (alphahist) delete alphahist;
	    if (lengthhist) delete lengthhist;
	    if (sizehist) delete sizehist;
	    if (disthist) delete disthist;
	    if (widthhist) delete widthhist;
	    if (rawrate) delete rawrate;
	    if (max2hist) delete max2hist;
	    if (centroids) delete centroids;
	    if (musoalhist) delete musoalhist;
	    if (musizehist) delete musizehist;
	    if (mugainhist) delete mugainhist;
	    if (mulensizehist) delete mulensizehist;
	    if (muonnesshist) delete muonnesshist;
	}

    };

    void processOne( RunInfo &, std::string &id );
    void updateAngles( double , const RawHeaderRecord & );
    void derotatePoint(double, Coordinate_t &);
    void generatePlots( std::string &, std::string &, std::string &, int );   
    void mapSkyBrightness( const std::string &id, 
			   const std::vector<Pedestal> &onpeds, 
			   const std::vector<Pedestal> &offpeds, 
			   const Array_t &gain,
			   const RawHeaderRecord &header, 
			   double avg_gpstime, bool derotation, Camera *cam,
			   int telescope_id);

    //    Camera *_cam;
    GaussianDeviate *_gaussdev;

    std::string _wupdir;

    bool _is_verbose_enabled;
    bool _is_display_enabled;
    bool _is_interactive;
    bool _is_zenithcorrect_enabled;
    bool _is_first_run;
    bool _is_overwrite_enabled;
    bool _is_cleanup_enabled;
    
    double _latitude; 	                // latitude of telescope
    double _longitude;			// longitude of telescope

    double _theta;                      // derotation angle (see updateAngles)
    double _az;                         // current azimuth (see updateAngles)
    double _el;                         // currentelevation (see updateAngles)

    int _event_number;	                // global event number



};

double get_time();

#endif
