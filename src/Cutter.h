// Cutter.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef CUTTER_H
#define CUTTER_H

#include <string>
#include <vector>
#include <list>
#include "Histogram.h"
#include "DataRecords.h"
#include "DataWriter.h"
#include "StarCatalog.h"
#include "ImageAnalyzer.h"

struct RunInfo;
class Image2D;
class PlotMaker;

struct RunStatistics {
    double significance;
    double excess;
};


/**
 * Cuts parameterized data, calculates statistics, and generates final
 * histograms.
 *
 * \todo: this is a badly designed class.  Cutter should be a small
 * object which applies a specific set of cuts (subclasses for other
 * cutting methods). All the analysis should be done elsewhere.
 *
 * 
 * \todo: For multi-telescope analysis, make a telescope_id argument
 * or field and only process events from the specified
 * telescope. 
 */
class Cutter {

 public:

    Cutter();
    ~Cutter();
    
    void enableOutput( bool val=1){_output_is_enabled=val;}
    void enableScaledParameterOutput( bool val=1 ){_output_scaled_parameters=val;}
    void enableFastImageProcessing(bool val=1){
	_fast_image_processing_is_enabled= val;
    }
    void enable2D( bool val=1 ){_is_2d_enabled=val;}
    void enablePointingCheck( bool val=true ){_pointing_check_is_enabled=val;}
    void enableDisplay( bool val=true ){_display_is_enabled=val;}
    void process( RunInfo &ri );
    
    RunStatistics getTrackingStatistics( CutRecord &run, double ratio );
    RunStatistics getPairStatistics( CutRecord &onrun, CutRecord &offrun );
    RunStatistics getTotalPairStatistics();

    void outputStatistics();
    double energyEstimator( double size, double dist );

    CutRecord getTotal( std::vector<CutRecord> &);

    double maxLikelihoodSignif( double n_on, double n_off, double alpha );

    CutRecord cut( RunInfo &, const std::string &, char );
    void generateImage(CutRecord &on ,CutRecord &off, std::string dir,
		       std::string title, bool finalimage=false);

    void radiallyBinPoints(std::vector<Coordinate_t>& poolist,Image2D& destimage, 
			   double radius); 

    void setOutputDir( std::string dirname ){ _outdir = dirname+"/"; }
    void printCutRecordFields( std::ostream &stream );
    void checkPointing();

    void checkDiagnostics();

    double getElongationFactor(double zenith);

    void setTelescopeID( int num ){_telescope_id = num;}


    void clear() {

	_onruns.clear();
	_offruns.clear();
        _trackruns.clear();
	_corr_centroid.clear();
	_sourcenames.clear();
	_onstars.clear();
	_offstars.clear();

    }

    enum types {ON,OFF,TRACK};


 private:
   
    void outputStars();
    void writeStarList( std::list<Star> &starlist, std::string filename );

    Histogram *_alpha_on;
    Histogram *_alpha_off;
    Histogram *_energy_on;
    Histogram *_energy_track;
    Histogram *_energy_off;
    Histogram *_size_on;
    Histogram *_size_off;

    std::vector<CutRecord> _onruns;
    std::vector<CutRecord> _offruns;
    std::vector<CutRecord> _trackruns;
    std::vector< std::vector<Coordinate_t> > _corr_centroid;

    PlotMaker *_xplotmaker;

    std::list<std::string> _sourcenames; 
    std::vector<double> _on_ra;
    std::vector<double> _on_dec;
    std::vector<double> _off_ra;
    std::vector<double> _off_dec;

    std::string _outdir;

    StarCatalog _starcatalog;
    std::list<Star> _onstars;
    std::list<Star> _offstars;
    
    int _image_xdim;
    int _image_ydim;
    double _degperpix;
    int _numruns_on;
    int _numruns_off;
    int _numruns_trk;
    double _track_ratio;
    double _ra,_dec;
    bool _is_2d_enabled;
    bool _pointing_check_is_enabled;
    bool _display_is_enabled;
    bool _output_scaled_parameters;
    bool _fast_image_processing_is_enabled;
    std::string hline;
    int _centroid_range;

    bool _output_is_enabled;
    ParamDataWriter *_writer_on;
    ParamDataWriter *_writer_off;
    
    double _radial_smoothing_factor;


    int _telescope_id;


};

std::ostream& operator<<( std::ostream &stream, CutRecord &c );

#endif
