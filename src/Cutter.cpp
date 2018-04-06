// Cutter.cpp
//
// Karl Kosack <kosack@hbar.wustl.edu>
// Checked Li and Ma formula; some slight modifications JB 030617

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <list>
#include <cmath>
#include <sys/types.h>         	// For POSIX mkdir
#include <sys/stat.h>		// For POSIX mkdir
#include <unistd.h>
#include <gsl/gsl_math.h>

#include "Exceptions.h"
#include "Config.h"
#include "DataReader.h"
#include "DataWriter.h"
#include "ParamDataReader.h"
#include "Histogram.h"
#include "ProgressBar.h"
#include "Cutter.h"
#include "Image2D.h"
#include "StarCatalog.h"
#include "PlotMaker.h"
#include "Log.h"
#include "Types.h"
#include "EnergySpectrum.h"
#include "SuperCutter.h"
#include "ZCutter.h"
#include "EZCutter.h"
#include "SpectralCutter.h"

using namespace std;
using namespace Cuts;


Cutter::Cutter() 
    : _image_xdim(39), _image_ydim(39), _degperpix(0.1), 
      _is_2d_enabled(false), _pointing_check_is_enabled(true), 
      _display_is_enabled(true),
      _outdir("Totals/"), _centroid_range(4), _output_is_enabled(false), 
      _output_scaled_parameters(false), 
      _fast_image_processing_is_enabled(false),_telescope_id(0)
{

    double minx, miny, maxx,maxy;

    hline = "==============================================================================";

    // JB 030617 Note that the way histogramming currently is implemented,
    // minx is the lower bound on the minimum x-bin, NOT THE CENTER x-value
    // OF THAT BIN.  Similarly maxx is the upper bound on the largest
    // x-bin.
    minx = -_image_xdim*_degperpix/2.0;
    maxx = _image_xdim*_degperpix/2.0;
    miny = -_image_ydim*_degperpix/2.0;
    maxy = _image_ydim*_degperpix/2.0;

    _energy_on = _energy_off = _energy_track = NULL;

    _alpha_on = new Histogram(18,0,90,"alpha-on");
    _alpha_off = new Histogram(18,0,90,"alpha-off");

    _size_on = new Histogram( 15,1,4,"log10 Size-ON");
    _size_off = new Histogram( 15,1,4,"log10 Size-OFF");

    // Set up output directory:

    if (mkdir(_outdir.c_str(), S_IRUSR|S_IWUSR|S_IXUSR|S_IXGRP|S_IRGRP )) {
	if (errno != EEXIST) {
	    perror("mkdir");
	    throw CriticalAnalysisException("Couldn't create "+_outdir);
	}
    }

    // Load star catalog

    _starcatalog.load( getSupportDir()+"YBS.edb" );
    _starcatalog.load( getSupportDir()+"SAO.edb" );
    _onstars.clear();
    _offstars.clear();

    _radial_smoothing_factor = 2.5*_degperpix; // 2.5*_degperpix
    //    _radial_smoothing_factor = 0.325; // 2.5*degperpix


    _xplotmaker = NULL;


    // initialize the writer.  Eventually, there should be an option
    // to use text vs ntuple, but that may not matter (really should
    // implement an ntuple2text converter). 

    _writer_on = new ParamDataWriterNtuple( _outdir+
					    "events_passing_cuts-on.ntuple",
					    true );

    _writer_off = new ParamDataWriterNtuple( _outdir+
					     "events_passing_cuts-off.ntuple",
					     true );
    

}

Cutter::~Cutter() {

    if (_energy_on != NULL) delete _energy_on;
    if (_energy_off != NULL) delete _energy_off;
    if (_energy_track != NULL) delete _energy_track;
    delete _alpha_on;
    delete _alpha_off;
    delete _size_on;
    delete _size_off;
    if (_xplotmaker) delete _xplotmaker;

    RawHeaderRecord header;
    _writer_on->writeHeader( header );
    _writer_off->writeHeader( header );

    delete _writer_on;
    delete _writer_off;
    
}


CutRecord
Cutter::
getTotal( vector<CutRecord> &vec ) {

    CutRecord total;

    for (int i=0; i<(int)vec.size(); i++) {
	total.addValuesFrom( vec[i] );
    }

    return total;

}
 
/**
 * Returns the maximum likelihood significance (see Li and Ma, 1983)
 * \param n_on number of on-source counts
 * \param n_off number of off-source counts
 * \param alpha 1/(tracking ratio) for tracking runs, or on/off duration 
 *        for on/off pairs.
 *
 * NOTE alpha is defined as 1/(tracking ratio) where the
 * tracking ratio is the ratio of (number of off-source bins)/(on-source bins)
 * So, e.g., for an "alpha" cut of alpha<15 for on-source and 20 < alpha < 65
 * for off-source, the tracking ratio = 3.0 and alpha = 0.333.
 * Defining alpha in this way, excess = n_on - alpha * n_off.
 */
double 
Cutter::
maxLikelihoodSignif( double n_on, double n_off, double alpha ) {
    double sum = n_on + n_off;
    double arg1 = (1.0+alpha)/alpha * ( n_on/sum);
    double arg2 = (1.0+alpha)*(n_off/sum);
    double sigsqr;

    // If more off than on, fall back to approximation to Li and Ma:
    if (alpha*n_off > n_on) {
	if (sum>1e-5) 
	    return (n_on-alpha*n_off)/sqrt(alpha*sum);
	else
	    return 0;
    }

    // Calculate Li and Ma signif:

    if (arg1 > 1e-10 && arg2 > 1e-10) {
	sigsqr = fabs(2.0*(n_on*log(arg1) + n_off*log(arg2)));
    }
    else {
	if ( sum > 1.0e-5 )
	    return (n_on-alpha*n_off)/sqrt(alpha*sum);
	else
	    return 0.0;
    }

    if ((n_on - alpha*n_off) > 0.0) 
	return sqrt(sigsqr);
    else
	return -sqrt(sigsqr);
    
}




RunStatistics 
Cutter::
getTrackingStatistics( CutRecord &run, double track_ratio ) {

    RunStatistics stats;
    
    stats.excess = run.orientation - run.trackoff/track_ratio;
    stats.significance = maxLikelihoodSignif( (double)run.orientation,
						 (double)run.trackoff,
						 1.0/track_ratio );

    return stats;

}


RunStatistics 
Cutter::
getPairStatistics( CutRecord &onrun, CutRecord &offrun ) {

    RunStatistics stats;
    double a = onrun.duration / offrun.duration;
    // TODO: maybe put option to use a =
    // onrun.trackoff/offrun.trackoff instead here to be used if there
    // is no matching off source data

    stats.excess = onrun.orientation - a*offrun.orientation;
    stats.significance = maxLikelihoodSignif( (double)onrun.orientation,
					      (double)offrun.orientation, a );
    return stats;

}


RunStatistics
Cutter::
getTotalPairStatistics() {

    RunStatistics prstats;
    CutRecord total_on,total_off;

    total_on  = getTotal( _onruns );
    total_off = getTotal( _offruns );

    // Generate 2d plots

    prstats = getPairStatistics( total_on, total_off );

    if (_is_2d_enabled) {
// 	cout << endl<<endl;
// 	cout << "Generating total 2D image..." <<endl;
	    generateImage( total_on, total_off, _outdir, 
			   _sourcenames.back(),true);
    }

    return prstats;

}


/**
 * Print out the totals for the currently processed runs.
 */
void
Cutter::
outputStatistics() {

    double excess, trexcess, alltrexcess;
    RunStatistics stats;

    Logger *logger = Logger::instance();

    if (_onruns.size() != _offruns.size()) {
	logger->printf("Number of ON runs does not match OFF runs!");
    }

    // Print out/save stars in the Field of View:
    outputStars();
    cout << endl;
    
    // Check for pointing offsets by cross-correlating the sky
    // brightness maps!
    
    if (_is_2d_enabled && _pointing_check_is_enabled) 
	checkPointing();

    // Accumulate total statistics

    CutRecord total_on, total_off, total_trk, total_all;

    total_on  = getTotal( _onruns );
    total_off = getTotal( _offruns );
    total_trk = getTotal( _trackruns );
    total_all.addValuesFrom(total_on);
    total_all.addValuesFrom(total_trk);

    // Per Run summary:

    cout << "Pair Analysis: " << endl;
    if (_onruns.size() >0) {
	cout << hline << endl;
	printCutRecordFields(cout);
	cout << hline << endl;
	for (int i=0; i<(int)_onruns.size(); i++) {
	    cout << _onruns[i]  << endl;
	    cout << _offruns[i] << endl;
	}
	//PFR: changed elevation to zenith
	cout << hline << endl;
	cout << "PAIR                       EXCESS         SIGNIF     <ZENITH>"<<endl;
	cout << hline << endl;
       	for (int i=0; i<(int)_onruns.size(); i++) {
	    stats = getPairStatistics( _onruns[i], _offruns[i] );
	    cout << _onruns[i].runid <<"/"<<_offruns[i].runid 
		 << setw(15-_onruns[i].runid.length()-_offruns[i].runid.length()) << " "
		 << setw(16) << stats.excess
		 << setw(16) << stats.significance 
		 << setw(16) << 90.0 - _onruns[i].average_elevation*180.0/M_PI
		 << endl;
	}

    }

    cout << endl<<endl;
    cout << "Tracking Analysis: "<< endl;

    if ( _trackruns.size() > 0 ) {
	cout << hline << endl;
	printCutRecordFields(cout);
	cout << hline << endl;
	for (int i=0; i<(int)_trackruns.size(); i++) {
	    cout << _trackruns[i] << endl;
	}

    }
    //PFR: changed elevation to zenith
    cout << hline << endl;
    cout << "RUN                        EXCESS         SIGNIF     <ZENITH>"
	 << endl;
    cout << hline << endl;
    for (int i=0; i<(int)_trackruns.size(); i++) {
	stats = getTrackingStatistics( _trackruns[i], _track_ratio );
	cout << _trackruns[i].runid << setw(16-_trackruns[i].runid.length())<< " "
	     << setw(16) << stats.excess 
	     << setw(16) << stats.significance 
	     << setw(16) << 90.0 - _trackruns[i].average_elevation*180.0/M_PI
	     << endl;
    }
    
    cout << "ON runs analyzed as tracking:"<<endl;
    for (int i=0; i<(int)_onruns.size(); i++) {
	stats = getTrackingStatistics( _onruns[i], _track_ratio );
	cout << _onruns[i].runid << setw(16-_onruns[i].runid.length())<< " "
	     << setw(16) << stats.excess 
	     << setw(16) << stats.significance 
	     << setw(16) << 90.0 - _onruns[i].average_elevation*180.0/M_PI
	     << endl;
    }
    
    checkDiagnostics();

    // print out summary

    cout << endl;
    cout << "Combined Analysis:" << endl;
    cout << "============================================="
	 << "==================================" <<  endl
	 << "(TOTALS)                    ON            OFF          TRACK   ALL-AS-TRACK"
	 << endl
	 << "============================================="
	 << "==================================" <<  endl;
    cout << "RUNS         : " 
	 << setw(15) << (int) _onruns.size()
	 << setw(15) << (int) _offruns.size()
	 << setw(15) << (int) _trackruns.size()
	 << setw(15) << (int) _trackruns.size() + (int)_onruns.size()
	 << endl;

    cout << "DURATION     : " 
	 << setw(15) << total_on.duration
	 << setw(15) << total_off.duration
	 << setw(15) << total_trk.duration
	 << setw(15) << total_all.duration
	 << endl;

    cout << "TOTAL EVENTS : " 
	 << setw(15) << total_on.total
	 << setw(15) << total_off.total
	 << setw(15) << total_trk.total 
	 << setw(15) << total_all.total 
	 << endl;

    cout << "VALID EVENTS : " 
	 << setw(15) << total_on.valid
	 << setw(15) << total_off.valid  
	 << setw(15) << total_trk.valid 
	 << setw(15) << total_all.valid 
	 << endl;

    cout << "PASS TRIG    : " 
	 << setw(15) << total_on.trigger  
	 << setw(15) << total_off.trigger  
	 << setw(15) << total_trk.trigger  
	 << setw(15) << total_all.trigger  
	 << endl;

    cout << "PASS SHAPE   : " 
	 << setw(15) << total_on.shape    
	 << setw(15) << total_off.shape    
	 << setw(15) << total_trk.shape     
	 << setw(15) << total_all.shape     
	 << endl;

    cout << "PASS ORIENT  : " 
	 << setw(15) << total_on.orientation      
	 << setw(15) << total_off.orientation
	 << setw(15) << total_trk.orientation       
	 << setw(15) << total_all.orientation       
	 << endl;

    cout << endl
	 <<"TRACK OFF    : " 
	 << setw(15) << total_on.trackoff      
	 << setw(15) << total_off.trackoff
	 << setw(15) << total_trk.trackoff
	 << setw(15) << total_all.trackoff
	 << endl;

    
    cout << endl
	 << "============================================="
	 << "==================================" <<  endl
	 << "(STATISTICS)            ON/OFF       TRACKING   ALL-AS-TRACK"
	 << endl
	 << "============================================="
	 << "==================================" <<  endl;

    RunStatistics allstats, prstats, trstats;

    allstats = getTrackingStatistics( total_all, _track_ratio );
    trstats = getTrackingStatistics( total_trk, _track_ratio );
    prstats = getPairStatistics( total_on, total_off );
  
    cout << "EXCESS       : "
	 << setw(15) << prstats.excess
	 << setw(15) << trstats.excess
	 << setw(15) << allstats.excess
	 <<endl;

    cout << "SIGNIFICANCE : "
	 << setw(15) << prstats.significance
	 << setw(15) << trstats.significance
	 << setw(15) << allstats.significance
	 <<endl;

    cout.precision(8);

    cout << endl
	 << "Tracking results used a tracking ratio of " 
	 << _track_ratio << endl 
	 << "On/Off results used a on/off duration ratio of " 
	 << total_on.duration/total_off.duration << endl
	 << endl;

    if (_offruns.size() > 0) {
	cout << "The derived tracking ratio from " << _offruns.size()
	     << " off runs is "
	     << (double)total_off.trackoff/(double)total_off.orientation 
	     << endl;
	cout << "'off-alpha' for ON runs : "<<total_on.trackoff << endl;
	cout << "'off-alpha' for OFF runs: "<<total_off.trackoff << endl;
	cout << "'off-alpha' ON/OFF ratio: "
	     <<(double)total_on.trackoff/(double)total_off.trackoff
	     <<endl;
       }


    _alpha_on->save(_outdir+"alpha-total-on");
    _alpha_off->save(_outdir+"alpha-total-off");
    _energy_on->save(_outdir+"energy-total-on");
    _energy_track->save(_outdir+"energy-total-track");
    _energy_off->save(_outdir+"energy-total-off");
    _size_on->save(_outdir +"size-total-on");
    _size_off->save(_outdir +"size-total-off");

    // Generate 2d plots

    if (_is_2d_enabled) {
	cout << endl<<endl;
	cout << "Generating total 2D image..." <<endl;
	if (_onruns.size()>0) 
	    generateImage( total_on, total_off, 
			   _outdir, _sourcenames.back(),true);

	if (_trackruns.size()>0) {

	    CutRecord fake;
	    
	    fake.duration = total_trk.duration;

	    generateImage( total_trk, fake, 
			   _outdir, 
			   _sourcenames.back()+" (tracking)",true);
	    
	}
    }


    // Sanity checks

    checkDiagnostics();

}

/**
 * Cross-correlate the skybrightness maps to check for pointing errors
 */
void
Cutter::
checkPointing() {

    int numruns = _onruns.size();
    Image2D *xcor;
    Image2D *source;
    int ix,iy;
    int ci[numruns][numruns];
    int cj[numruns][numruns];

    double xsum,ysum;
    int ixp, iyp;
    int n;
    double denom;
    double weight;
    char filename[128];

    if (numruns <= 0) return;

    _corr_centroid.resize(numruns);
    for (int i=0; i<numruns; i++)
	_corr_centroid[i].resize(numruns);
    

    cout << "Cross-correlating the skybrightness maps..." << endl;
    ProgressBar prog(((numruns-1)*(numruns-1))/2);
    int count=0;

    for (int i=0; i<numruns; i++){
	for(int j=0; j<numruns; j++){
	    _corr_centroid[i][j].x = 0;
	    _corr_centroid[i][j].y = 0;
	}
    }

    // Loop over all possible pairs of runs...

    for (int i=0; i<numruns-1; i++) {
	for (int j=i+1; j<numruns; j++) {

	    prog.print(count);
	    count++;
	    
	    source = &(_onruns[i].skybright);
	    xcor = source->crossCorrelateWith( _onruns[j].skybright );
	    xcor->getIndexOfMax( ix,iy );

	    if (ix<xcor->getXDim() && iy<xcor->getYDim()) {
		ci[i][j] = ix;
		cj[i][j] = iy;
	    }
	    else {
		cout << "ERROR: ix="<<ix<<" iy="<<iy<< endl;
		ci[i][j] = 0;
		cj[i][j] = 0;
	    }

	    // find the centroid about the max
	    xsum=0;
	    ysum=0;
	    n=0;
	    denom=0;
	    for (int k=-_centroid_range; k<=_centroid_range; k++) {
		for (int l=-_centroid_range; l<=_centroid_range; l++) {

		    ixp = (ix + k + xcor->getXDim() ) % xcor->getXDim();
		    iyp = (iy + l + xcor->getYDim() ) % xcor->getYDim();

		    weight = xcor->getPixel( ixp,iyp );

		    xsum += weight * xcor->getXCoord(ixp);
		    ysum += weight * xcor->getYCoord(iyp);
		    denom += weight;
  
		    n++;
	
		}
	    }

	    // store the centroid positions

	    _corr_centroid[i][j].x = xsum/denom;
	    _corr_centroid[i][j].y = ysum/denom;
	    _corr_centroid[j][i].x = xsum/denom;
	    _corr_centroid[j][i].y = ysum/denom;


	    if (i==0) {
		sprintf(filename, "%s/cross-corr-%s-%s", _outdir.c_str(),
			_onruns[i].runid.c_str(), _onruns[j].runid.c_str());
		xcor->save(filename);
	    }

	    delete xcor;
	}
    }

    prog.printClear();

    string fname = _outdir+"cross-corr-centers.dat";
    ofstream centfile( fname.c_str() );
    centfile << "Cross-correlation center pixels: "<< endl;
    for (int i=0; i<numruns; i++) {
    	for (int j=i+1; j<numruns; j++) {
	    centfile << "("<<ci[i][j]<<","<<cj[i][j]<<")\t";
	}
	centfile << endl;
    }

    centfile << "Cross-correlation center positions: "<< endl;
    for (int i=0; i<numruns; i++) {
    	for (int j=i+1; j<numruns; j++) {
	    centfile << "("<<_corr_centroid[i][j].x<<","
		 <<_corr_centroid[i][j].y<<")\t";
	}
	centfile << endl;
    }
    centfile.close();

    // figure out which runs seem mismatched and output warning
    // messages

    double offs;
    float badruns[numruns];
    int nbadruns[numruns];
    int bestrun;
    float minbad=99999;
    for(int i=0; i<numruns; i++) {
	badruns[i] = 0;
	nbadruns[i] =0;
    }

    for (int i=0; i<numruns; i++) {
    	for (int j=i+1; j<numruns; j++) {
	    offs = hypot(_corr_centroid[i][j].x, _corr_centroid[i][j].y);
	    if (offs > 0.1){
		badruns[i]+=offs;
		badruns[j]+=offs;
		nbadruns[i]++;
		nbadruns[j]++;
	    }
	}
    }

    for (int i=0; i<numruns; i++) {
	if (badruns[i] < minbad) {
	    minbad = badruns[i];
	    bestrun=i;
	}
    }

    cout << "Run with least pointing offset: "<<_onruns[bestrun].runid<<endl;
        
    cout << "Pointing offsets relative to "<<_onruns[bestrun].runid
	 << " and  mismatch level:" << endl;

    for (int i=0; i<numruns; i++) {
	
	offs = hypot(_corr_centroid[bestrun][i].x,
		     _corr_centroid[bestrun][i].y);

	cout << _onruns[i].runid << " "
	     << setw(10)<< offs
	     << setw(10)<< badruns[i]
	     << setw(10)<< nbadruns[i];
	if (offs > 0.1){
	    cout << " ** BAD?";
	}
	cout << endl;
    }

    // Output first row of cross correlations to a file
    fname = _outdir+"cross-correlation.dat";
    cout << "See "<< fname << " for detailed pointing offset data"<< endl;
    ofstream corrfile( fname.c_str() );
    corrfile << "Comparison with run "<<_onruns[bestrun].runid <<endl;
    for (int j=0; j<numruns; j++) {
	corrfile << "# "<<_onruns[j].runid  <<endl
		 <<"XAlign=" <<_corr_centroid[bestrun][j].x <<endl
		 <<"YAlign=" << _corr_centroid[bestrun][j].y <<endl
		 << endl;
    }
    corrfile.close();

    // output again in a machine readable format for alignment
    
    ofstream alignfile;
    
    for (int i=0; i<numruns; i++) {
	fname = _outdir+"alignment_to_"+_onruns[i].runid+".align";
	alignfile.open( fname.c_str() );
	for (int j=0; j<numruns; j++) {
	    alignfile << _onruns[j].runid << " "
		      <<_corr_centroid[i][j].x << " "
		      << _corr_centroid[i][j].y 
		      << endl;
	}
	alignfile.close();
    }

}

void
Cutter::
outputStars(){

    list<Star>::iterator star;
    string prevname;
    int oncount=0;
    int offcount=0;

    // First, lets sort the star list by name so we can deal with
    // duplicates easier.

    _onstars.sort();
    _offstars.sort();

    // Mark duplicates by setting their catalog_id to something > 0
    // (Stars with catalog_id==0 are plotted with names in the IDL
    // plotting script)

    for (star=_onstars.begin(); star!=_onstars.end(); star++) {
	if (star->name != prevname || star == _onstars.begin()) {
	    prevname = star->name;
	    oncount++;
	}
	else {
	    star->catalog_id++;
	}
    }

    for (star=_offstars.begin(); star!=_offstars.end(); star++) {
	if (star->name != prevname || star == _offstars.begin()) {
	    prevname = star->name;
	    offcount++;
	}
	else {
	    star->catalog_id++;
	}
    }


    // Now, write the star positions to a text file, for easy plotting

    writeStarList( _onstars, _outdir+"nearbystars-on.dat" );
    writeStarList( _offstars, _outdir+"nearbystars-off.dat" );


    // Remove duplicates and write star names to screen

    _onstars.unique();
    _offstars.unique();

    cout << "BRIGHT STARS IN ON-SOURCE FOV:"<< endl;
    for (star=_onstars.begin(); star!=_onstars.end(); star++) {
	if (star->magnitude < 6) {
	    cout << "\t'"<<star->name
		 << "', "<<"magnitude "<< star->magnitude
		 << ", "<<"at "<<setprecision(4)
		 << hypot(star->tx,star->ty)<<" deg from "
		 << "camera center"
		 << endl;
	}
    }

    cout << "BRIGHT STARS IN OFF-SOURCE FOV:"<< endl;
    for (star=_offstars.begin(); star!=_offstars.end(); star++) {
	if (star->magnitude < 6) {
	    cout << "\t'"<<star->name
		 << "', "<<"magnitude "<< star->magnitude
		 << ", "<<"at "<<setprecision(4)
		 << hypot(star->tx,star->ty)<<" deg from "
		 << "camera center"
		 << endl;
	}
    }

}

void 
Cutter::
writeStarList( list<Star> &starlist, string filename ) {

    list<Star>::iterator star;

    ofstream starfile( filename.c_str() );
    if (starfile.fail())
	throw MildAnalysisException("Couldn't write star file: "+filename);

    if (starlist.size()==0) return;

    starfile << starlist.size() << endl;

    for (star=starlist.begin(); star!=starlist.end(); star++) {

	starfile << *star << endl;

    }
    
    starfile.close();

}



/**
 * Save a datafile (image2d-total.dat) containing the 2D excess and
 * significance values for all of the data that has been currently
 * cut.
 *
 * The format of the output file is 6 columns per line:
 *   - column 1: x camera coordinate bin (in degrees)
 *   - column 2: y camera coordinate bin (in degrees)
 *   - column 3: excess counts 
 *   - column 4: significance
 *   - column 5: Right Ascension of bin
 *   - column 6: Declination of bin
 */
void
Cutter::
generateImage(CutRecord &on, CutRecord &off, string dir, string title,
	      bool finalimage){

    int i,j;
    double x,y;
    Image2D excess2d(_image_xdim,_image_ydim);
    Image2D signif2d(_image_xdim,_image_ydim);
    Image2D ra2d(_image_xdim,_image_ydim);
    Image2D dec2d(_image_xdim,_image_ydim);

    double minx,maxx,miny,maxy;
    list<Star>::iterator star;
    Image2D onimage = on.im2d;   // store on and off images so we're
    Image2D offimage = off.im2d; // no smoothing the original
    double a = on.duration / off.duration; 
    double non, noff;

    char telid[100];
    sprintf( telid, "%d", _telescope_id );



    // set up xplotmaker if it is not already there

    if (_display_is_enabled) {
	if (_xplotmaker == NULL) {
	    _xplotmaker = new PlotMaker( "X", "x.out" );
	}
    }

    // set initial values
    
    minx = -_image_xdim*_degperpix/2.0;
    maxx = _image_xdim*_degperpix/2.0;
    miny = -_image_ydim*_degperpix/2.0;
    maxy = _image_ydim*_degperpix/2.0;
 
    excess2d.setCoordinateBox( minx,miny  , maxx,maxy );
    signif2d.setCoordinateBox( minx,miny  , maxx,maxy );
    ra2d.setCoordinateBox( minx,miny  , maxx,maxy );
    dec2d.setCoordinateBox( minx,miny  , maxx,maxy );

    // If fast processing isn't enabled, rebin all of the points of
    // origin in the poolist vector using the more accurate
    // Image::addHistRadially() function.  Otherwise, just add up the
    // binned images and smooth afterward. The second method is much
    // faster, but works poorly when the smoothing radius is on the
    // order of the image's pixel spacing.

    if (_fast_image_processing_is_enabled == false) {

	Image2D tempon(_image_xdim,_image_ydim);  
	Image2D tempoff(_image_xdim,_image_ydim);  

	tempon.setCoordinateBox( minx,miny  , maxx,maxy );
	tempoff.setCoordinateBox( minx,miny  , maxx,maxy );
	vector<Coordinate_t>::iterator coord;
	double dist;

	for (coord=on.poolist.begin();coord!=on.poolist.end();coord++) {
	    tempon.addHistRadially(coord->x,coord->y,1.0,
				   _radial_smoothing_factor);
	}
	for (coord=off.poolist.begin();coord!=off.poolist.end();coord++) {
	    tempoff.addHistRadially(coord->x,coord->y,1.0,
				    _radial_smoothing_factor);
	}
	
	onimage = tempon;
	offimage=tempoff;
    }
    else {
	// do it the fast, but less accurate way
	onimage.applyRadialSmoothing( _radial_smoothing_factor );
	offimage.applyRadialSmoothing( _radial_smoothing_factor );
    }

    // Do smoothing
    
//     onimage.applyRadialSmoothing( _radial_smoothing_factor );
//     offimage.applyRadialSmoothing( _radial_smoothing_factor );

    // Calculate significance and excess:

    for (i=0; i<_image_xdim; i++) {
	for (j=0; j<_image_ydim; j++) {

	    non = onimage.getPixel(i,j);
	    noff = offimage.getPixel(i,j);

	    // Calculate the excess
	    excess2d.setPixel( i,j, non-noff*a);

	    // calculate significance for each point 
	    signif2d.setPixel(i,j, maxLikelihoodSignif(non,noff,a) );

	}
    }


    // Calculate the RA/DEC coordinates of each gridpoint.  This
    // allows simple contour plotting to show ra/dec lines:

    for (i=0; i<_image_xdim; i++) {
	for (j=0; j<_image_ydim; j++) {

	    x = ra2d.getXCoord(i) * M_PI/180.0; // these must be in radians
	    y = ra2d.getYCoord(j) * M_PI/180.0;
	    
	    dec2d.setPixel(i,j, (atan( ((sin(_dec) + y*cos(_dec))/
				       (cos(_dec) - y*sin(_dec))) *
				      cos( atan(-1.0*x/
						(cos(_dec)-y*sin(_dec)))) ))
			   * 180.0/M_PI );
			   
	    ra2d.setPixel(i,j, (_ra + atan(x/(y*sin(_dec)-cos(_dec))))*
			  12.0/M_PI );

	}
    }

    // Output the 2d significance and excess:

    ra2d.save( dir+"ra-mesh" );
    dec2d.save( dir+"dec-mesh" );

    excess2d.save( dir+"excess");
    signif2d.save( dir+"signif");

//     cout << "DEBUG: "<<dir<<": 2d CENTER EXC,SIG: "
// 	 <<excess2d.getPixel(19,19)<<" "
// 	 <<signif2d.getPixel(19,19) << endl;

    // Plot the final image to screen and a postscript file:

    excess2d.expand(); // make hi-res
    signif2d.expand(); 

    vector<PlotMaker*> pm;
    
    if (_display_is_enabled) {
	pm.resize(2);
	if (!finalimage) {
	    pm[0] = _xplotmaker;
	    pm[1] = new PlotMaker( "PS", dir+"image.eps" );
	}
	else {
	    pm[0] = new PlotMaker( "X", "x.out" );
	    pm[1] = new PlotMaker( "PS", dir+"image.eps" );
	}
    }
    else {
	pm.resize(1);
	pm[0] = new PlotMaker( "PS", dir+"image.eps" );
    }
    
    
    for (i=0; i<pm.size(); i++) {

	pm[i]->setAxisBox( excess2d );
	pm[i]->plot( excess2d );
	pm[i]->setContourColorName( "black" );
	pm[i]->setContourLineStyle( PlotMaker::LINE_SOLID );
	pm[i]->contourPlot( signif2d, 1.0, 50.0, 1.0 );
	pm[i]->setColorMap( PlotMaker::COLORMAP_BLUES );
	pm[i]->setContourColorName( "orange" );
	pm[i]->setContourLineStyle( PlotMaker::LINE_DASHED );
	pm[i]->plotRADecGrid( ra2d, dec2d );
	pm[i]->plotAxes();
	pm[i]->plotTitle( "Significance & Excess");
	pm[i]->plotSubtitle( title+" (T"+telid+")" );

	// plot the stars:
	for (star=_onstars.begin(); star!=_onstars.end(); star++) {
	    pm[i]->plot( *star );
	}
	for (star=_offstars.begin(); star!=_offstars.end(); star++) {
	    pm[i]->plot( *star, 4 );
	}

	pm[i]->flush();

	if (pm[i] != _xplotmaker) 
	    delete pm[i];
	else 
	    pm[i]->newPage();
    }


    // Plot the brightness Maps need to generate new copies of them,
    // so the old ones aren't overwritten

    Image2D skybrightness = on.skybright;
    Image2D tubeoffness = on.tubeoff;
    tubeoffness.addImage( off.tubeoff );


    if (finalimage==false || _display_is_enabled==false) {
 	pm.resize(1);
       	pm[0] = new PlotMaker( "PS", dir+"skybrightness.eps" );
    }
    else {
	pm.resize(2);
	pm[0] = new PlotMaker( "X", ".y.out" );
	pm[1] = new PlotMaker( "PS", dir+"skybrightness.eps" );
    }

    skybrightness.expand();     // Make higher resolution and interpolate
    tubeoffness.expand();

    for (i=0; i<pm.size(); i++) {
	pm[i]->setColorMap( PlotMaker::COLORMAP_INVERSEGREY );
	pm[i]->setAxisBox( skybrightness );
	pm[i]->plot( skybrightness );
	pm[i]->setContourColorName( "LightGreen" );
	pm[i]->contourPlot( tubeoffness,
			    tubeoffness.minValue(),
			    tubeoffness.maxValue(),
			    (tubeoffness.maxValue()-
			    tubeoffness.minValue())/10.0 );
	pm[i]->setContourColorName( "orange" );
	pm[i]->setContourLineStyle( PlotMaker::LINE_DASHED );
	pm[i]->plotRADecGrid( ra2d, dec2d );
	pm[i]->plotAxes();
	pm[i]->plotTitle( "Sky Brightness"  );
	pm[i]->plotSubtitle( title+" (T"+telid+")" );

	for (star=_onstars.begin(); star!=_onstars.end(); star++) {
	    pm[i]->plot( *star );
	}
	for (star=_offstars.begin(); star!=_offstars.end(); star++) {
	    pm[i]->plot( *star, 4 );
	}

	pm[i]->flush();
	
	if (pm[i] != _xplotmaker) delete pm[i];
	else  pm[i]->newPage();
    }

}

/**
 * Cut the data in the specified run, and build up statistics for all
 * runs. Each run's cut information is stored in the
 * onruns/offruns/trackruns vectors.
 */ 
void 
Cutter::
process( RunInfo &ri ) {
    
    Logger *logger = Logger::instance();
    
    _radial_smoothing_factor = ri.cuts.smoothing_radius;
    
    // do some checks

    if (ri.cuts.radial_analysis && _is_2d_enabled==false) {
	logger->printf("You requested radial analysis, "
		       "but 2D analysis is disabled! Enabling 2D...");
	enable2D(true);
    }

    // now process the run.

    try {
	if (ri.type == Config::ONOFF) {

	    _onruns.push_back(cut( ri, ri.onid, ON ));
	    _offruns.push_back(cut( ri, ri.offid, OFF ));

	    if (_is_2d_enabled) {

		generateImage( _onruns[_onruns.size()-1],
			       _offruns[_offruns.size()-1],
			       ri.onid+"/",ri.onid+"/"+ri.offid );
	    }

	}
	else {
	    
	    _trackruns.push_back(cut( ri, ri.onid, TRACK ));

	    if (_is_2d_enabled) {
		
		// make a fake off run so that generateImage works

		CutRecord fake;
		fake.duration = _trackruns[_trackruns.size()-1].duration;

 		generateImage( _trackruns[_trackruns.size()-1],
 			       fake,
 			       ri.onid+"/",ri.onid+" (tracking)" );
	    }


	}

    } 
    catch (MildAnalysisException &e) {
	cout << "*** WARNING: pair ";
	if ( ri.type == Config::ONOFF) 
	    cout << ri.onid <<"/"<<ri.offid<<": ";
	else 
	    cout << ri.onid;

	cout << e.what() << endl;
	
    }

    
}

/**
 * Read the parameterized data for the specified run id and count the
 * number of events that pass the various cuts. Also generates
 * histograms of rate, size, energy, etc.  This function is called by
 * Cutter::process for each run in the configuration file, so there is
 * usually no need to call it directly.
 *
 * \returns a CutRecord containing all the totals for the specified run.
 *
 * \todo: Implement 2-D analysis for tracking runs (currently, only
 * ON/OFF pairs are used.
 *
 * \todo: separate total alpha plot for tracking runs.
 *
 */
CutRecord
Cutter::
cut( RunInfo &ri, const string &id, char type ) {

    HillasParameterization param,oldparam;
    ParamDataReader *reader;
    HeaderRecord header;
    string efilename = id+"/"+id+"-energy.dat";
    ofstream energyfile(efilename.c_str());

    double x;

    Image2D centroids(39,39);
    centroids.setCoordinateBox( -_image_xdim*_degperpix/2.0, 
			   -_image_ydim*_degperpix/2.0,
			   _image_xdim*_degperpix/2.0,
			   _image_ydim*_degperpix/2.0 );    

    EnergyEstimatorFactory *efactory = EnergyEstimatorFactory::instance();
    EnergyEstimator *energy = efactory->getEstimator(ri.energyestimator);

    // allocate the energy spectrum histograms if needed:

    if (_energy_on == NULL) 
	_energy_on = new EnergySpectrum(ri,"EnergyEstimator-on");
    if (_energy_off == NULL) 
	_energy_off= new EnergySpectrum(ri,"EnergyEstimator-off");
    if (_energy_track == NULL) 
	_energy_track = new EnergySpectrum(ri,"EnergyEstimator-track");
      
    // write out cuts that are to be used:

    string cutinfofile = id+"/"+id+"-cuts.txt";
    ofstream cutinfo(cutinfofile.c_str());
    if (!cutinfo.fail()) {
	cutinfo << ri.cuts << endl;
	cutinfo.close();
    }
    else {
	throw AnalysisException("Couldn't write to "+cutinfofile);
    }

    // open data

    reader = new ParamDataReader(id+"/"+id+"-param.ntuple");
    reader->getHeaderRecord( header );

    // throw error if the user specified a telescope which is not
    // specified in the header:
    if (_telescope_id >= header.num_telescopes) {
	cout << "NOTE: you specified telescope number "<<_telescope_id
	     << " but there are only " << header.num_telescopes <<endl
	     << "      telescopes in the array (the numbers start at 0)!" 
	     << endl;
	delete reader;
	delete energy;
	throw CriticalAnalysisException( "Invalid telescope number");
    }

    
    if (type == ON || type == TRACK) {
	_ra = header.ra;
	_dec = header.dec;
    }

    // add source name to list for sanity check:
    _sourcenames.push_back( header.sourcename );

    if (type == ON) {
	_on_ra.push_back(header.ra);
	_on_dec.push_back(header.dec);
    }
    else if (type == OFF) {
	_off_ra.push_back(header.ra);
	_off_dec.push_back(header.dec);
    }


    _track_ratio = ri.cuts.tracking_ratio;

    //============================================================
    // Initialize Histograms

    Histogram relrate( 30, 0,30, "Rate" );
    Histogram alpha( 18, 0,90, "alpha" );
    Histogram sizehist( 30,0,1000, "Size" );
    Histogram ehist( 18,-1.0,2.0, "EnergyEstimator");
    
    //============================================================
    // Determine the lightcurve binning start and end times based 
    // on the preset timebase.
    
    const double MJD_TO_MIN = 24*60;
    double tstartmin = header.starttime * MJD_TO_MIN;
    double tendmin   = header.endtime * MJD_TO_MIN;
    double tbasemin  = ri.utbase * MJD_TO_MIN;
    double diff1     = tstartmin - tendmin;
    double delta1    = diff1 - floor(diff1);
    double tbinstart = tstartmin - delta1;
    double diff2     = tendmin - tbinstart;
    double delta2    = ceil(diff2) - diff2;
    double tbinend   = tendmin + delta2;

    int nminutes = (int) (tbinend - tbinstart);

    if (nminutes<1) nminutes = 1;

    tbinstart /= MJD_TO_MIN;
    tbinend  /= MJD_TO_MIN;

    if (tbinstart>tbinend) {
	double tmp = tbinstart;
	tbinstart = tbinend;
	tbinend = tmp;
    }

    // Now tbinstart and tbinend are the start end end times
    // normalized to the base time (ri.utbase) assuming 1 minute bins.
    // The first and last bins are going to be partially filled, by
    // (1-delta1) and (1-delta2) minutes respectively.
    
    Histogram rate( nminutes, tbinstart, tbinend, "Absolute Rate" );
    Histogram toffrate( nminutes, tbinstart, tbinend, 
			      "TrackOff Rate" );
//     rate.setBinFillFraction( 0, 1.0-delta1 );
//     rate.setBinFillFraction( nminutes-1, 1.0-delta2 );
//     toffrate.setBinFillFraction( 0, 1.0-delta1 );
//     toffrate.setBinFillFraction( nminutes-1, 1.0-delta2 );
//     toffrate.setScaleFactor( 1.0/ri.cuts.tracking_ratio );


    //============================================================
    // Loop over events and cut the data

    bool first = true;
    double elapsed  = 0;
    double starttime = 0;
    bool pass_radial;
    double dista,distb;
    int n_too_close=0;

    CutRecord thisrun;		// initialized to 0 in constructor
    thisrun.runid=id;	  
    thisrun.average_elevation = header.average_elevation;

    // set up the requested cut method:

    SuperCutter *cutter;
    CutFactory *cutfactory = CutFactory::instance();
    cutter = cutfactory->newCutter( ri );
    
    cutter->setCuts( ri.cuts );

    if (ri.cuttype == Config::ZCUTTER || ri.cuttype == Config::EZCUTTER) {
	((ZCutter*)(cutter))->setCamera(header.nadc[_telescope_id],atoi(ri.utdate.c_str()));
    }

    // loop over data

    for (int i=0; i<reader->size(); i++) {

	reader->getNextEventRecord( param );

	// skip events which aren't from the specified telescope
	// (wucut is not designed for array analysis)
	if (param.telescope_id != _telescope_id) continue;

	// Sotre info about the first event for time calculations later...
	if (first==true) {
	    starttime = param.gpstime;
	    first = false;
	}

	elapsed = param.gpstime-starttime;

	thisrun.total++;
	thisrun.valid++;

	// Apply elongation and other corrections to the parameters
	// (specific to the type of cutter used). Store original
	// uncorrected parameters for later use. 

	oldparam=param; 
	cutter->applyCorrections( param );

	// shift points of origin for alignment if requested.  This is
	// used for stacking a bunch of runs that may have pointing
	// errors, and of course only works with 2D analysis (alpha
	// cut is unaffected)

	if (ri.cuts.alignment == true) {
	    
	    param.point_of_origin_a.x -= ri.cuts.align_offset.x;
	    param.point_of_origin_b.x -= ri.cuts.align_offset.x;

	    param.point_of_origin_a.y -= ri.cuts.align_offset.y;
	    param.point_of_origin_b.y -= ri.cuts.align_offset.y;
	}
	
	// ----------------------------------------------------
	// Trigger and Muon Cut:
	
	if (cutter->pass( param, SIZE|MAX|FRAC|LENSIZE )) {
	    
	    thisrun.trigger++;

	    // ----------------------------------------------------
	    // Shape Cut:
	    if (cutter->pass( param, DISTANCE|LENGTH|WIDTH )) {
		
		thisrun.shape++;
		
		// Update "alpha total" plot
		
		if (type == ON || type == TRACK)
		    _alpha_on->increment(param.alpha*180.0/M_PI);
		else if (type == OFF)
		    _alpha_off->increment(param.alpha*180.0/M_PI);
		else 
		    cout << "** Cutter: unknown run type " << type << endl;
		
		alpha.increment(param.alpha*180.0/M_PI);
		
		// For tracking analysis, count events in the
		// "off" alpha region
		
		if (param.alpha*180.0/M_PI >= 20.0 &&
		    param.alpha*180.0/M_PI <= 65.0) {
		    thisrun.trackoff++;
		    toffrate.increment( param.gpstime );
		}
		
		
		// update 2D analysis...
		
		if (_is_2d_enabled){
		    

		    // Accept the preferred point of origin "a". It's
		    // preferred in the sense that the asymmetry seems
		    // to point that way (see HillasImageAnalyzer)

		    if (_fast_image_processing_is_enabled) {
			thisrun.im2d.addHist(param.point_of_origin_a.x, 
					     param.point_of_origin_a.y, 1);
		    }
		    else {
			thisrun.poolist.push_back(param.point_of_origin_a);
		    }

		    // if it's not within the valid asymmetry region,
		    // then accept the second point of origin too
		    // since we can't trust asymmetry

		    if ( cutter->pass(param,ASYMM) == false) {

			// Since this method could lead to double
			// counting, we need to check if the points of
			// origin are too close to one another:
			
			if ( hypot(param.point_of_origin_a.x-
				   param.point_of_origin_b.x,
				   param.point_of_origin_a.y-
				   param.point_of_origin_b.y) 
			     > 2.0*ri.cuts.smoothing_radius) {

			    if (_fast_image_processing_is_enabled) {
				thisrun.im2d.addHist(param.point_of_origin_b.x,
						     param.point_of_origin_b.y,
						     1);
			    }
			    else {
				thisrun.poolist.push_back(param.point_of_origin_b);
			    }
			}
			else {
			    n_too_close++;
			}


		    }

		    // If an radial-analysis is enabled, see if
		    // distance beteween point of origin and offset
		    // position is < smoothing_radius (actually
		    // functioning as a radius cut, not an alpha cut).
		    // otherwise, do standard alpha cut

		    pass_radial = false;

		    if (ri.cuts.radial_analysis) {

			dista = hypot( param.point_of_origin_a.x 
				      - ri.cuts.radial_offset.x, 
				      param.point_of_origin_a.y 
				      - ri.cuts.radial_offset.y);

			if (cutter->pass(param,ASYMM) == true) {

			    // asymmetry is valid, so just look at A

			    if (dista <= ri.cuts.smoothing_radius) 
			    pass_radial = true;
			}
			else {

			    // otherwise, check both A and B

			    distb = hypot( param.point_of_origin_b.x 
					   - ri.cuts.radial_offset.x, 
					   param.point_of_origin_b.y 
					   - ri.cuts.radial_offset.y);
			    
			    if ((dista <= ri.cuts.smoothing_radius) || 
				(distb <= ri.cuts.smoothing_radius)) 
				pass_radial = true;
			}
			
		    }

		    
		}
		
		// ----------------------------------------------------
		// final Orientation cut.  This is either an alpha cut
		// or a radius cut depending on whether "offset
		// analysis" is enabled.
		
		if ( (ri.cuts.radial_analysis==true && pass_radial) || 
		     (ri.cuts.radial_analysis==false && 
		      cutter->pass( param, ALPHA )  )) {
		    
		    thisrun.orientation++;

		    // update the energy estimator histograms

		    // WARNING: tracking runs may not work for
		    // spectral analysis, but they're implemented
		    // anyway!

		    // NOTE: When we decided to use the older ZCuts
		    // over EZCuts for the SgrA spectrum, we realized
		    // that we had derived the energy estimator
		    // function using EZCuts (which modifies size but
		    // not distance), while ZCuts modifies both.  To
		    // allow the energyestimator to work in both
		    // cases, I pass oldparam.distance (the unscaled
		    // distance) to the estimator instead of
		    // param.distance!

 		    x = energy->getEstimate(param.size,oldparam.distance,
					    param.zenith);

		    // x is log10(Energy)

// 		    cout << "DEBUG: E= "<<param.sim.primary_energy
// 			 <<" Est= "<<pow(10,x)
// 			 <<" error= "
// 			 <<(log10(param.sim.primary_energy) - x)
// 			 <<endl;
		    
		    // store the energy estimate in case output is
		    // enabled (then it can be plotted against true
		    // energy for spectral resolution plots)
		    
		    param.energy_estimate = 
			oldparam.energy_estimate = pow(10,x);

		    if (type == ON)
			_energy_on->increment( x );
		    if (type == TRACK)
			_energy_track->increment( x );
		    else if (type == OFF)
			_energy_off->increment( x );

		    // update the diagnostic histograms

		    relrate.increment( elapsed*24*60 );
		    rate.increment( param.gpstime );
		    sizehist.increment( param.size );

		    if (type==ON) 
			_size_on->increment( log10(param.size) );
		    else if (type == OFF)
			_size_off->increment( log10(param.size) );

		    energyfile << setw(20) << x
			       << setw(20) << param.size
			       << setw(20) << param.distance
			       << endl;
		    ehist.increment( x );
		    centroids.addHist( param.centroid.x, 
				       param.centroid.y,1);
		    
		    // write out the parameters (if requested).
		    // Writes out the unscaled parameters by default,
		    // or the parameters scaled by the particular cut
		    // technique if the
		    // enableCorrectedParameterOutput(true) function
		    // has been called.
		    if (_output_is_enabled) {
			if (type==ON || type==TRACK){
			    if (_output_scaled_parameters)
				_writer_on->writeParameterization( param );
			    else
				_writer_on->writeParameterization( oldparam );
			}
			else if (type==OFF) {
			    if (_output_scaled_parameters) 
				_writer_off->writeParameterization( param );
			    else
				_writer_off->writeParameterization( oldparam );
			}
		    }
		    
		    
		}
	    }
   	}
	
    } // end loop over events

    
    //============================================================
    // Run duration is total elapsed time 

    thisrun.duration = elapsed*24*60;

    //============================================================
    // Update the sky brightness and tubeoffness maps

    if ((type==ON || type==TRACK)) {
	try {
	    
	    char telid[100];
	    sprintf( telid, "-%03d", _telescope_id );
	    thisrun.skybright.load(id+"/"+id+telid+"-sky_brightness.im2d");
	    thisrun.tubeoff.load(id+"/"+id+telid+"-tubeoffness.im2d");
	}
	catch (MildAnalysisException &e) {
	    
	    // clear the image if file doesn't exist (probably a sim run)
	    thisrun.skybright.clear();
	    thisrun.tubeoff.clear();

	}
    }

    //============================================================
    // Save histograms
    relrate.save( id+"/"+id+"-relative-rate" );
    sizehist.save( id+"/"+id+"-size" );
    rate.save( id+"/"+id+"-rate" );
    toffrate.save( id+"/"+id+"-rate-trackoff" );
    alpha.save( id+"/"+id+"-alpha" );
    ehist.save( id+"/"+id+"-energy" );

    // Write out centroids plot:

    centroids.save(id+"/"+id+"-cut-centroids");
    PlotMaker pm( "PS", id+"/"+id+"-cut-centroids.eps" );
    pm.setAxisBox( centroids );
    pm.setColorMap( PlotMaker::COLORMAP_INVERSEGREY );
    pm.plot( centroids );
    pm.plotAxes();
    pm.plotTitle( "Centroid Positions for " + id );
    pm.flush();

    //============================================================
    // Locate and save field stars


// KPK removed precession because it's not doing the right thing (the
// algorithm is fine, but the RA and DEC are already in j2000 epoch,
// which is what the stars are in so precessing doesn't make sense).
// Only need to precess if the star catalog is different from j2000.
// So really, should just precess to j2000 here always!

    _starcatalog.precessToDate( atoi(ri.utdate.c_str()) );
    
    if (type == ON || type == TRACK) {
	_onstars.clear();
	_starcatalog.findNearbyStars( header.ra, header.dec, 2.0, _onstars );
    }
   
    else{
	_offstars.clear();
	_starcatalog.findNearbyStars( header.ra, header.dec, 2.0, _offstars );
    }
    energyfile.close();
    delete reader;
    delete cutter;
    delete energy;
    return thisrun;

}


/**
 * Return the energy estimate for the given size and distance values.
 * \returns x = log(Energy)
 */
double 
Cutter::
energyEstimator( double size, double dist ) { 

    double x;
    double lnsize = log(size);
    static const double A = -8.11;
    static const double B = 2.56;
    static const double C = -9.25;
    static const double D = 0.120;
    static const double E = 6.26;
    static const double F = 0.0105;
    
    x  = (A + B*lnsize + C*dist - D*lnsize*lnsize +
	  E*dist*dist + F*lnsize*dist);

    return x;

}


/**
 * \returns zenith angle dependent elongation factor for the camera
 */
double
Cutter::getElongationFactor( double zen ) {

    const double A=0.912865;
    const double B=0.717872;
    
    return A*cos(zen) + B;

    // old constant value was 1.65    
    //     return 1.65;
   
}


/**
 * Run some diagnostics to look for bad things...
 */
void
Cutter::checkDiagnostics() {

    Logger *logger = Logger::instance();

    // As sanity check, make sure there aren't multiple source 
    // names listed!

    _sourcenames.sort();
    _sourcenames.unique();

    if (_sourcenames.size() > 1 ) {
	logger->printf("You may have combined multiple sources!");
	logger->printf("         the following were used:");
		list<string>::iterator it;
	for (it = _sourcenames.begin(); it!=_sourcenames.end(); it++){
	    logger->printf("\t%s", it->c_str() );
	}
    }
    
    // Look at the OFF-ALPHA values for on and off runs in a pair and
    // check that they are somewhat close to eachother, otherwise warn
    // the user that there may be a mismatch

    double percentdiff;

    for (int i=0; i<_onruns.size(); i++) {
	
	percentdiff = std::fabs((_onruns[i].trackoff - _offruns[i].trackoff) 
	    / (double)(_onruns[i].trackoff));
	
	if (percentdiff > 0.2) {
	    logger->printf("%s has an OFF-ALPHA value %2.0f%% different from "
			   "%s - check they are a valid pair",
			   _onruns[i].runid.c_str(),
			   percentdiff*100.0,
			   _offruns[i].runid.c_str()); 
	}
	
    }


    // check that all of the RAs and DECs make sense:

    for (int i=0; i<_on_ra.size(); i++) {
	double delta = fabs(_on_ra[i] - _off_ra[i]);
	
	delta *= 12.0/M_PI*60.0;
	
	if ( fabs(delta-30.0)>2.0 ) {
	    logger->printf("%s and %s have RAs which are different by %.1dmin! (should be 30min)",
			   _onruns[i].runid.c_str(),
			   _offruns[i].runid.c_str(),
			   delta );
	}
	
    }


}


ostream& operator<<( ostream &stream, CutRecord &c ) {

    stream << setw(15-c.runid.size()) << c.runid << ": "  
	   << setw(10) << c.total
	   << setw(10) << c.valid
	   << setw(10) << c.trigger
	   << setw(10) << c.shape
	   << setw(10) << c.orientation
	   << setw(10) << c.trackoff;

    return stream;

}

void
Cutter::printCutRecordFields( ostream &stream ) {
    stream << "ID"  << setw(8) << " "
	   << setw(10) << "TOTAL"
	   << setw(10) << " VALID"
	   << setw(10) << " TRIGGER"
	   << setw(10) << " SHAPE"
	   << setw(10) << " ORIENT" 
	   << setw(10) << " OFFALPHA" 
	   << endl;
}
