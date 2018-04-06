//
//
// Parameterizer.cpp
// Karl Kosack <kosack@hbar.wustl.edu>
//
// Modified 030201 by JB to correct some bugs in 2-d analys
// save some computation time and correct a bug where zenith
// corrections were applied twice when display was enabled
// 
//
// Modified July 2005 by KPK - Adding functionality to work with
// multiple telescopes (i.e. VERITAS).  Since the original code really wan't designed to do that, it's somewhat of a messy situation

#include <iostream>
#include <cmath>
#include <iomanip>
#include <list>
#include <cctype>
#include <sys/types.h>         	// For POSIX mkdir
#include <sys/stat.h>		// For POSIX mkdir
#include <gsl/gsl_rng.h>	// Random num generators 
#include <gsl/gsl_randist.h>	// Random distributions 
#include <gsl/gsl_sf_trig.h>

#include "Camera.h"
#include "DataReader.h"
#include "DataRecords.h"
#include "DataWriter.h"
#include "ImageAnalyzer.h"
#include "MuonImageAnalyzer.h"
#include "ImageCleaner.h"
#include "PedestalFinder.h"
#include "GainFinder.h"
#include "ProgressBar.h"
#include "Parameterizer.h"
#include "Histogram.h"
#include "GaussianDeviate.h"
#include "SuperCutter.h"
#include "ZCutter.h"
#include "SpectralCutter.h"
#include "Image2D.h"
#include "PlotMaker.h"
#include "Log.h"
#include "AngleConverters.h"
#include "ImageFilter.h"

using namespace std;
using namespace Cuts;

Parameterizer::
Parameterizer()
    : _is_verbose_enabled(true), _is_display_enabled(false),
      _is_interactive(false),_is_first_run(true) , _is_cleanup_enabled(true),
      _latitude(.553065751136), _longitude(1.935190) , _event_number(0) {

    _wupdir = getSupportDir();

}

Parameterizer::
~Parameterizer() {

    delete _gaussdev;
    
}

/**
 * Parametrize the specified run.
 */
void
Parameterizer::
process(  RunInfo &ri ) {

    PedestalFinder ped;
    RunInfo ri2 = ri;

    cout << endl;

    if (_is_first_run) {
	_gaussdev = new GaussianDeviate(ri.cachedir);
    }

    switch ( ri.type ) {

      case Config::ONOFF:

	  if (_is_verbose_enabled) {
	      cout <<"----------------------------------"
		   <<"--------------------------------------"<< endl;
	      cout <<"---- ON/OFF PAIR ---- " <<ri.onid<<"/"<<ri.offid<<" ";
	      for(int k=0; k< 70-(int)(ri.onid.size()
				       +(int)ri.offid.size()+1+21); k++)
		  cout << "-";
	      cout << endl;
	      cout <<"----------------------------------"
		   <<"--------------------------------------"<< endl;
	  }

	  ri2 = ri;

	  // Process on padded with off note padpeds is really
	  // off-source pedestals and must be loaded even if padding
	  // is off (for the skybrightness map)

	  //	  ped.getOffSourcePeds( ri, _padpeds );
	  processOne( ri, ri.onid );
	  
	  // Process off padded with on
	  //	  ped.getOnSourcePeds( ri, _padpeds );
	  processOne( ri2, ri2.offid );
	  break;
	  
      case Config::TRACK:
	  if (_is_verbose_enabled) {
	      cout <<"----------------------------------"
		   <<"--------------------------------------"<< endl;
	      cout <<"---- TRACKING RUN ---- " <<ri.onid<<" ";
	      for(int k=0; k< 70-(int)(ri.onid.size()+21+1); (int)k++)
		  cout << "-";
	      cout << endl;
	      cout <<"----------------------------------"
		   <<"--------------------------------------"<< endl;
	  }
	  
	  processOne( ri, ri.onid );
	  break;
	  
      default:
	  cerr << "*** Parameterizer: unknown run type: " << ri.type<<endl;
	  break;
    }
    
    _is_first_run = false;
    

}


/**
 * Print out run summary information in human readable format.
 */ 
void
Parameterizer::printSummary( ostream &stream, RunInfo &ri, 
			     string &id, RawHeaderRecord &header,
			     ParamDataWriter *writer, RawDataReader *data, 
			     TelescopeArray &array ){

    stream <<"Run information: " << endl;
    stream << "\tUTDATE      : " << ri.utdate << endl;
    stream << "\tRUN ID      : " << id << endl;
    stream << "\tNITROGEN ID : " << ri.n2id << endl;
    stream << "\tPADDING     : ";
    if (ri.padding) {
	if (ri.type == Config::ONOFF) {
	    if (id == ri.onid)
		stream << ri.offid << endl;
	    else
		stream << ri.onid << endl;
	}
	else {
	    stream << ri.padlevel << " (fixed level)"<< endl;
	}
    }
    else 
	stream << "disabled" << endl;
    
    stream << "\tINPUT FILE  : " << data->getFilename() << endl;
    stream << "\tINPUT TYPE  : " << data->getTypeString()<< endl;
    stream << "\tOUTPUT      : " << writer->getFilename() << endl;
    stream << "\tOUTPUT TYPE : " << writer->getTypeString() << endl;
    stream << "\tSOURCENAME  : \"" << header.sourcename << "\""<< endl
	   << "\tRA          : " << getHMSString(header.ra) <<endl
	   << "\tDEC         : " << getDMSString(header.dec) <<endl
	   << "\tSTART TIME  : " << setprecision(10) 
	   << header.starttime << endl
	   << "\tPICTURE     : " << ri.picthresh << endl
	   << "\tBOUNDARY    : " << ri.bndthresh << endl;
    stream << "\tDEROTATION  : ";
    if(ri.derotation) 
	stream << "enabled" << endl;
    else
	stream << "disabled" << endl;
    
    
    int ntel = array.getNumTelescopes();
    stream << "\tTELESCOPES  : " << ntel <<endl;

    stream << "\tNADC        : ";
    for (int i=0; i<header.nadc.size(); i++) 
	stream <<  header.nadc.at(i) << " ";
    stream << endl;

    stream << "\tCAMERA ID   : ";
    for (int i=0; i<ntel; i++) 
	stream << array.getCamera(i)->getCameraID() <<" ";
    stream << endl;


    stream << "\tPIXELS      : ";
    for (int i=0; i<ntel; i++) 
	stream << array.getCamera(i)->getNumPixels() <<" ";
    stream << endl;

    stream << "\tSIGMA_PIX   : ";
    for (int i=0; i<ntel; i++) 
	stream << array.getCamera(i)->getSigmaPixel() <<" ";
    stream << endl;

    stream << "\tPE TO DC    : ";
    for (int i=0; i<ntel; i++) 
	stream << array.getCamera(i)->getPEToDC() <<" ";
    stream << endl;
    
    stream << "\tMASKED TUBES: ";
    vector<int> masked;
    for (int i=0; i<1024; i++) {
	if (ri.tubemask[i] == true)
	    masked.push_back(i);
    }
    if (masked.size() >0) {
	for (int i=0; i<(int)(masked.size()); i++)
	    stream << masked[i] << " ";
    }
    else {
	stream << "none";
    }
    stream << endl;

    stream << "\tCAM OFFSET  : ";
    if(ri.camera_offset_analysis==true) {
	stream << ri.camera_offset.x << " , "<<ri.camera_offset.y<<endl;
    }
    else {
	stream << "disabled " << endl;
    }

    stream << "\tIMAGE FILTER: ";
    if(ri.image_filter==true) {
	stream << "amt="<<ri.filter_amount 
	       << " sub="<<ri.filter_subdivide
	       << " pct="<<ri.filter_percent
	       << " thr="<<ri.filter_threshold
	       << endl; 
    }
    else {
	stream << " disabled " << endl;
    }

}


/**
 * Process one run (called by process()). This is the main
 * parameterization routine.
 *
 * \todo: This should be broken down into smaller pieces - it's huge
 * and clunky.
 *
 * \todo: maybe make the display run in its own thread so that it can
 * update at a slower rate as the parameterization is proceeding.
 *
 * \bug: need to fix the RA/Dec override to deal with OFF-ON runs.
 * Should specify time difference (30 min) and ON-OFF order
 *
 * \bug: multiple cameras are not correctly handled - right now, all
 * cameras are assumed to be the same (when TelescopeArray is
 * initialized, it takes a single value for the number of
 * pixels). There needs to be a more intelligent mechanism. Probably,
 * there should be a number of telescopes field in the header for when
 * array data is used (set to 1 by default)
 */
void
Parameterizer::
processOne(  RunInfo &ri, string &id ) {

    register int i;


    PedestalFinder ped;
    GainFinder gf;
    RawDataReader *data =NULL;
    ParamDataWriter *writer=NULL;
    ImageCleaner *threshcleaner = NULL;
    HillasParameterization param;
    MuonParameterization mparam;
    RawHeaderRecord header;
    RawEventRecord  event;
    Array_t adc;
    ofstream muonfile;
    ofstream goodmuons;

    vector<TelescopeData> tel(1);

    vector<int> cleanpixels;
    cleanpixels.reserve(100);

    string filename;
    string fname;
    double start_el;

    int datasize=0;
    string tempstr;
    string answer;
    double gpstime_start=0;
    double gpstime_last=0;
    double deviate;
    int dtype;
    double zen;

    SimShowerRecord cursimshower;

    CutFactory *cutfactory = CutFactory::instance();
    SuperCutter *cutter;
    cutter = cutfactory->newCutter( ri );    
    cutter->setCuts( ri.cuts );

    Logger *logger = Logger::instance();

    // for interactive mode:
    bool interactive_cuts=false;
    bool interactive_triplet_centers=false;
    bool interactive_hold=false;
    bool interactive_muon=false;
    bool interactive_rotate=false;
    bool interactive_correct=true;
    bool interactive_clean = true;
    bool interactive_params = true;
    PlotMaker::CameraPlotType interactive_disptype=PlotMaker::CAMERA_IMAGE;

    int ison;
    if (ri.onid == id) ison = 1;
    else ison=0;

    // DEBUGGING: write out raw data
//     ofstream rawfile( string(id+"-rawdata.txt").c_str() );
    
    
    //=========================================================================
    // Set up output directory:

    string dir = id + "/";
    if (mkdir(dir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR | S_IXGRP | S_IRGRP )) {
	if (errno != EEXIST) {
	    perror("mkdir");
	    throw CriticalAnalysisException("Couldn't create output directory!");
	}
    }

    //=========================================================================
    // Open the output files:

    try {

	switch ( ri.outtype ) {
	  case Config::NTUPLE:
	      writer = new ParamDataWriterNtuple( dir+id+"-param.ntuple",
						  _is_overwrite_enabled);
	      break;
	  case Config::TEXT:
	      writer = new ParamDataWriterText( dir+id+"-param.txt",
						_is_overwrite_enabled );
	      break;
	}

	if (ri.muoncalibrate) {
	    
	    string muonfilename = dir+id+"-muons.txt";
	    muonfile.open( muonfilename.c_str() );

	    muonfilename = dir+id+"-goodmuons.txt";
	    goodmuons.open( muonfilename.c_str() );

	}


    }
    catch (MildAnalysisException &e) {
	cerr << "Skipping run since "<<e.what() << endl;
	return;
    }

    //=========================================================================
    // Open the input file:


    RawDataReaderFactory *datafactory = RawDataReaderFactory::instance();
    
    try {

	data = datafactory->getReader( ri, id );

    }
    catch (AnalysisException &e) {
	delete data;
	throw (e);
    } 

    dtype = data->getType();
    
    //=========================================================================
    // Get some information: 

    datasize = data->size();
    if (datasize == 0) {
	delete data;
	throw AnalysisException("There is no data in "+id
				+"! Check that the file is not corrupt, "
				+"or remove it from the conf file");
    }


    
    // fetch the header from the data

    data->getHeaderRecord( header );

    cout << "DEBUG: header.nadc.size() = "<<header.nadc.size()<<endl;

    // Fill in unknown values from data file  and check for misinformation:

    removeWhiteSpace(header.sourcename);
    removeWhiteSpace(ri.sourcename);
    transform( header.sourcename.begin(),header.sourcename.end(), 
	       header.sourcename.begin(), (int(*)(int))tolower);

    if (ri.sourcename != "") {
	if (ri.sourcename != header.sourcename) {
	    logger->printf("Sourcename may be bad for '%s': conf says '%s', "
			   "data says '%s'", id.c_str(),
			   ri.sourcename.c_str(), header.sourcename.c_str() );
	    
	}
	header.sourcename = ri.sourcename;
    }

    if (ri.ra != 0.0)  {	// if there's an ra override
	if (ri.ra != header.ra) {
	    logger->printf("RA may be incorrect for '%s': conf says '%f',"
			   " data says '%f'",
			   id.c_str(), ri.ra, header.ra);	    
	    logger->printf("RA Override currently assumes ON before OFF "
			  "and 30minute offset");
	}
	if ( id == ri.offid) 
	    header.ra = ri.ra + (30*((M_PI/12.0)/60.0));  // offset 30 minutes
	else
	    header.ra   = ri.ra;
    }
    if (ri.dec != 0.0) {	// if there is a dec override
	if (ri.dec != header.dec) {
	    logger->printf( "DEC may be incorrect for %s: "
			    "conf says '%f', data says '%f'", 
			    id.c_str(), ri.dec, header.dec);
	}
	header.dec  = ri.dec;
    }
    if (ri.ntubes != 0)      header.nadc[0] = ri.ntubes;

    // warn if ra and dec are 0 (usually means there was a problem with
    // the tracking computer when the run was being taken:

    if (header.ra==0 && header.dec==0) {
	logger->printf("RA and Dec are 0 for %s! "
		       "Maybe you need to specify them?", id.c_str());
    } 

    // set up progress bar

    ProgressBar progress( datasize,id );

    //=========================================================================
    // Now that the header has been loaded, we can get the array
    // information.  The array object will automatically load the
    // correct cameras.  Right now, there's just one "nadc" value
    // stored in the header, so I assume all cameras are the same.
    // Eventually, there needs to be a mechanism to detect different
    // cameras in the data. In principle, there's no reason that there
    // can't be mixed cameras.

    TelescopeArray array(ri.arrayconfig, 
			 header.nadc, atoi(ri.utdate.c_str()) );

    int numtel = array.getNumTelescopes();
    int curtel = 0; // current telescope

    tel.resize(array.getNumTelescopes());

    // KPK: NOTE: header.nadc[k] is the number of adc values that the data
    // file is reporting while array.getCamera(k)->getNumPixels() is
    // the number of pixels in the camera according to the camera
    // definition file.  These numbers really should be the same, but
    // in some cases there nadc is larger than numPixels when extra
    // adc channels are installed on the telescope for monitoring
    // other data. Hopefully, I've used the correct one in each case
    // to resize arrays, but if not it may cause things to crash in
    // certain circumstances.


    if (numtel != header.num_telescopes) 
	logger->printf( "The data file (%s) contains %d telescopes, while "
			"the array config (%s) specifies %d. "
			"Are you sure you are using the correct array "
			"configuration?", 
			id.c_str(), header.num_telescopes, 
			ri.arrayconfig.c_str(), numtel);


    //=========================================================================
    // Fetch the pedestal/gain information for each telescope in the array

    for (int k=0; k<numtel; k++) {

	ped.getPeds( ri, id, tel[k].peds, k );
	gf.getGains( ri, tel[k].gains, k );
    
	if (ri.type == Config::ONOFF) {
	    if (ison) 
		ped.getOffSourcePeds( ri, tel[k].padpeds,k );
	    else 
		ped.getOnSourcePeds( ri, tel[k].padpeds,k );
	}
	else if (ri.type == Config::TRACK && ri.pad_track_with_run==true) {
	    ped.getOffSourcePeds( ri, tel[k].padpeds,k );
	}
	
    

	// Pad to fixed value if a tracking run and no padding run is specfied
	if (ri.padding == true && ri.type == Config::TRACK 
	    && ri.pad_track_with_run == false) {
	    
	    tel[k].padpeds.resize( tel[k].peds.size() );
	    for (i=0; i<(int)(tel[k].padpeds.size()); i++) {
		tel[k].padpeds[i].dispersion = ri.padlevel;
		tel[k].padpeds[i].pedestal = 0;
	    }
	}

    }




    //=========================================================================
    //  HillasImageAnalyzer, ImageCleaner, etc.

    Camera *curcam = array.getCamera(0);

    // if the user requests analysis centered around a position other
    // than 0,0, we need to shift the point of origin of the camera before
    // parameterization:

    if (ri.camera_offset_analysis == true) {
	
	for (int k=0; k<numtel; k++) {
	    array.getCamera(k)->shiftCameraCoordinates( ri.camera_offset.x, 
							ri.camera_offset.y );
	}
    }

    
    // stuff for high-resolution camera generated by the ImageFilter routines
    // right now, image filtering is only implemented for a single telescope.
    Camera* hires_cam=NULL;
    Array_t hires_x,hires_y, hires_adc;
    Array_t camrotx(curcam->xCoords());
    Array_t camroty(curcam->yCoords());
    vector<Pedestal> hires_peds;
    double s_max;
    vector<int> hires_cleanpixels;
    hires_cleanpixels.reserve(500);
    
    if (ri.image_filter == true) {
	// for image filter method, set up the analyzer to use the
	// high-res version of camera generated by the filter program.

	// de-rotate the camera (must be aligned with axis):
	
	double angl = atan2( camroty[1]-camroty[0],camrotx[1]-camrotx[0] );
	Coordinate_t pt;

	if (angl>M_PI/2.0) angl = -(M_PI - angl);

	cout << "ImageFilter: De-rotating the camera by "<<angl*180.0/M_PI
	     <<" degrees"<<endl;


	for (i=0; i<camrotx.size(); i++) {
	    pt.x = camrotx[i];
	    pt.y = camroty[i];
	    derotatePoint( angl, pt );
	    camrotx[i] = pt.x;
	    camroty[i] = pt.y;
	}
	
	filter_init( curcam->getLastPixel()+1, // actual number of pixels
		     camrotx, camroty, 
		     hires_x, hires_y, 
		     ri.filter_subdivide,  // msub
		     ri.filter_amount, // afilter 
		     ri.filter_iter, // nc
		     ri.filter_smoothness);  //spas );

	hires_cam = new Camera( hires_x.size(), hires_x, hires_y );
	tel[0].analyzer = new HillasImageAnalyzer( hires_cam );

	hires_peds.resize(hires_x.size());
	for ( i=0; i<hires_peds.size();i++){
	    hires_peds[i].pedestal=0.0;
	    hires_peds[i].dispersion=0.0;
	}
    }
    else {
	//  standard method: set up analyzer to use the camera
	for (int k=0; k<numtel; k++) {
	    tel[k].analyzer = new HillasImageAnalyzer( array.getCamera(k) );
	}
    }  

    for (int k=0; k<numtel; k++) {
	tel[k].muonanalyzer = new MuonImageAnalyzer( array.getCamera(k) );
    }

    // set up the plotters for interactive mode:

    PlotMaker *pm=NULL;
    PlotMaker *pspm=NULL;
    PlotMaker *hires_pm=NULL;
    if (_is_interactive || _is_display_enabled ) {
	pm = new PlotMaker( "X", "x.out",
 			    curcam->getMinX() - 0.4, 
 			    curcam->getMinY() - 0.4, 
 			    curcam->getMaxX() + 0.4, 
 			    curcam->getMaxY() + 0.4 
 			    );
	if (ri.image_filter==true) {
	    hires_pm = new PlotMaker( "X", "x2.out",
				      hires_cam->getMinX() - 0.4, 
				      hires_cam->getMinY() - 0.4, 
				      hires_cam->getMaxX() + 0.4, 
				      hires_cam->getMaxY() + 0.4 
				      );
	}
    }

    
    //=========================================================================
    // Set up the histograms... These get deleted automatically when
    // the tel[] array is destructed

    for (int k=0; k<numtel; k++) {

	curcam = array.getCamera(k);
	int numtubes = curcam->getLastPixel() - curcam->getFirstPixel();
	cout << "DEBUG: "<<k<<": numtubes="<<numtubes<<endl;

	tel[k].lensize = new Histogram(100, 0.0001, 0.002, "length/size");
	tel[k].pictubes = new Histogram( array.getCamera(k)->getNumPixels(), 
					 0, array.getCamera(k)->getNumPixels(), 
					 "tubes in picture");
	tel[k].tubehits= new Histogram( numtubes, 
					curcam->getFirstPixel(),
					curcam->getLastPixel(), "tube hits");
	tel[k].phihist= new Histogram( 180, 0.0, 2.0*M_PI,  "phi angle");
	tel[k].alphahist= new Histogram( 18, 0.0, 90.0, "raw alpha" );
	tel[k].lengthhist= new Histogram( 50, 0.0, M_PI_2, "length" );
	tel[k].sizehist= new Histogram( 50, 0.0, 3000, "size" );
	tel[k].disthist= new Histogram( 50, 0.0, M_PI_2, "distance" );
	tel[k].widthhist= new Histogram( 50, 0.0, M_PI_2, "width" );
	tel[k].rawrate= new Histogram( 30,0,30, "Raw Rate" );
	tel[k].deltat= new Histogram(100,0,0.00001,"DeltaT");
	tel[k].max2hist= new Histogram( 1000,0,1000,"Max2" );
	tel[k].centroids= new Image2D(39,39);
	
	// for muon calibration
	tel[k].musoalhist= new Histogram( 100,0,4000,"signal/arclength" );
	tel[k].musizehist= new Histogram( 50,0.0,2000.0, "Muon size");
	tel[k].muonnesshist= new Histogram( 50,0.0,1.0, "Muonness");
	tel[k].mugainhist= new Histogram( 50,0,500, "Muon gain factor");
	tel[k].mulensizehist= new Histogram( 100,0.0001,0.002, "Muon Length/Size");
	
	tel[k].centroids->setCoordinateBox(-2,-2,  2,2);

    }

    vector<ofstream>  anglefile(numtel);
    string anglefilename;
    for (int k=0; k<numtel; k++) {
	char telid[100];
	sprintf( telid, "-%03d", k );
	string anglefilename = dir+id+telid+"-angles.dat";
	anglefile[k].open(anglefilename.c_str());
    }

    //=========================================================================
    // Print run information to file and screen (if requested)

    if (_is_verbose_enabled){
	printSummary( cout, ri, id, header, writer, data, array );
    }
    ofstream runinfo( string(dir+id+"-run_info.txt").c_str() );
    printSummary( runinfo, ri, id, header, writer, data, array );
    
    //=========================================================================
    // if Padding is enabled, precalculate the 
    // max pedestal dispersions for the run pair

    for (int k=0; k<numtel; k++) {
	
	tel[k].combinedpeds.resize( tel[k].peds.size() );
	tel[k].sigmaped.resize( tel[k].peds.size() );

	if (ri.padding == true ) {
	    
	    for (i=0; i<tel[k].combinedpeds.size(); i++) {
		
		if (tel[k].peds[i].dispersion > tel[k].padpeds[i].dispersion) 
		    tel[k].combinedpeds[i] = tel[k].peds[i];
		else
		    tel[k].combinedpeds[i] = tel[k].padpeds[i];
		
		// Note that if either padding or current run has a tube
		// off or a star, the attributes of combined peds take on
		// that of the current run.
		
		if (tel[k].padpeds[i].type != Pedestal::GOOD ){
		    tel[k].combinedpeds[i].type = tel[k].padpeds[i].type;
		}
		
		if (tel[k].peds[i].type != Pedestal::GOOD ){
		    tel[k].combinedpeds[i].type = tel[k].peds[i].type;
		}
		
		// Precalculate the pedestal sigmas
		
		tel[k].sigmaped[i] = pow(tel[k].combinedpeds[i].dispersion,2) 
		    - pow(tel[k].peds[i].dispersion,2);
		
		if (tel[k].sigmaped[i] > 0) tel[k].sigmaped[i] = sqrt(tel[k].sigmaped[i]);
		else tel[k].sigmaped[i] = -1.0e-5;
		
	    }
	    
	}
	else {
	    tel[k].combinedpeds = tel[k].peds;
	}
    }

    // calculate average pedestal dispersion. This is used when image
    // filtering is enabled to make the pedestal cut.
    // TODO: fix this to work with multiple telescopes
    double avgpeddisp=0;
    int navg=0;
    for (i=0; i<tel[0].sigmaped.size(); i++) {
	if (tel[0].combinedpeds[i].type==Pedestal::GOOD){
	    avgpeddisp += tel[0].combinedpeds[i].dispersion;
	    navg++;
	}	
    }
    avgpeddisp /= (double)navg;

    // Apply the user-specified tube mask from the config file:
    // TODO: implement for multiple telescopes

    for ( i=0; i<tel[0].combinedpeds.size(); i++) {
	if (ri.tubemask[i] == 1) {
	    tel[0].combinedpeds[i].type = Pedestal::TUBEOFF;
	    tel[0].peds[i].type = Pedestal::TUBEOFF;
	    tel[0].padpeds[i].type = Pedestal::TUBEOFF;
	}
    }

    // Initialize the imagecleaner using the combined pedestal values

    for (int k=0; k<numtel; k++) {
	curcam = array.getCamera(k);
	tel[k].cleaner = new DefaultImageCleaner( *curcam,tel[k].combinedpeds,
						  ri.picthresh,ri.bndthresh );
    }

    
    if (ri.image_filter == true ) {
	threshcleaner = new ThresholdImageCleaner( *hires_cam, 0.0 );
    }
    

    //=========================================================================
    // Loop over events in the data and parameterize...

    vector<int> n(numtel);
    vector<int> n_trig(numtel);
    vector<int> n_other(numtel);
    int count =0;

    for (int k=0; k<numtel; k++) {
	n[k] = n_trig[k] = n_other[k] = 0;
    }
    
    

    while ( data->isDone() == false ) {

	// ---------------------------------------------------------
	// Get events until no more are found 
	// ---------------------------------------------------------

	try {

	    if (interactive_hold == false) 
		data->getNextEventRecord( event );  
	    else 
		cout << "INTERACTIVE MODE: holding current event"<<endl;

	}
	catch (MildAnalysisException e) {
	    // Mild exceptions should be noted, but processing should
	    // still continue.  Anything stronger will cause this to
	    // return and will be caught at a higher level.
	    cerr << endl << "WARNING: " << e.what() << endl;
	    break;    
	}
	
	if (n[0] == 0) {
	    // it's the first event, so store the start time
	    gpstime_start = event.gpstime;
	}

	// store a copy of the adc values
	if (adc.size() != event.adc.size()) adc.resize(event.adc.size());
	adc = event.adc;

	// --------------------------------------------------------
	// Parameterize each event if it's the right type:
	// ---------------------------------------------------------
	if (event.type == RawDataReader::EVENT) {

	    // Figure out which telescope to use:
	    if (event.telescope_id >= numtel) {
		cout << "WARNING: "
		     << "Skipping event with telescope id "<<event.telescope_id
		     << " (there are only "<<numtel
		     <<" telescopes in the array!)"
		     << endl;
		continue;
	    }

	    curcam = array.getCamera(event.telescope_id);
	    curtel = event.telescope_id;

	    // ----------------------------------------- 
	    // Subtract pedestals and pad if requested:
	    // -----------------------------------------

	    for ( i=0; i<(int)event.adc.size(); i++){

		// If the tube is on (according to the pedestal),

		if (tel[curtel].peds[i].type == Pedestal::GOOD) { 

		    if (ri.padding == true) {

			// "Pad" to bring the noise up to the maximum
			// of the on and off peddisps.  Currently, add
			// in noise based on a gaussian distribution,
			// but a poisson distribution may be better
			// for low statistics. That would require
			// switching gsl_ran_gaussian() to
			// gsl_ran_poisson() and changing the formula
			// a bit.


			if (tel[curtel].sigmaped[i] > 0.0) {

			    deviate = _gaussdev->getDeviate();
			    adc[i] += deviate*tel[curtel].sigmaped[i];

			}

		    }

		    adc[i] -= tel[curtel].peds[i].pedestal;
		    
		}

		// Otherwise mark the tube as off

		else {
		    adc[i] = 0; 
		}

	    }

	    // ----------------------------------------- 
	    // Apply gains:
	    // ----------------------------------------- 
	    
	    adc *= tel[curtel].gains;

	    // ----------------------------------------- 
	    // Clean the image
	    // and generate a linked list of pixels that
	    // are in the picture and boundary.  
	    // ----------------------------------------- 

	    if (ri.image_filter == false) {
	       
		// standard cleaning method:
		 
		cleanpixels = tel[curtel].cleaner->getCleanPixels( adc );

		// add the clean pixels to a histogram to look for
		// misfiring tubes 
		
		// TODO: speed this up by not using
		// increment() since it's integers anyway. Right now, I
		// add 0.5 to make sure the correct bin is incremented
		// (since the value lies on the bin-boundary, which seems
		// to cause random behavior due to rounding)
		
		for (i=0; i<cleanpixels.size(); i++) {
		    tel[curtel].tubehits->increment(cleanpixels[i]+0.5);
		}

	    }
	    else {
		// Image filter method: smooth image with hexagonal
		// FFT filtering.  The resulting image has more pixels
		// than the original, since the resolution of the
		// camera is modified.

		// first, get the regular cleaning for comparison:
		cleanpixels = tel[curtel].cleaner->getCleanPixels( adc );

		// now filter the image:
		s_max = filter_image( adc,hires_adc );


		// clean the image with a threshold based on a
		// percentage of the maximum value after subtracting
		// the picture threshold:
		((ThresholdImageCleaner*)(threshcleaner))
		    ->setThreshold( ri.filter_threshold*avgpeddisp 
				    + (s_max-ri.filter_threshold*avgpeddisp)
				    *ri.filter_percent );
		hires_cleanpixels = threshcleaner->getCleanPixels( hires_adc );

	    }



	    // ----------------------------------------- 
	    // Get the image parameterization
	    // ----------------------------------------- 

	    updateAngles( event.gpstime, header );
	    zen =  M_PI_2 - _el; 	   // zenith angle

	    if ( ri.image_filter == true )
		tel[curtel].analyzer->parameterize( hires_adc, hires_cleanpixels );
	    else
		tel[curtel].analyzer->parameterize( adc, cleanpixels );




	    param = tel[curtel].analyzer->getHillasParameters();
	    param.osctime = event.osctime;
	    param.livetime = event.livetime;
	    param.gpstime = event.gpstime;
	    param.event_number = _event_number;
	    param.zenith =  zen; 
	    param.telescope_id = event.telescope_id;

	    // if there's a zenith override
	    if (ri.zenithoverride > 1e-9) {
		param.zenith = ri.zenithoverride;
	    }

	    // ----------------------------------------- 
	    // Get the muon parameters if requested
	    // ----------------------------------------- 

	    if ( ri.muoncalibrate ) {

		// muons usually have a size less than 2000, so lets
		// only calculate muon parameters when this occurs.
		// Should define a muon size cut in the Config file!

 		if ( param.size < 3000 ) {
		    
		    try {
			tel[curtel].muonanalyzer->parameterize( adc,cleanpixels );
			mparam = tel[curtel].muonanalyzer->getMuonParameters();
		    }
		    catch (MildAnalysisException e) {
			
		    }
		    
		    tel[curtel].muonnesshist->increment( mparam.muonness );

		    if ( mparam.muonness > ri.muon_threshold ) {
			tel[curtel].musoalhist->increment( mparam.soal );
			tel[curtel].musizehist->increment( param.size );
			tel[curtel].mugainhist->increment( mparam.mugain );
			tel[curtel].mulensizehist->increment( param.length_over_size );
		    }
		}
		else {
		    mparam.clear();
		}
		
		// write out the parameters
		
		muonfile << param.event_number << " " << mparam << endl;
				
 	    }
	    
	    
	    // ----------------------------------------- 
	    // Display if interactive:
	    // ----------------------------------------- 

	    if ( _is_display_enabled ) {
		
		if (_is_interactive == true) {


		    cout << "PARAMS   : "<<param << endl;

		    if (interactive_correct) {
			cutter->applyCorrections(param);
			cout << "CORRECTED: " << param << endl;
		    }

  
		    cout << "Clean Pixels: ";
		    for (i=0; i<cleanpixels.size(); i++)
			cout << cleanpixels[i] << " ";
		    cout << endl;
		    if (ri.muoncalibrate)
			cout << "MUON PARAMETERS: "<< mparam << endl;


		    tempstr = "";
		    while (tempstr != "n") {

			if (interactive_triplet_centers) 
			    tel[curtel].muonanalyzer->enablePlotPoints(true);
			else 
			    tel[curtel].muonanalyzer->enablePlotPoints(false);
			
			if (interactive_muon) {
			    // seek to next muon
			    if (mparam.muonness > ri.muon_threshold) { 
				interactive_muon=false;
			    }
			    else {
				break;
			    }
			    
			}

			
			if (interactive_cuts) {
			    cout << "Seeking: "<<param.event_number 
				 << "     \r" << flush;

			    if (cutter->pass(param,(SIZE|MAX|FRAC|LENSIZE
						       |DISTANCE|LENGTH
						       |WIDTH))){
				cout << "Event " << param.event_number
				     << " Passed cuts    " << endl;
			    }
			    else break;
			} // End if cuts are enabled in interactive mode.


			// Plot the camera

			pm->setAxisBox( *curcam );
			if (interactive_rotate){
			    pm->pushState();
			    pm->rotate( _theta );
			}

			if (interactive_clean == false) {
			    cleanpixels.clear();
			}


			if (ri.image_filter==true) {
			    hires_pm->setAxisBox( *hires_cam );
			    hires_pm->plot( *hires_cam, hires_adc, 
				      hires_peds,hires_cleanpixels,
				      interactive_disptype);
			    hires_pm->plot( param );
			    hires_pm->plotAxes();
			    hires_pm->plotTitle("Filtered Image");
			    hires_pm->newPage();
			}

			pm->plot( *curcam, adc, tel[curtel].combinedpeds, 
				  cleanpixels, interactive_disptype );
			
			
			if (interactive_params) 
			    pm->plot( param );

			if (ri.muoncalibrate) pm->plot( mparam );
			pm->plotAxes();
			switch (interactive_disptype) {
			case PlotMaker::CAMERA_IMAGE:
			    pm->plotTitle("Camera image - "+header.sourcename);
			    break;
			case PlotMaker::CAMERA_PEDS:
			    pm->plotTitle("Peds - "+header.sourcename);
			    break;
			case PlotMaker::CAMERA_PEDDISPS:
			    pm->plotTitle("sigma_ped - "
					  +header.sourcename);
			    break;
			case PlotMaker::CAMERA_TUBENUMS:
			    pm->plotTitle("Tube numbers - "
					  +header.sourcename);
			    break;
			default:
			    pm->plotTitle("Camera image - "+header.sourcename);
			    break;
			}
		    
			if (interactive_rotate)pm->popState();
			pm->newPage();
		       			
			cout << "-----------------------------------------"
			     << endl
			     << "[n] Next event         [s] Save as ps    "
			     << "      [p] set Plot type" << endl
			     << "[m] next Muon          [r] Redraw        "
			     << "      [i] set dIsplay scale" <<endl
			     << "[,] mark good muon     [c] toggle Cuts (";
			if (interactive_cuts) cout << "*"; else cout <<" ";
			cout <<")     [C] toggle corrections (";
			if (interactive_correct) cout << "*"; else cout <<" ";
			cout << ")" << endl;
			cout << "[h] Hold event (";
			if (interactive_hold) cout << "*"; else cout <<" ";
			cout <<")     [R] toggle Rotation (";
			if (interactive_rotate) cout << "*"; else cout <<" ";
			cout <<") [t] toggle triplet cntrs (";
			if (interactive_triplet_centers) 
			    cout << "*"; else cout <<" ";
			cout <<")"<<endl;
			cout << "CHOICE: ";
			cin >> tempstr;
			cout << tempstr[0] << endl;

			switch (tempstr[0]) {
			case 'c':
			    interactive_cuts = !interactive_cuts;
			    break;
			case 'C':
			    interactive_correct = !interactive_correct;
			    break;
			case 'h':
			    interactive_hold = !interactive_hold;
			    break;
			case 'l':
			    interactive_clean = !interactive_clean;
			    break;
			case 'P':
			    interactive_params = !interactive_params;
			    break;
			case 't':
			    interactive_triplet_centers = 
				!interactive_triplet_centers;
			    break;
			case 's':
			    cout << "Enter filename: ";
			    cin >> fname;
			    pspm = new PlotMaker( "PS", fname.c_str() );
			    pspm->setAxisBox( *curcam );
			    pspm->plot( *curcam, adc, tel[curtel].combinedpeds,
					cleanpixels,interactive_disptype);
			    if (interactive_params)
				pspm->plot( param );
			    if (ri.muoncalibrate) pspm->plot( mparam );
			    pspm->plotAxes();
			    switch (interactive_disptype) {
			    case PlotMaker::CAMERA_IMAGE:
				pm->plotTitle("Camera image - "+
					      header.sourcename);
				break;
			    case PlotMaker::CAMERA_PEDS:
				pm->plotTitle("Peds - "+header.sourcename);
				break;
			    case PlotMaker::CAMERA_PEDDISPS:
				pm->plotTitle("sigma_ped - "
					      +header.sourcename);
				break;
			    case PlotMaker::CAMERA_TUBENUMS:
				pm->plotTitle("Tube numbers - "
					      +header.sourcename);
				break;
			    default:
				pm->plotTitle("Camera image - "
					      +header.sourcename);
				break;
			    }
			    delete pspm;

			    break;
			case 'S':
			    if (ri.image_filter==true) {
				cout << "Enter filename: ";
				cin >> fname;
				pspm = new PlotMaker( "PS", fname.c_str() );
				pspm->setAxisBox( *hires_cam );
				pspm->plot( *hires_cam, hires_adc, 
					    hires_peds, hires_cleanpixels,
					    PlotMaker::CAMERA_IMAGE);
				if (interactive_params)
				    pspm->plot( param );
				if (ri.muoncalibrate) pspm->plot( mparam );
				pspm->plotAxes();
				pm->plotTitle("Filtered Camera image - "+
					      header.sourcename);
				delete pspm;
			    }
			    else cout << "Image Filter is disabled"<<endl;
			    break;

			case 'd':
			    for (int k=0; k<100; k++) {
				adc[k] -= 0.3*adc[k];
			    }
			    break;
			case 'r':
			    break;
			case 'i':
			    cout << "ENTER NEW MAX ADC VALUE (-1 for auto):";
			    cin >> answer;
			    pm->setPixelScale(atoi(answer.c_str()));
			    break;
			case 'm':
			    if (ri.muoncalibrate==true){
				interactive_muon = true;
				tempstr = "n";
			    }
			    else {
				cout << "MuonCalibration is currently "
				     << "disabled! " << endl;
			    }
			    break;
			case ',':
			    goodmuons << param.event_number << " " 
				      <<mparam << endl;
			    cout << "Saved event to goodmuons.txt"<<endl;
			    break;
			case 'p':
			    cout << "CHOOSE PLOT TYPE: 1) image 2) peds "
				 << "3) ped disps 4) tube nums: ";
			    cin >> answer;
			    switch ( atoi(answer.c_str())){
			    case 1:
				interactive_disptype=PlotMaker::CAMERA_IMAGE;
				break;
			    case 2:
				interactive_disptype=PlotMaker::CAMERA_PEDS;
				break;
			    case 3:
				interactive_disptype=
				    PlotMaker::CAMERA_PEDDISPS;
				break;
			    case 4:
				interactive_disptype=
				    PlotMaker::CAMERA_TUBENUMS;
				break;
			    default:
				cout << "INVALID OPTION: "<<answer<<"!"<<endl;
			    }
			    break;
			case 'R':
			    interactive_rotate = !interactive_rotate;
			    break;
			case 'n':
			    break;  // don't print a warning
			default:
			    cout << "Unknown command '"<< tempstr << "'"<< endl;
			}
			    
		    } // End while

		} // End if interactive
		else {
		    // non-interactive, just display
		    pm->setAxisBox( *curcam );
		    pm->plot( *curcam, adc, tel[curtel].combinedpeds, cleanpixels,
			      interactive_disptype);
		    pm->plot( param );
		    if (ri.muoncalibrate) pm->plot( mparam );
		    pm->plotAxes();
		    pm->newPage();
		} // End else not interactive
			

	    } // End if diplay enabled

	    // -----------------------------------------
	    // derotate the point of origin and centroid
	    // -----------------------------------------
	    
	    // write out the angles every so often for diagnostics:
	    if (n[curtel] % 256==0) {
		anglefile[curtel].precision(15);
		anglefile[curtel] << (event.gpstime-gpstime_start)*24*60<< "\t" 
				  << _theta *180.0/M_PI  << "\t"
				  << _el * 180.0/M_PI << "\t"
				  << _az * 180.0/M_PI 
				  << endl;
	    }

	    if (ri.derotation == true) {

		derotatePoint( _theta, param.point_of_origin_a ); 
		derotatePoint( _theta, param.point_of_origin_b ); 
		derotatePoint( _theta, param.centroid );

	    }

	    param.on=ison;

	    // ----------------------------------------- 
	    // Write parameters to disk and update
	    // diagnostic histograms
	    // ----------------------------------------- 

            if (data->getType()==RawDataReader::SIM) {
                cursimshower=((RawDataReaderSim*)(data))->getSimShowerRecord();
		param.sim = cursimshower;
	    }
	    else {
		memset( &(param.sim), 0, sizeof(SimShowerRecord));
	    }
	    
	    writer->writeParameterization( param );
	    	    
	    if (param.invalid==0) {
		tel[curtel].lensize->increment( param.length_over_size );
		tel[curtel].pictubes->increment( param.pixels_in_picture );
		tel[curtel].phihist->increment( param.phi );
		tel[curtel].alphahist->increment( param.alpha*180.0/M_PI );
		tel[curtel].lengthhist->increment( param.length );
		tel[curtel].widthhist->increment( param.width );
		tel[curtel].rawrate->increment((param.gpstime-gpstime_start)*24*60 );
		tel[curtel].disthist->increment( param.distance );
		tel[curtel].sizehist->increment( param.size );
		tel[curtel].deltat->increment( param.gpstime-gpstime_last );
		tel[curtel].max2hist->increment( param.max[2] );
		tel[curtel].centroids->addHist( param.centroid.x,param.centroid.y,1);
	    }

	    gpstime_last = param.gpstime;
	    n_trig[curtel]++;

	    if (n_trig[curtel] == 1) {
		// this is the first trigger event, so store the starting
		// zenith angle
		// TODO: make multi-telescope!
		start_el = _el;
	    }

	}
	else {
	    // not a trigger event
	    n_other[curtel]++;
	}
       

	// --------------------------------------------------------
	// Output progress info if requested
	// --------------------------------------------------------
	if (_is_verbose_enabled && (_is_interactive == false) ) {
	    if (count%256 == 0)
		progress.print( count );
	    if (ri.image_filter==true && count%16==0)
		progress.print(count);
	}
	
	if (interactive_hold == false) {
	    n[curtel]++;
	    count++;
	    _event_number++;
	}

    }

    // DEBUGGING:
//     rawfile.close() ;

    progress.printClear();

    // store the end time in the header to pass on to the cutter
    // (this is needed for rate binning)
    
    header.endtime = event.gpstime;
    header.average_elevation = (start_el + _el) /2.0;

    if (_is_verbose_enabled) {
	cout << "Summary: " << endl
	     << "\tTRIG EVENTS : " ;
	for (int k=0; k<numtel;k++) cout << n_trig[k] << " ";
	cout << endl
	     << "\tOTHER EVENTS: " ;
	for (int k=0; k<numtel;k++) cout << n_other[k] << " ";
	cout << endl
	     << "\tTOTAL EVENTS: " ;
	for (int k=0; k<numtel;k++) cout << n[k] << " ";
	cout << endl
	     << "\tLIVETIME    : " << setprecision(4)
	     <<  event.livetime/60.0 
	     << " min"<< endl
	     << "\tELAPSED GPS : " << setprecision(10) 
	     << (header.endtime-gpstime_start)*24*60<<" min"<< endl
	     << "\tDEAD TIME   : " 
	     << (header.endtime-gpstime_start)*24*60- event.livetime/60.0 
	     << " min "
	     <<endl<<endl;

	runinfo << "Summary: " << endl
	     << "\tTRIG EVENTS : " ;
	for (int k=0; k<numtel;k++) runinfo << n_trig[k] << " ";
	runinfo << endl
	     << "\tOTHER EVENTS: " ;
	for (int k=0; k<numtel;k++) runinfo << n_other[k] << " ";
	runinfo << endl
	     << "\tTOTAL EVENTS: " ;
	for (int k=0; k<numtel;k++) runinfo << n[k] << " ";
	runinfo << endl
	     << "\tLIVETIME    : " << setprecision(4)
	     <<  event.livetime/60.0 
	     << " min"<< endl
	     << "\tELAPSED GPS : " << setprecision(10) 
	     << (header.endtime-gpstime_start)*24*60<<" min"<< endl
	     << "\tDEAD TIME   : " 
	     << (header.endtime-gpstime_start)*24*60- event.livetime/60.0 
	     << " min "
	     <<endl<<endl;

	    
    }
    
    runinfo.close();


    //======================================================
    // generate diagnostic plots and brightness map
    //======================================================

    for (int k=0; k<numtel; k++) {
	char telid[100];
	sprintf( telid, "-%03d", k );
	tel[k].lensize->save(dir+id+telid+"-raw-lensize");
	tel[k].pictubes->save(dir+id+telid+"-raw-pictubes");
	tel[k].tubehits->save(dir+id+telid+"-tubehits");
	tel[k].phihist->save(dir+id+telid+"-raw-phi");
	tel[k].alphahist->save(dir+id+telid+"-raw-alpha");
	tel[k].lengthhist->save(dir+id+telid+"-raw-length");
	tel[k].widthhist->save(dir+id+telid+"-raw-width");
	tel[k].rawrate->save( dir+id+telid+"-raw-rate" );
	tel[k].sizehist->save( dir+id+telid+"-raw-size" );
	tel[k].disthist->save( dir+id+telid+"-raw-dist" );
	tel[k].deltat->save( dir+id+telid+"-deltat" );
	tel[k].centroids->save( dir+id+telid+"-centroids" );
	tel[k].max2hist->save( dir+id+telid+"-raw-max2" );

	if (ri.muoncalibrate) {
	    tel[k].musoalhist->save(dir+id+telid+"-muon-soal");
	    tel[k].musizehist->save(dir+id+telid+"-muon-size");
	    tel[k].mulensizehist->save(dir+id+telid+"-muon-lensize");
	    tel[k].mugainhist->save(dir+id+telid+"-muon-gain");
	    tel[k].muonnesshist->save( dir+id+telid+"-raw-muonness");
	}

	// TODO: last arg is telescope_id
	generatePlots( id, ri.n2id, ri.cachedir, k );
    }

    // generate sky brightness/tubeoffness maps if it's an on or track run
    // (no need to do for off, since it will be negative of on)

    if (id == ri.onid && dtype != RawDataReader::SIM) {
	for (int k=0; k<numtel; k++) {
	    cout << "Generating sky brightness & tubeoffness maps for " 
		 << "telescope "<< k << " run "<< id 
		 <<"..."<< endl;
	    mapSkyBrightness( id, tel[k].peds, tel[k].padpeds, 
			      tel[k].gains, header, 
			      (event.gpstime+gpstime_start)/2.0,
			      ri.derotation,array.getCamera(k),k );

	}
    }
    
    // output a warning if too many tubes are turned off since the
    // skybrightness map will be problematic in the pointing
    // check later on

    for (int k=0; k<numtel; k++) {
	int ntubesoff = 0;
	curcam = array.getCamera(k);
	for (i=curcam->getFirstPixel(); i<=curcam->getLastPixel(); i++) {
	    if (tel[k].combinedpeds[i].type == Pedestal::TUBEOFF) 
		ntubesoff++;
	}
	
	cout << "Tubes physically turned off: "<< ntubesoff << endl;
	if (ntubesoff >= 10 ) {
	    logger->printf("Pointing checks may be invalid for telescope %d "
			   "run %s since "
			   "%d tubes were physically turned off", k,
			   id.c_str(),ntubesoff);
	} 
    }


    //======================================================
    // clean up
    //======================================================

    writer->writeHeader( header );

    if (_is_verbose_enabled) {
	cout << "========================================"
	     << "========================================"
	     << endl<<endl;
    }

    // Clean out the cached file if requested

    if (_is_cleanup_enabled)
	datafactory->clearCache(ri.cachedir);


    for (int k=0; k<numtel; k++) 
	anglefile[k].close();
    
    if (ri.muoncalibrate){
	muonfile.close();
	goodmuons.close();
    }

    if (pm) delete pm;
    if (hires_pm)delete hires_pm;
    if (hires_cam) delete hires_cam;

    delete writer;
    delete data;

}

/**
 * Generate plots using gnuplot
 */
void
Parameterizer::
generatePlots( string &id, string &n2id, string &cachedir, int telescope_id ) {

    string command;
    int ret;
    char telid[100];

    sprintf( telid, "%03d", telescope_id );

    cout << "Generating diagnostic plots for " << id << "..." << endl;

    command = _wupdir + "/gen-diag-plots.sh "+_wupdir+" "
	+id+" "+n2id+" "+cachedir+" "+telid;
    ret = system(command.c_str());

    if (ret) {
	cout << "** Plotting failed.  Make sure you have GNUPlot installed"
	     << endl
	     << "** and the plotting scripts are installed in either " 
	     << _wupdir << endl 
	     << "** or the directory specified by the WUPARAMDIR "
	     << "environment variable" << endl
	     << endl;
    }


}


/**
 *  Do a pointing check using the pedestal variances.  Generates a 2-D
 *  grid that matches the grid used by the Cutter for the 2-D
 *  significance/excess plots which contains the sky brightness at
 *  each point in RA/DEC space.  Also generates a 2-D "tubeoffness"
 *  map of PMT's that were turned off due to extrememe pedestal
 *  variances.
 *
 *  \param onpeds on-source Pedestal array
 *  \param offpeds off-source Pedestal array
 *  \param grid 2d array of doubles where the brightnesses should be stored.
 *  \param header the RawHeaderRecord for the run (specifying RA & DEC)
 *  \param avg_gpstime the average gps time of the run.
 */
void 
Parameterizer::
mapSkyBrightness( const string &id, const vector<Pedestal> &onpeds, 
		  const vector<Pedestal> &offpeds, const Array_t &gain,
		  const RawHeaderRecord &header, double avg_gpstime,
		  bool derotation, Camera *cam, int telescope_id) {
    
    const double DEGPERPIX=0.1;	// degrees per pixel (assuming square pixels)
    const double SIGR2 = 0.02;	// something like the RMS^2 of the PSF
    
    Image2D brightmap(39,39);
    Image2D tubeoffmap(39,39);

    Array_t xcoord, ycoord;
    Array_t pedvardiff;
    Coordinate_t coord;
    int i,j,k;
    double sky_brightness, tubeoffness;
    double r2, gauss_weight;
    int nxbins = brightmap.getXDim();
    int nybins = brightmap.getYDim();

    // calculate the correct derotation angle (theta)

    updateAngles( avg_gpstime, header );

    // precalculate the pedestal variance difference (yes, variances
    // not dispersions)

    pedvardiff.resize(onpeds.size());
    for (i=0; i<onpeds.size(); i++) {

	if ( (onpeds[i].type == Pedestal::TUBEOFF || 
	      offpeds[i].type == Pedestal::TUBEOFF)
	     ||
	     (onpeds[i].type == Pedestal::STAR && 
	      offpeds[i].type == Pedestal::STAR)) {

	    pedvardiff[i] = 0;

	}
	else {
	    pedvardiff[i] = pow(onpeds[i].dispersion,2) 
		- pow(offpeds[i].dispersion,2);
	}	

    }
    
    // Get the camera coordinates and derotate them into tangential coords!

    xcoord.resize(cam->xCoords().size());
    ycoord.resize(cam->yCoords().size());
    xcoord = cam->xCoords();
    ycoord = cam->yCoords();

    if (derotation == true) {
	// not sure this is right: 
	// should we rotate by -theta instead?
	for (i=0; i<(int)xcoord.size(); i++) {
	    coord.x = xcoord[i];
	    coord.y = ycoord[i];
	    derotatePoint( _theta, coord );
	    xcoord[i] = coord.x;
	    ycoord[i] = coord.y;
	}
    }

    // Set up the coordinate grid:

    brightmap.setCoordinateBox( -nxbins*DEGPERPIX/2.0,-nybins*DEGPERPIX/2.0,
				nxbins*DEGPERPIX/2.0,nybins*DEGPERPIX/2.0);
    tubeoffmap.setCoordinateBox( -nxbins*DEGPERPIX/2.0,-nybins*DEGPERPIX/2.0,
				 nxbins*DEGPERPIX/2.0,nybins*DEGPERPIX/2.0);
    
    // Accumulate the sky brightness and tubeoffness:

    for (i=0; i<nxbins; i++) {
	for (j=0; j<nybins; j++) {
	    
	    sky_brightness = 0;
	    tubeoffness = 0;

	    for (k=0; k<(int)xcoord.size(); k++) {
			
		// tube is ON or contains a star in the ON
		// run, or just good in the OFF, so update sky
		// brightness map
		
		r2 = (pow(xcoord[k] - brightmap.getXCoord(i),2) +
		      pow(ycoord[k] - brightmap.getYCoord(j),2));
		gauss_weight = exp(-r2/SIGR2);
		sky_brightness += pedvardiff[k]*gauss_weight
		    *pow(gain[k],2); // flat field
		
		
		if ((onpeds[k].type == Pedestal::TUBEOFF) ||
		    (offpeds[k].type == Pedestal::TUBEOFF) ||
		    (onpeds[k].type == Pedestal::STAR) || 
		    (offpeds[k].type == Pedestal::STAR)) {

		    // tube is OFF, so update tubeoffness map
		    r2 = (pow(xcoord[k] - tubeoffmap.getXCoord(i),2) +
			  pow(ycoord[k] - tubeoffmap.getYCoord(j),2));
		    gauss_weight = exp(-r2/SIGR2);
		    tubeoffness  += 1.0*gauss_weight;

		}
		
	    }
		
		
	    brightmap.addToPixel(i,j, sky_brightness);
	    tubeoffmap.addToPixel( i,j, tubeoffness );
	    
	}
    }
    
    // save the maps...

    char telid[100];
    sprintf( telid, "-%03d", telescope_id );

    brightmap.save( id+"/"+id+telid+"-sky_brightness" );
    tubeoffmap.save( id+"/"+id+telid+"-tubeoffness" );

}


void
Parameterizer::
setLatLonInRadians( double lat, double lon ) {

    _latitude = lat;
    _longitude = lon;

}


/**
 * Recalculate angles - sets _el, _az, _theta to the correct values
 *
 * \param gpstime time in MJD format
 * \param header  contains the RA and DEC information
 *
 * \todo: should 24 in TIME_TO_RADIANS be length of sidereal day?
 */
void
Parameterizer::
updateAngles( double gpstime, const RawHeaderRecord &header ) {
    
    double hourangle;
    double tmp1,tmp2;
    double cosdec, sindec, sinlat, coslat, sinhr, coshr;
    double sk,ck;

    const double TIME_TO_RADIANS = 2*M_PI/24.0;  
    double tt, tt0;
    int imjd = (int)(gpstime);
    double sidereal;
    double uttime;

    // interval of time measured in julian centuries elapsed since 1
    // Jan 2000 at 12h UT:
    tt = (imjd - 51544.5)/36525.0;

    // fractional time in hours since UT midnight:
    // Note: mjd is in units of DAYS
    uttime = (gpstime-imjd) * 24.0; 

    tt0=6.697374558 + (2.400051336e3 * tt)+(2.5862e-5 *  tt*tt);
    sidereal = fmod( uttime * 1.002737909 + tt0 - 7.392333333, 24.0 );
    if (sidereal < 0.0) sidereal += 24.0;
    sidereal *= TIME_TO_RADIANS;

    // the hour angle in radians

    hourangle = sidereal - header.ra;
    
    // precalculate some stuff for speed

    sinlat = sin(_latitude);
    coslat = cos(_latitude);
    sindec = sin(header.dec);
    cosdec = cos(header.dec);
    sinhr  = sin(hourangle);
    coshr  = cos(hourangle);

    // calculate the elevation
    
    tmp1 = sindec*sinlat + cosdec*coslat*coshr;
    if (tmp1 > 1.0) tmp1 = 1.0;
    else if (tmp1 < -1.0) tmp1 = -1.0;
    _el = asin(tmp1);

//     if (_el < 0) {
// 	cout << "DEBUG: Hey! elevation is negative: " << _el << endl;
// 	cout << "\tgpstime = " << gpstime << endl;
//     }

    // calculate azimuth

    tmp1 = -cosdec*coslat*sinhr;
    tmp2 = sindec - sinlat*sin(_el);
    
    if (tmp1>1.0) tmp1=1.0;
    else if (tmp1<-1.0) tmp1 = -1.0;
    if (tmp2>1.0) tmp2=1;
    else if (tmp2<-1.0) tmp2 = -1.0;
    
    _az = atan2(tmp1,tmp2);

    // get theta (derotation angle)
    
    sk = -coslat*sin(_az);
    ck = cos(_el)*sinlat-sin(_el)*coslat*cos(_az);

    if (sk > 1.0) sk = 1.0;
    else if( sk<-1.0) sk = -1.0;
    if (ck > 1.0) ck = 1.0;
    else if( ck<-1.0) ck = -1.0;
    
    _theta = -atan2(sk,ck);

//     gsl_sf_angle_restrict_pos_e( &_el );
//     gsl_sf_angle_restrict_pos_e( &_az );

}


/**
 * Derotate the given point from camera coordinates to tangential coordinates.
 *
 * Note: x-coord is in local direction of RA, y is in local direction
 * of DEC
 *
 * \todo: check whether sign of x is same as sign of RA
 */
void  
Parameterizer::
derotatePoint(double theta, Coordinate_t &point) {

    double sintheta, costheta;
    Coordinate_t newpoint;

    sintheta = sin(theta);
    costheta = cos(theta);

    newpoint.x = point.x*costheta + point.y*sintheta;
    newpoint.y = -point.x*sintheta + point.y*costheta;
    
    point.x = newpoint.x;
    point.y = newpoint.y;
    
}


