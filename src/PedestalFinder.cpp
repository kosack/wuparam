// PedestalFinder
// Karl Kosack <kosack@hbar.wustl.edu>

#include <iostream>
#include <fstream>
#include <string>
#include <valarray>
#include <list>
#include <vector>
#include <algorithm>

#include "Exceptions.h"
#include "DataReader.h"
#include "Camera.h"
#include "Config.h"
#include "ProgressBar.h"
#include "PedestalFinder.h"
#include "PlotMaker.h"
#include "Log.h"

using namespace std;

//#define DEBUG_PEDS

void
PedestalFinder::
getOnSourcePeds( RunInfo &ri, vector<Pedestal> &peds, int telescope_id ) {
    getPeds( ri,ri.onid,peds,telescope_id );
}

void 
PedestalFinder::
getOffSourcePeds( RunInfo &ri, vector<Pedestal> &peds, int telescope_id ) {
    getPeds( ri,ri.offid,peds,telescope_id );
}


/**
 * Set the thresholds to detect turned off tubes or stars.
 *
 * \param lower threshold to detect turned off tubes (default is 0.6)
 * \param upper threshold to detect stars (default is 1.5)
 */
void
PedestalFinder::
setTubeOffThresholds( double lower, double upper ) {

    _upthresh = upper;
    _lowthresh = lower;

}


/**
 * Returns the pedestal values as an array of Pedestal structs.
 *
 * \todo: Simulations only contain one pedestal record, so their
 * dispersions are always 0. Due to this, the ImageCleaner accepts too
 * many pixels, since it uses the picthresh*peddisp as the threshold!
 * Need to fix somehow...
 *
 * \todo: implement lockfiles to prevent more than one process from
 * writing/reading to the database (necessary for the MPI-enabled or a
 * multi-threaded version)
 */
void
PedestalFinder::getPeds( RunInfo &ri,  const string &id, 
			 vector<Pedestal> &peds, int telescope_id ) {

    register int i;
    RawDataReader *data = NULL;
    char telstring[50];
    sprintf( telstring, "%03d", telescope_id );
    string pedsfilename = ri.cachedir+id+"-"+telstring+".peds";
    ifstream infile;
    int nadc=0, goodnadc=0,n=0, count=0, progcount=0;
    Array_t peds2;
    int ntubesoff=0;
    RawHeaderRecord header;
    RawEventRecord event;
    list<int> startubes;
    list<int> offtubes;
    double d2;
    bool badeventflag;
    int badeventcount=0;
    char temp;
    vector<int> badpixel;
    int dtype;
    int windowsize;  //added charge integration window size, for calculating _badpedthresh.
    Logger *logger = Logger::instance();

    //===========================================================
    // Check if the pedestals already exist. 

    infile.open( pedsfilename.c_str(), ios::in  );
    if (infile) {

	//---------------------------------------------
	// READ IN PEDESTALS, SINCE THEY EXIST

	readPeds( pedsfilename, peds );

	infile.close();	
	return;
	
    }
    
    else { 

	//---------------------------------------------
	// CALCULATE PEDESTALS SINCE THEY DON'T EXIST:


	cout << "PedestalFinder: Calculating telescope "
	     <<telescope_id<<" pedestals for " << id 
	     << " in " << ri.cachedir 
	     << "..."<< endl;

	// Open the datafile

	RawDataReaderFactory *datafactory = RawDataReaderFactory::instance();
	
	try {
	    
	    data = datafactory->getReader( ri, id );
	    
	}
	catch (AnalysisException &e) {
	    delete data;
	    throw (e);
	} 
	
	int datasize = data->size();
	ProgressBar progress( datasize, "Pedestal" );

	
	dtype = data->getType();
	data->getHeaderRecord( header );
	if (telescope_id >= header.num_telescopes ) {
	    cout << "PedestalFinder: WARNING: telescope "<<telescope_id
		 << " is not in the data file's header!"<<endl;
	    cout << "Since the number is pixels is unknown, I'll try 490."
		 << endl;
	    nadc = 490;
	}
	else {
	    nadc = header.nadc[telescope_id];
	}

	windowsize = header.windowsize;

	goodnadc = nadc;
	if (nadc == 492) goodnadc = 379; // hack to not include outer tubes!
	if (nadc == 120) goodnadc = 109; // hack for 109 camera+extra tubes


	if (nadc == 499) _badpedthresh = 30*windowsize; 
	// hack to increase bad ped threshold for T1 camera.
	


	peds.resize(goodnadc);
	peds2.resize(goodnadc);
	badpixel.resize(goodnadc);
	Pedestal zero = {0,0,0};
	
	for (i=0; i<goodnadc; i++){
	    peds[i] = zero;
	    peds2[i] = 0.0;
	    badpixel[i] = 0;
	}


#ifdef DEBUG_PEDS	    		
	// for debugging:
 	//if (nadc == 499) Camera cam( "veritas_499.camera" );
	Camera cam( "490.camera" );
	vector<int> cleanpixels;
 	vector<Pedestal> dummypeds;
 	dummypeds.resize( nadc );
 	for ( i=0; i<nadc; i++){
 	    dummypeds[i].pedestal = 0.0;
 	    dummypeds[i].dispersion = 0.0;
 	    dummypeds[i].type = Pedestal::GOOD;
 	}
 
 	PlotMaker *pm  = new PlotMaker( "X", "x.out",
 					cam.getMinX() - 0.4, 
 					cam.getMinY() - 0.4, 
 					cam.getMaxX() + 0.4, 
 					cam.getMaxY() + 0.4 
 					);
 	pm->setAxisBox( cam );
 	pm->setPixelScale( 100 );
	// end
#endif

	// Go through each event in data.  If event is a pedestal,
	// accumulate it...

	while (data->isDone() == false) {

	    try {
		data->getNextEventRecord( event );
	    }
	    catch (EOFException e) {
		break;
	    }

#ifdef DEBUG_PEDS	    
	    // display for debugging:
 	    pm->plot( cam, event.adc, dummypeds, cleanpixels );
#endif	    

	    if (event.type == RawDataReader::PEDESTAL 
		&& event.telescope_id==telescope_id) {

#ifdef DEBUG_PEDS	    
 		pm->plotAxes();
#endif

		// check if any pixel in the event exceeds a
		// threshold, if so skip it! (this gets rid of
		// accidental nitrogen pulses that slip into the
		// pedestal events
		
		badeventflag=false;

		for (i=0; i<goodnadc; i++) {

		    if (event.adc[i] > _badpedthresh){
			badeventflag = true;
			badpixel[i]++;
		    }

		}

		// Accumulate pedestal if it's a good event

		if (badeventflag == false ){
		    for (i=0; i<peds.size(); i++) {
			peds[i].pedestal += event.adc[i];
			peds2[i] += (event.adc[i]*event.adc[i]);
		    }
		    n++;
		}
		else {
		    badeventcount++;
		}


	    }

	    // update progressbar
	    if (progcount == 512 || progcount == 0) { 
		progress.print( count );
		progcount = 0;
	    }
	    progcount++;

	    count++;

	    if (count > ri.max_ped_events) {
		cout << "PedestalFinder: bailed out after "<<ri.max_ped_events
		     <<" events"<< endl;
		break;
	    }


#ifdef DEBUG_PEDS	    
 	    cout << "DEBUG: ped event "<< n << " total events " 
 		 << count <<endl;
 	    pm->newPage();
 	    cout << "NEXT EVENT?";
 	    cin >> temp;
#endif	    

    	    
	}
	
	progress.printClear();

	// Output warning messages when there may be bad pixels
	cout << "PedestalFinder: events with pixels above 'bad' threshold: "
	     << badeventcount << endl;
	if (badeventcount>0){
	    for (i=0; i<goodnadc; i++) {
		if (badpixel[i] > 3) {
		    logger->printf("PedestalFinder (%s): pixel %d had high peds in %d events (it may be bad)", id.c_str(), i, badpixel[i]); 
		}
	    }
	}


	if (n>0) {
	    for (i=0; i<peds.size(); i++) {
		peds[i].pedestal /= (double)n;
		peds2[i] /= (double)n;
	    }
	}
	else {

	    logger->printf("PedestalFinder (%s): There were no pedestal"
			   " events in the data!", id.c_str());
	    
	}

	delete data;
	
	for (i=0; i<peds.size(); i++){
	    d2 = peds2[i] - pow(peds[i].pedestal,2);
	    peds[i].dispersion = d2>=0? sqrt(d2) : 0;
	}
	
	// now peds contains the correct pedestals and dispersions
	
	//----------------------------------------------------
	// Check for tubes which are off or malfunctioning:
	// (Simulations have no pedvars as of yet, so don't turn tubes
	// off for them)
    
	// First we want to calculate the median, 1/4 median, and 3/4
	// median points of the pedestal dispersions.  To do this, we
	// need to sort the list of dispersions.  upperlimit will be
	// the 3/4 median value times the _upthresh and lowerlimit is
	// the 1/4 median value times the _lowthresh

	double median, lowerlimit, upperlimit;
	vector<double> disps;
	int midpoint;
	disps.resize( peds.size() );
	for (i=0; i<peds.size(); i++)  disps[i] = peds[i].dispersion;
	sort( disps.begin(), disps.end() );

	midpoint = static_cast<int>(rint((static_cast<double>(disps.size())/2.0)));

	median = disps[midpoint];
	lowerlimit = median * _lowthresh;
	upperlimit = median * _upthresh;
	cout << "PedestalFinder: median pedestal dispersion is "
	     << median << endl;
	
	// Now, go through the tubes and mark them as good, or off

	if (dtype != RawDataReader::SIM ) {
	    cout << "PedestalFinder: threshold range = ("<<lowerlimit<< ","
		 <<upperlimit <<")" << endl;
	    
	    for (i=0; i<(int)peds.size(); i++ ) {
		
		if (peds[i].dispersion < lowerlimit) {
		    offtubes.push_back(i);
		    peds[i].type = Pedestal::TUBEOFF; // mark as off
		    ntubesoff++;
		}  
		else if (peds[i].dispersion >= upperlimit) {
		    startubes.push_back(i);
		    peds[i].type = Pedestal::STAR; // mark as star
		    ntubesoff++;
		}
		else {
		    peds[i].type == Pedestal::GOOD;
		}

	    }
	}

	// Print out the off tubes:
	
	startubes.sort();
	offtubes.sort();
	
	cout << endl;
	list<int>::iterator tube;
	cout << "The following tubes were shut off due to excess brightness"
	     << " (i.e. a star):"<<endl;
	cout << "-----------------------------------------"
	     << "-----------------------------" << endl;
	for (tube=startubes.begin(); tube != startubes.end(); tube++) {
	    cout << *tube << " ";
	}
	cout << endl << endl;
	cout << "The following tubes appear to have been turned off "
	     << "manually: "<<endl;
	cout << "-----------------------------------------"
	     << "-----------------------------" << endl;
	for (tube=offtubes.begin(); tube != offtubes.end(); tube++) {
	    cout << *tube << " ";
	}
	cout << endl << endl;


	cout << endl<< "PedestalFinder: " << ntubesoff << " of " << nadc 
	     << " tubes were shut off." << endl;

	//---------------------------------------------
	// write pedestals to disk:
	
	writePeds(pedsfilename, peds);

	cout << "PedestalFinder: Done calculating pedestals for " 
	     << id << endl;
	cout << "PedestalFinder: total pedestal events: "<< n << endl;
	cout << endl;

    }

}


void
PedestalFinder::writePeds( string &filename, vector<Pedestal> &peds ) {

    ofstream outfile;

    outfile.open(filename.c_str(), ios::out);
    if (!outfile) {

	cout << "*** Couldn't open " << filename << " for writing!" << endl;

    }
    else {

	outfile << peds.size() << endl;
	
	for(int i=0; i<(int)peds.size(); i++) {
	    outfile << i 
		    << "\t" << peds[i].pedestal 
		    << "\t" << peds[i].dispersion 
		    << "\t" << peds[i].type
		    << endl;
	}
	
	outfile.close();

    }


}


void
PedestalFinder::readPeds( string &filename, vector<Pedestal> &peds ) {

    ifstream infile(filename.c_str());
    int n;
    int id;
    Logger *logger = Logger::instance();

    if (!infile) {
	
	// shouldn't ever get here since we already tried opening...
	cout << "PedestalFinder::readPeds() Couldn't read pedestal file: "
	     << filename<< endl;
	exit(1);

    }


    infile >> n;
    peds.resize(n);

    for(int i=0; i<n; i++) {
	
	infile >> id 
	       >> peds[i].pedestal 
	       >> peds[i].dispersion
	       >> peds[i].type;

	if (id != i) {
	    logger->printf("PedestalFinder (%s): tube %d seems to be "
			   "misnumbered in the cache file", 
			   filename.c_str(), id);
		
	}

    }
    
    infile.close();

}



