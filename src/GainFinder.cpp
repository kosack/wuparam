//
// GainFinder.cpp
//
// Karl Kosack <kosack@hbar.wustl.edu>

#include <list>
#include <vector>

#include "Exceptions.h"
#include "DataReader.h"
#include "Config.h"
#include "PedestalFinder.h"
#include "ProgressBar.h"
#include "Camera.h"
#include "PlotMaker.h"

#include "GainFinder.h"

using namespace std;

/**
 * Set gain cutoff thresholds. The default values are used if this 
 * method is not called explicitly.
 * @param min - minimim adc value to include in gain calculation (default is 50)
 * @param max - maximum adc value to include in gain calculation (default 1000)
 * @param percent - percent of camera above which the image is considered "saturated"
 * and should be ignored for gain calculations (default 0.75)
 */
void
GainFinder::setThresholds( double min, double max, double percent ) {

    _minadc = min;
    _maxadc = max;
    _minpercent = percent;

}

/**
 * Fetch the gains for the specified nitrogen run.
 * @param ri - Run information struct from Config
 * @param id - nitrogen gain id
 * @param gains - uninitialized Array_t where gains will be stored
 *
 * \todo: make gains database have a code column which marks bad gains
 * \todo: implement lock file!
 */
void
GainFinder::getGains( RunInfo &ri, Array_t &gains, int telescope_id ) {


    PedestalFinder ped;
    RawDataReader *data = NULL;
    char telstring[50];
    sprintf( telstring, "%03d", telescope_id );
    string n2filename = ri.cachedir+ri.n2id+"-"+telstring+".gains";
    vector<Pedestal> peds;
    Array_t adc;
    Array_t ngains;
    double meanadc=0.0;
    int nadc = 0;
    int nlow=0, nhigh=0, nmean=0;
    RawHeaderRecord header;
    RawEventRecord event;
    list<int> badgains;
    list<int>::iterator tube;
    int n=0, cnt=0;
    int skipped = 0;
    string temp;
    
    //=======================================================
    // First, we need to fetch the pedestals for the N2 run need to
    // change the offset thresholds since N2 runs have very few events
    // and thus large variances.
    
    ped.setTubeOffThresholds( 0.3, 3.0 );
    ped.getPeds( ri, ri.n2id, peds ); 
    
    //=======================================================
    // Check if nitrogen gains exist in database...

    ifstream infile( n2filename.c_str() );
    if (infile) {
	// Gains exist, so read them in

	readGains( n2filename, gains );

	infile.close();
	return;

    }
    else {

	// Gains don't exist, so calculate them

	//---------------------------------------------
	// Open the datafile

	RawDataReaderFactory *datafactory = RawDataReaderFactory::instance();
	
	try {
	    
	    data = datafactory->getReader( ri, ri.n2id );
	    
	}
	catch (AnalysisException &e) {
	    delete data;
	    throw (e);
	} 

	cout << "GainFinder: Calculating telescope " << telescope_id 
	     << " gains for " << ri.n2id  
	     << " in " << ri.cachedir
	     <<"..." << endl;
	
	//======================================================= 
	// Check if this is simulation data and just return an array of
	// 1.0's (kind of a hack, so should be implemented better
	// later...)
	
	if (data->getType() == RawDataReader::SIM) {
	    cout << "GainFinder: simulation data, using 1.0 for all"<< endl;
	    gains.resize(490);
	    for(int i=0; i<(int)gains.size(); i++)
		gains[i] = 1.0;
	    writeGains( n2filename, gains );
	    delete data;
	    return;
	}
	

	ProgressBar progress( data->size(), "Gains" );

	//---------------------------------------------
	// Go through each event in data...

	data->getHeaderRecord( header );

	if (telescope_id >= header.num_telescopes ) {
	    cout << "GainFinder: WARNING: telescope "<<telescope_id
		 << " is not in the data file's header!"<<endl;
	    cout << "Since the number is pixels is unknown, I'll try 490."
		 << endl;
	    nadc = 490;
	}
	else {
	    nadc = header.nadc[telescope_id];
	}


	gains.resize(nadc);
	ngains.resize(nadc);
	adc.resize(nadc);
	
	for(int i=0; i<nadc; i++) {
	    ngains[i] = 0;
	    gains[i] = 0.0;
	}

	if (nadc == 492) nadc = 379; // hack to ignore outter tubes 

	// for debugging:
// 	Camera cam( "490.camera" );
// 	vector<int> cleanpixels;
// 	PlotMaker *pm  = new PlotMaker( "X", "x.out",
// 					cam.getMinX() - 0.4, 
// 					cam.getMinY() - 0.4, 
// 					cam.getMaxX() + 0.4, 
// 					cam.getMaxY() + 0.4 
// 					);
// 	pm->setAxisBox( cam );
// 	pm->setPixelScale( 1024 );
	

	while ( data->isDone() == false ) {

	    // get next event

	    try {
		data->getNextEventRecord( event );
	    }
	    catch (MildAnalysisException &e) {
		break;
	    }

	    
	    if ( event.type == RawDataReader::EVENT 
		 && event.telescope_id == telescope_id) {
		

		// debugging:
// 		pm->plot( cam, event.adc, peds, cleanpixels );
// 		pm->newPage();
// 		cout << "NEXT EVENT?";
// 		cin >> temp;


		// Subtract pedestals

		Array_t &a = event.adc;
		
		for (int i=0; i<nadc; i++) {
		    
		    //check for tube off
		    if (peds[i].type == Pedestal::GOOD ) {
			adc[i] = a[i] - peds[i].pedestal;
		    }
		    else adc[i] = 0.0;
		}


		// debugging:
// 		pm->plot( cam, adc, peds, cleanpixels );
// 		pm->newPage();
// 		cout << "NEXT EVENT?";
// 		cin >> temp;


		// Get the mean adc value

		meanadc = 0;
		nmean = 0;
		nlow = 0; 
		nhigh = 0;
		for (int i=0; i<nadc; i++) {
		    if (adc[i] < _minadc) nlow++;
		    else if (adc[i] > _maxadc) nhigh++;
		    else {
			meanadc += adc[i];
			nmean++;
		    }
		}
		
		// Accumulate the adc values

		if (nmean > (int)(nadc*_minpercent)) {
		    // Tubes are unsaturated, calculate
		    // the flat field values
		    
		    meanadc /= (double) nmean;

		    for (int i=0; i<nadc; i++) {
			if ((adc[i] > _minadc) && (adc[i] < _maxadc)) {
			    gains[i] += (meanadc / adc[i]);
			    ngains[i] += 1.0;
			}
		    }
		    	    
		}
		else {
		    // tubes are saturated or off
		    skipped++;
		}
	    }

	    // Update progressbar
	    n++;
	    cnt++;
	    if (cnt == 128 || cnt==0) {
		cnt = 0;
		progress.print( n );
	    }


	}

	delete data;

	progress.printClear();
	cout << "GainFinder: Skipped " << skipped 
	     << " events with too small a signal" << endl;

       
	//---------------------------------------------
	// now calculate the gains:

	for (int i=0; i<nadc; i++) {

	    // normalize gains
	    if (ngains[i] > 0.5) {
		gains[i] /= ngains[i];
	    }
	    else {
		// KPK: changed this from =1.0 to =0.0 to avoid
		// tubes with bad gains
		gains[i] = 0.0;
		badgains.push_back(i);
	    }
	}

	if (badgains.size() > 0) {
	    
	    cout <<endl;
	    cout << "There wasn't enough data to determine gains for the "
		 << "following tubes:" << endl;
	    cout << "-----------------------------------------"
		 << "-----------------------------" << endl;
	    for (tube=badgains.begin(); tube!=badgains.end(); tube++) 
		cout << *tube << " ";
	    cout << endl;
	}

	//--------------------------------------
	// Write the gains to disk

	writeGains( n2filename, gains );

	cout << "GainFinder: Done calculating gains for " 
	     << ri.n2id << endl;	

    }

}


void
GainFinder::readGains( string &filename, Array_t &gains ) {

    ifstream infile( filename.c_str() );
    int n, id;

    if (!infile) {
	
	// shouldn't ever get here since we already tried opening,
	// but just to be safe...
	cout << "GainFinder::readGains() Couldn't read gains file:"
	     << filename << endl;
	exit(1);

    }
    

    infile >> n;
    gains.resize(n);
    
    for (int i=0; i<n; i++) {

	infile >> id >> gains[i];
	
	if (id != i) {
	    cout << "\t*** WARNING: tube " << id 
		 << " seems to be misnumbered in " 
		 << filename << endl;    
	}

    }

    infile.close();

}

void
GainFinder::writeGains( string &filename, Array_t &gains ) {

    ofstream outfile( filename.c_str() );

    if (!outfile) {
	cout << "** Warning: Couldn't write to gains file: " 
	     << filename << endl;
    }
    else {
	
	outfile << gains.size() << endl;

	for (int i=0; i< (int)gains.size(); i++){
	    outfile << i << "\t" << gains[i] << endl;
	}

	outfile.close();
	
    }

}


