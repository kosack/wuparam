//
// Program to do array analysis for multi-telescope data.  
// 
// 050910, KPK - right now it's just a stub of a program which does
// nothing in particular other than read the data.  Vicky will be
// doing the implementation.

#include <iostream>
#include <iomanip>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <exception>
#include <stdexcept>
#include <unistd.h>
#include <signal.h>

#include "Types.h"
#include "Exceptions.h"
#include "Config.h"
#include "Parameterizer.h"
#include "ProgressBar.h"
#include "ImageAnalyzer.h"
#include "Cutter.h"
#include "Camera.h"
#include "Log.h"
#include "ParamDataReader.h"
#include "PlotMaker.h"

const char ESC = 0x1b;

using namespace std;

/**
 * Function prototypes and globals for wuarray. I put these in their
 * own namespace so they don't conflict with functions in wuparam,
 * wucut, etc.
 */
namespace WUArray {

    /**
     * Describes the type of run that is being processed
     */
    enum RunType {
	TRACK,
	ON,
	OFF
    };

    void showUsage();
    void getCommandLineOptions(int argc , char **argv);
    void shutdown();
    void trap_ctrl_c(int);
    void process( RunInfo &ri );
    void analyzeRun( RunInfo &ri, RunType runtype );

    PlotMaker *pm=NULL;    

};

using namespace WUArray;


/**
 * Process the config file and do the analysis.
 */
int main( int argc, char* argv[] ) {


    // ===================================================================
    // Set up some things for handling ctrl-c and proper shutdown:
    // ===================================================================

    // set up ctrl-c handler  
    struct sigaction sa;
    
    sigfillset( &sa.sa_mask );
    sa.sa_flags = 0;
    sa.sa_handler = trap_ctrl_c;
    if (sigaction(SIGINT, &sa, NULL) < 0) {
        perror( "sigaction SIGINT" );
        exit(1);
    }
    
    // turn off cursor because it looks nicer! Uses the VT100 escape
    // sequence for "cursor off"    
    cerr << ESC << "[25l";

    // call the shutdown function whenever an exit() is encountered
    atexit(shutdown);

    // set up the logger (which is used to write out important error messages)
    Logger *logger = Logger::instance();
    logger->info("<<< This is WUArray >>>");

    // ===================================================================
    // Process the command-line options:
    // ===================================================================

    getCommandLineOptions( argc, argv );

    if (argc-optind < 1) {
	showUsage();
	exit(1);
    }

    // ===================================================================
    // Read the configuration and process the data! This is all
    // wrapped in a "try/catch block to handle any errors that get
    // thrown. Most wuparam functions throw an AnalysisException (or a
    // subtype like CriticalAnalysisException or
    // MildAnalysisException) when they encounter an error.  These are
    // caught and processing stops. I also catch standard template
    // library runtime_error exceptions, just in case something throws
    // one.
    // ===================================================================
    
    try {
	
	// ---------------------------------------------------------------
	// Open and read the config file...
	// ---------------------------------------------------------------


	Config conf( argv[optind] );

	if (conf.getRunCount() == 0) {
	    throw AnalysisException("No runs were specified in "
				    "the config file!");
	}
	
	// ---------------------------------------------------------------
	// Loop through all the runs in the config file:
	// ---------------------------------------------------------------

	RunInfo ri;

	while (conf.isDone() == false) {

	    ri = conf.getNextRun();

	    // process the run:

	    process( ri ); 

	}

	



    }
    catch (AnalysisException &e) {
	cout << "*** Couldn't process configuration '" << argv[optind] 
	     << "' because:  " << e.what() 
	     << endl;
	    
	exit(1);	    
    }
    catch (runtime_error &e) {
	cout << "*** RUNTIME ERROR: "<< e.what() << endl;
	exit(1);
    }
    catch (...) {
	cout << "*** Caught unknown or unhandled exception!"<<endl;
	exit(1);
    }

    
    exit(0);  // note this automatically calls shutdown() 
    

}




/**
 * Process a single run. This is called for each run in the config file.
 */
void WUArray::process( RunInfo &ri ) {

    cout << "WUArray: processing  ";
    switch( ri.type ){
	
    case Config::TRACK:
	cout << "tracking run: " << ri.onid << endl; 
	analyzeRun( ri, TRACK );
	break;
    case Config::ONOFF:
	cout << "ON/OFF pair: " << ri.onid<<" / "<<ri.offid<<endl; 
	analyzeRun( ri, ON );
	analyzeRun( ri, OFF );
	break;
    default:
	cout << "UNKNOWN RUN TYPE "<<ri.type<<endl;

    }
    
    
    



}


/**
 * Analyze a single run.  This is called by process.  Right now, it
 * just reads the data and demonstrates some simple things you can do.
 *
 * \param ri the full run-info struct \param runtype - should be
 * process the on, off, or tracking run? These cases may be treated
 * separately.
 */
void
WUArray::analyzeRun( RunInfo &ri, WUArray::RunType runtype ) {

    string id;

    if (runtype == TRACK || runtype == ON)
	id = ri.onid;
    if (runtype == OFF)
	id = ri.offid;
	        

    //===============================================================
    // Open the data file for reading:
    //===============================================================

    ParamDataReader reader(id+"/"+id+"-param.ntuple");
    HeaderRecord header;

    reader.getHeaderRecord( header );
    
    // the header contains some info about the parameterized data (see
    // DataRecords.h) For example:
    cout << "DATA HEADER: "<<endl
	 << "-------------------------"<<endl;
    cout << "SOURCENAME: " << header.sourcename << endl;
    cout << "NUM TELESCOPES: " << header.num_telescopes << endl;
    for (int i=0; i<header.num_telescopes;i++) {
	cout << "TELESCOPE "<< i<<" has "<< header.nadc[i] << " pixels"<<endl;
    }
    cout << endl;

    //===============================================================
    // Set up the array (which automatically loads the cameras)
    //===============================================================
    
    TelescopeArray array( ri.arrayconfig, header.nadc, atoi(ri.utdate.c_str()));

    cout << "There are "<<array.getNumTelescopes()<<" telescopes in the array"
	 << endl;
    cout << "The locations are: "<<endl;

    for (int i=0; i< array.getNumTelescopes(); i++) {
	
	cout << "\t(" << setw(6) << array.getLocationX(i) 
	     << ", "  << setw(6) << array.getLocationY(i) 
	     << ", "  << setw(6) << array.getLocationZ(i)  << ")"<<endl;

    }
    

    //===============================================================
    // for debugging purposes, you could plot the array, for
    // example. When constructing a plotmaker, the first argument is
    // the output type "X" for x-windows display, "ps" for
    // postscript. The second arg is the output file (ignored for
    // X-windows).  Note that I didn't really design plotmaker to work
    // nicely with array plots, so you may need to modify it to work
    // better.
    // ===============================================================

    if(pm==NULL) {
	pm = new PlotMaker( "X", "x.out" ,
			    array.getMinX()-10, array.getMinY()-10,
			    array.getMaxX()+10, array.getMaxY()+10 );
	cout << "DEBUG: Created new PlotMaker"<<endl;
    }

    cout << "DEBUG: plotting the array for example..."<<endl;
    pm->setScaleFactor( 20.0 );
    pm->setTicLevel( 2 );
    pm->setAxisBox( array.getMinX(), array.getMinY(),
			 array.getMaxX(), array.getMaxY() );
    pm->plot(array);
    pm->plotAxes();
    pm->plotTitle("The Array: ");
    pm->plotSubtitle("configuration: "+ri.arrayconfig);
    pm->newPage(); // newpage creates a new page for eps file, or
		   // swaps the display buffers for X plots
    sleep(2);

    // example of plotting a single camera:
    cout << "DEBUG: plotting camera 0 for example..."<<endl;
    Array_t testimage(array.getCamera(0)->getNumPixels());
    vector<Pedestal> testpeds(array.getCamera(0)->getNumPixels());
    vector<int> testcleanpix;
    Camera *cam = array.getCamera(0);

    pm->pushState(); // save the old space and settings
    pm->setTicLevel(1);
    pm->setScaleFactor(1);
    pm->setSpace(cam->getMinX()-0.4, cam->getMinY()-0.4,
		 cam->getMaxX()+0.4, cam->getMaxY()+0.4);
    pm->setAxisBox( *cam );
    pm->plot( *cam, testimage, testpeds,testcleanpix );
    pm->plotAxes();
    pm->plotTitle("Camera View");
    pm->plotSubtitle("Camera 0");
    pm->popState(); // restore old space and settings
    pm->newPage();
    sleep(2);

    // ===============================================================
    // Now, process the data in the data file: Here is where you'll
    // have to combine each sequence of events into an "array event",
    // which stores array info and an array of telescope-events
    // (HillasParameterizations).  You'll probably want to define a
    // structure like:
    //
    //     struct ArrayEvent {
    //    
    // 	     // Telescope Parameters:
    // 	     HillasParameterization televent[MAXTELS];
    //    
    // 	     //Array Parameters:
    // 	     ...
    //     }
    //    
    // to store the events in, and write a routine to write them to a
    // data file (e.g. using GSL's n-tuple writer)
    // ===============================================================

    HillasParameterization param;

    // here's an example of a 1D histogram from 0 to 90 with 18 bins:
    // you'll probably want to do an array or vector of these, one for
    // each telescope or something.
    Histogram examplehistogram(18, 0,90, "Alpha" );

    // An example of a 2D histogram: 30x30 bins, with coordinates
    // spanning (-2,-2)-(2,2)
    Image2D exampleimage( 30,30 );
    exampleimage.setCoordinateBox( -2,-2,2,2 );

    // a progress bar, just for fun:
    ProgressBar progress( reader.size() );

    cout << "There are "<< reader.size()
	 << " single-telescope events in the data"<<endl;

    // --------------------------------------------------
    // The main loop over events:
    // --------------------------------------------------

    for (int i=0; i<reader.size(); i++) {
	
	reader.getNextEventRecord( param );
	
	// note HillasParameterization contains a field called
	// "telescope_id" which is the integer telescope number
	// (starting at 0), which should be used to construct array
	// events.  The data should have events in the following order: 
	//       event 0 for telescope 0
	//       event 0 for telescope 1
	//       event 0 for telescope 2
	//       ..
	//       event 0 for telescope N
	//       event 1 for telescope 0
	//       event 1 for telescope 1
	//       event 1 for telescope 2
	//       .. 
	//       event 1 for telescope N
	//    ...
	
	examplehistogram.increment( param.alpha * 180.0/M_PI );
	exampleimage.addHistRadially( param.point_of_origin_a.x,
				      param.point_of_origin_a.y,
				      1.0, 0.2 );

	// print progress bar every 1024 events (don't do to often, as
	// it takes CPU cycles!)
	if (i%1024==0) progress.print( i ); 

    } // end loop over events

    // Done reading data.
    progress.printClear();

    // save the histogram to a file
    examplehistogram.save( id+"/"+id+"-example-alpha" );
    
    // save the image to a file. The extension will be .im2d, these
    // can be plotted later with wuplot or the image object itself
    // directly with the PlotMaker class.
    exampleimage.save( id+"/"+id+"-example-image" );
    
}

/**
 * Display a usage message:
 */
void
WUArray::showUsage() {

    cout << "USAGE: wuarray <config file>" << endl;
	
}


/**
 * This specifies and processes the various command-line options.  See
 * the getopt manpage for more information (type 'man 3 getopt' to see
 * the C documentation, which is in manual section 3).
 */
void 
WUArray::getCommandLineOptions(int argc , char **argv) {

    int ret=0;
    int c;

    while ( ret != -1 ) {

	ret=getopt( argc, argv, "ab:");
	c = ret;

	switch (c) {
	case 'a':
	    // recognises -a
	    cout << "Example of a switch option" <<endl;
	    break;
	case 'b':
	    // recognises -b <argument>
	    cout << "Example of a  option with argument: " 
		 << optarg <<endl;
	    break;
	default:
	    break;
	}
	
    }
    
}



void WUArray::trap_ctrl_c( int n) {

    Logger *logger = Logger::instance();

    if ( n == SIGINT ){
	cout << endl<<endl<<"WUArray was interrupted by CTRL-C. " << endl
	     << "Array analysis was not completed."<<endl;
	logger->info("WUArray was interrupted by CTRL-C");
    }

    exit(1);

}

void WUArray::shutdown() {

    Logger *logger = Logger::instance();

    // turn on cursor
    cerr << ESC << "[25h";

    // ---------------------------------------------------------------
    // clean up any global variables...
    // ---------------------------------------------------------------
    if (pm != NULL) {
	cout << "DEBUG: deleting plotmaker..."<<endl;
	delete pm;
    }
    
    logger->printf("WUArray exited");

}


