///////////////////////////////////////////////////////////////////////////
//  w  u  c  u  t
//
//  Washington University data cutting and analysis program for
//  Whipple Gamma-Ray telescope data.
//
//  by Karl Kosack (2002-06-01)
///////////////////////////////////////////////////////////////////////////

#include <iostream>
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
#include "Log.h"

using namespace std;

const char ESC = 0x1b;

void getCommandLineOptions(int argc , char **argv);
void showUsage();
void trap_ctrl_c( int n );
void shutdown();

bool Enable2d = true;
bool NewOutputDir = false;
bool PointingCheck = true;
bool Display = true;
bool Output =false;
bool ScaledParameterOutput = false;
bool FastImageProcessing=false;
string Outdir;
int TelescopeID = 0;


int
main( int argc , char **argv ) {

    RunInfo ri;
    Cutter cutter;
    Config *conf;
    int numruns=0;
    bool errorflag=false;

    //-------------------------------------- 
    // Check command line args and output 
    // usage info or header

    cout << "wucut - data cutter v"<< VERSION
	 << " <kosack@hbar.wustl.edu> " << endl
	 << "Washington University Physics Department, St. Louis M0"
	 << endl
	 << "============================================================" 
	 << endl << endl;


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


    getCommandLineOptions( argc, argv );

    if (argc-optind < 1) {
	showUsage();
	exit(1);
    }

    if (Enable2d)
	cutter.enable2D(true);
    
    cutter.enablePointingCheck( PointingCheck );
    cutter.enableDisplay( Display );
    cutter.enableOutput( Output );
    cutter.enableScaledParameterOutput( ScaledParameterOutput );
    cutter.enableFastImageProcessing( FastImageProcessing );
    cutter.setTelescopeID( TelescopeID );

    if (NewOutputDir) 
	cutter.setOutputDir( Outdir );

    Logger *logger = Logger::instance();
    logger->info("<<< This is WUCut >>>");


    //--------------------------------------
    // Open and read config file...
	
    try {
	conf = new Config( argv[optind] );
    }
    catch (AnalysisException &e) {
	    
	cout << "*** Couldn't process configuration '" << argv[1] 
	     << "' because:  " << e.what() 
	     << endl;
	    
	exit(1);
	    
    }

    cout << endl;

    cout << "Cutting " << conf->getRunCount() << " runs..." << endl;
    cout << "Using Telescope: " << TelescopeID << endl;
    ProgressBar progress(conf->getRunCount(),"Cutting");

    //--------------------------------------
    // go through each run in config file...
	    
    while (conf->isDone() == false) {

	ri = conf->getNextRun();

	try {
	    progress.setName(ri.onid);
	    progress.print(numruns);
	    numruns++;

	    cutter.process( ri );

	}
	catch (AnalysisException &e) {
	    cout << "Processing failed: " << e.what() << endl;	    
	    errorflag=true;
	    break;
	}
	
    }

    progress.printClear();
    delete conf;
    
    //-------------------------------------
    // If no runs were found, output a warning
	
    if (numruns == 0) {
	cout << "** No runs were specified in the configuration file: "
	     << argv[1] << endl
	     << "** Please edit it and add some runs to process!" << endl
	     << endl;	
    }
	
    cout << endl;
	
    //-------------------------------------
    // output statistics
    
    if (errorflag==false)
	cutter.outputStatistics();
    else 
	cout << "An error occured. Please correct it and re-run wucut." 
	     << endl;

    cout << "Goodbye." << endl;


    int ecount = logger->getErrorCount();
    if ( ecount > 0 ) {
	cout << "======================================" << endl;
	cout << "ATTENTION: "<<ecount<<" warnings were logged"<< endl;
	cout << "please look at '"<<logger->getFileName()<<"'"<<endl;
	cout << "======================================" << endl;
    }

    return 0;


}


void 
getCommandLineOptions(int argc , char **argv) {

    char c=0;

    while ( c != -1 && c!= 255) {

	c=getopt( argc, argv, "1d:PDoSft:");

	switch (c) {

	case '1':
	    Enable2d = false;
	    break;

	case 'd':
	    NewOutputDir = true;
	    Outdir = optarg;
	    break;

	case 'P':
	    PointingCheck = false;
	    break;

	case 'D':
	    Display = false;
	    break;
	    
	case 'o':
	    Output=true;
	    break;

	case 'S':
	    ScaledParameterOutput = true;
	    break;

	case 'f':
	    FastImageProcessing = true;
	    break;
	case 't':
	    TelescopeID = atoi(optarg);
	    break;

	default:
	    break;

	}
	
    }

}

void showUsage() {

    cerr << "USAGE: wucut [-1dPoDS] [-d <out dir>] [-t <tel num>] <config file>" << endl
	 <<endl
	 << "\t-1  1D analysis only (disable 2D)"<< endl
	 << "\t-d  <dir> output to the specified directory (default: Totals/)"
	 << endl
	 << "\t-P  disable the 2D pointing check (cross-correlation)"
	 << endl
	 << "\t-o  output raw events which pass cuts to an ntuple"
	 << endl
	 << "\t-S  if output is enabled, output scaled parameters instead of raw" 
	 << endl 
	 << "\t-D  disable the X display" 
	 << endl
	 << "\t-f  enable fast (but less accurate) 2D image processing"<<endl
	 << "\t-t <telescope number> "<<endl
	 << "\t    process only events from specified"
	 <<      " telescope, (default: 0)"
	 << endl << endl  
	 << "All options are specified in the config file."<<endl
	 <<endl;

}

void trap_ctrl_c( int n) {

    Logger *logger = Logger::instance();

    if ( n == SIGINT ){
	cout << endl<<endl<<"WUCut was interrupted by CTRL-C. " << endl
	     << "Cutting was not completed. If you re-run "
	     << "wucut, it will start up " << endl
	     << "where it left off." <<endl;
	logger->info("WUCut was interrupted by CTRL-C");
    }

    exit(1);

}

void shutdown() {

    Logger *logger = Logger::instance();

    // turn on cursor
    cerr << ESC << "[25h";

    logger->printf("WUCut exited");

}
