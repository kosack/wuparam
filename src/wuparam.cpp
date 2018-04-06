///////////////////////////////////////////////////////////////////////////
//  w  u  p  a  r  a  m 
//
//  Washington University data parameterization program for Whipple
//  Gamma-Ray telescope data.
//
//  by Karl Kosack (2000-08-21)
///////////////////////////////////////////////////////////////////////////

/**
 * \mainpage  wuparam - WU Whipple Data Analysis Package
 * 
 *  wuparam calculates the Hillas Parameters for events stored in
 *  Whipple telescope raw data files.
 *
 *  - To compile and run wuparam, you must install the following packages:
 *     -# The GNU Scientific Library (GSL)
 *          (http://rpmfind.net/linux/rpm2html/search.php?query=gsl-devel)
 *     -# GNU plotutils package
 *          (http://rpmfind.net/linux/rpm2html/search.php?query=plotutils)
 *
 *
 */

#include <iostream>
#include <ctime>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <list>
#include <exception>
#include <stdexcept>
#include <unistd.h>
#include <signal.h>

#ifdef HAVE_LIBMPI
#include <mpi.h>
#endif

#include "Types.h"
#include "Exceptions.h"
#include "Config.h"
#include "Parameterizer.h"
#include "DataReader.h"
#include "wuparam.h" 
#include "Log.h"

void showUsage();
void getCommandLineOptions(int argc , char **argv);
void trap_ctrl_c( int n );
void shutdown();

using namespace std;

const char ESC = 0x1b;

int
main( int argc , char **argv ) {

    RunInfo ri;
    Parameterizer param;
    Config *conf;
    int numruns=0;
    int myrank, numprocs;

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

    //-------------------------------------- 
    // Check command line args and output 
    // usage info or header

    cout << "===================================="
	 << "====================================" << endl
	 << "WUParam - parameterizer v"<< VERSION
	 << " <kosack@hbar.wustl.edu> " << endl
	 << "Washington University Physics Department, St. Louis M0"
	 << endl;

#ifdef HAVE_LIBMPI
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
    MPI_Comm_size( MPI_COMM_WORLD, &numprocs );
#endif

    cout << "===================================="
	 << "===================================="
	 << endl << endl;

    getCommandLineOptions( argc, argv );

    if (argc-optind < 1) {
	showUsage();
	exit(1);
    }


    param.enableDisplay( Options.display );
    param.enableInteractive( Options.interactive );
    param.enableVerbose( Options.verbose );
    param.enableOverwrite( Options.overwrite );
    param.enableCleanup( Options.clearcache );

#ifdef HAVE_LIBMPI
    cout << "USING MPI WITH " << numprocs << " PROCESSORS"
	 << endl;
#endif    

    Logger *logger = Logger::instance();

    logger->info("<<< This is WUParam >>>");

    //--------------------------------------
    // Open and read config file...
	
    try {
	conf = new Config( argv[optind] );
    }
    catch (AnalysisException &e) {
	    
	cout << "*** Couldn't process configuration '" << argv[optind] 
	     << "' because:  " << e.what() 
	     << endl;
	    
	exit(1);
	    
    }


    //--------------------------------------
    // go through each run in config file...

    try {
	    
	while (conf->isDone() == false) {
	    ri = conf->getNextRun();
	    numruns++;
#ifdef HAVE_LIBMPI
	    if (numruns%numprocs == myrank) {
		cout <<endl<< "Rank " <<myrank<< ": processing job "
		     <<numruns<<" of "<< numruns+conf->getRunCount()
		     << endl;
		param.process( ri );
	    }
#else
	    cout << endl << "[Processing run "<<numruns
		 <<", "<<conf->getRunCount()
		 << " left]"<< endl;
	    param.process( ri );
#endif
	    

	}
	
    }
    catch (AnalysisException &e) {
	cout << "Processing failed: " << e.what() << endl;	    
	exit(1);
    }

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
    cout << "Wuparam exited normally. Goodbye!" << endl;

#ifdef HAVE_LIBMPI
    MPI_Finalize();
#endif

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

    int ret=0;
    char c;

    Options.verbose = true;
    Options.interactive = false;
    Options.overwrite = false;
    Options.clearcache = true;


    while ( ret != -1 ) {

	ret=getopt( argc, argv, "qidomk");
	c = (char)ret;

	switch (c) {
	
	case 'q':
	      Options.verbose = false;
	    break;
	    
	case 'i':
	    Options.interactive = true;
	    Options.display = true;
	    break;
	    
	case 'd':
	    Options.display=true;
	    Options.interactive=false;
	    break;

	case 'o':
	    Options.overwrite=true;
	    break;
	    
	case 'k':
	    Options.clearcache = false;
	    break;

	default:
	    break;

	}
	
    }

}

void showUsage() {

    cerr << "USAGE: wuparam [-idqok] <config file>" << endl
	 << "\t-i\tInteractive mode" <<endl
	 << "\t-d\tDisplay events as they are processed (non-interactive)" 
	 <<endl
	 << "\t-q\tQuiet (less verbose text output)" <<endl
	 << "\t-o\tOverwrite previously parameterized files "
	 << "(skipped otherwise)" <<endl
	 << "\t-k\tKeep uncompressed files in the cache (normally deleted)"
	 << endl;
    
    cout <<endl
	 << "All other options are specified in the config file."<<endl
	 <<endl;

}


void trap_ctrl_c( int n) {

    Logger *logger = Logger::instance();

    if ( n == SIGINT ){
	cout << endl<<endl<<"WUParam was interrupted by CTRL-C. " << endl
	     << "Parametrization was not completed. If you re-run "
	     << "wuparam, it will start up " << endl
	     << "where it left off. If you wish to re-parameterize "
	     << "all files, run with the "<< endl
	     << "'-o' option (overwrite)." << endl
	     <<endl;
	logger->info("WUParam was interrupted by CTRL-C");
    }

    exit(1);

}

void shutdown() {

    Logger *logger = Logger::instance();

    // turn on cursor
    cerr << ESC << "[25h";

    logger->printf("WUParam exited");

}
