//
// wu2txt  - convert WUParam ntuples to text format
//

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
#include "ImageAnalyzer.h"
#include "ParamDataReader.h"
#include "Cutter.h"
#include "Log.h"

using namespace std;

namespace WU2Txt {
    void getCommandLineOptions(int argc , char **argv);
    void showUsage();
}

int
main( int argc , char **argv ) {

    WU2Txt::getCommandLineOptions(argc,argv);
    
    ParamDataReader *reader;
    HillasParameterization param;

    if (argc-optind < 1) {
	WU2Txt::showUsage();
	exit(1);
    }

    string ntuplefilename;

    try {

	for (int i=0; i<argc-optind; i++) {
	    
	    try {
		
		ntuplefilename = argv[optind+i];
		reader = new ParamDataReader( ntuplefilename );
		
		while (1) {
		    reader->getNextEventRecord( param );
		    cout << param << endl;
		}
	       		
	    }
	    catch (EOFException &e){
		// end of file
	    }
	}
    }
    catch (AnalysisException &e){
	cout << "conversion failed: "<< e.what() << endl;
    }
    
}



void 
WU2Txt::getCommandLineOptions(int argc , char **argv) {

    char c=0;

    while ( c != -1 && c!= 255) {

	c=getopt( argc, argv, "");

	switch (c) {

	default:
	    break;
	}
	
    }

}

void WU2Txt::showUsage() {

    cerr << "USAGE: wu2txt <ntuple> [<ntuple> ...]" << endl
	 << "output is to standard out by default"<< endl <<endl
	 << "\t-n output the number of events in the ntuple file" << endl;
    
}


