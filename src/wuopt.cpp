///////////////////////////////////////////////////////////////////////////
//  w  u  o p t
//
//  Washington University cut optimization program for Whipple
//  Gamma-Ray telescope data.
//
//  by Karl Kosack (2003-11-17)
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <queue>
#include <exception>
#include <stdexcept>
#include <iomanip>

#include "Types.h"
#include "Exceptions.h"
#include "Config.h"
#include "Parameterizer.h"
#include "ProgressBar.h"
#include "ImageAnalyzer.h"
#include "Cutter.h"
#include "Log.h"

using namespace std;
enum BoundType { UPPER, LOWER };

void getWUOptCommandLineOptions(int argc , char **argv);
void showUsage();

void optimizeParameter( Config *, ofstream &, string , BoundType );
void optimize2DParameter( Config *, ofstream &, string, BoundType,
			  int xpix, int ypix );


string What_to_optimize="lensize";
bool NewOutputDir = false;
string Outdir;
BoundType Bound;

double Low=0;
double High=0;
int Steps=0;

int
main( int argc , char **argv ) {

    Config *conf;
    int numruns=0;

    //-------------------------------------- 
    // Check command line args and output 
    // usage info or header

    cout << "wuopt - cut optimizer v"<< VERSION
	 << " <kosack@hbar.wustl.edu> " << endl
	 << "Washington University Physics Department, St. Louis M0"
	 << endl
	 << "============================================================" 
	 << endl << endl;

    getWUOptCommandLineOptions( argc, argv );
    
    //--------------------------------------
    // Open and read config file...
	
    try {
	conf = new Config( argv[6] );
    }
    catch (AnalysisException &e) {
	    
	cout << "*** Couldn't process configuration '" << argv[1] 
	     << "' because:  " << e.what()    << endl;
	    
	exit(1);
	    
    }


    try {
	cout << endl;
	cout << "Optimizing ";
	if (Bound == LOWER) cout << " lower ";
	else
	    cout << " upper ";
	cout << What_to_optimize<< " using database of " 
	     << conf->getRunCount() << " runs..." << endl;

	cout << "Optimization Limits: ("<<Low<<" , "<<High<<") in "
	     <<Steps<<" steps"<< endl;
       
	string bstring;
	if (Bound == LOWER) bstring = "-lower";
	else bstring = "-upper";
	string outfilename="optimize-"+What_to_optimize+bstring+".txt";
	ofstream outfile( outfilename.c_str() );
	
	cout << "Output to: "<< outfilename << endl;

	if (What_to_optimize[0] == '2')
	    optimize2DParameter( conf, outfile, What_to_optimize, Bound,
				 19,19 );
	else
	    optimizeParameter( conf, outfile, What_to_optimize, Bound );

	outfile.close();
	delete conf;

    }
    catch (AnalysisException &e) {
	cout << "WUOpt failed: " << e.what() << endl;	    
	cout << "Please correct problem and re-run wuopt." 
	     << endl;
    }
    



    return 0;


}


void 
optimizeParameter( Config *conf, ofstream &outfile, string parameter, BoundType bound ) {

    RunInfo ri;
    double delta = (High-Low)/(double)Steps;
    ProgressBar progress(Steps,"Optimizing");

    deque<RunInfo> &runqueue = conf->getRunQueue();
    deque<RunInfo>::iterator iter;   

    double value;
    RunStatistics stats;

    Cutter *cutter = new Cutter;
    cutter->enableOutput(false);
    cutter->enable2D(false);
    cutter->enableDisplay(false);

    //--------------------------------------
    // loop over optimization parameter
    
    for (int i=0; i<=Steps; i++) {    

	value = Low + delta*i;
	
	// go through each run in the queue...
	
	cutter->clear();
	for (iter = runqueue.begin(); iter != runqueue.end(); iter++) {
	    
	    progress.print(i);
	    
	    if (parameter=="length"){
		if (bound == LOWER) 
		    iter->cuts.length.lower = value;
		else 
		    iter->cuts.length.upper = value;
	    }
	    else if (parameter== "width"){
		if (bound == LOWER) 
		    iter->cuts.width.lower = value;
		else 
		    iter->cuts.width.upper = value;
	    }
	    
	    else if (parameter== "size"){
		if (bound == LOWER) 
		    iter->cuts.size.lower = value;
		else 
		    iter->cuts.size.upper = value;
	    }

	    else if (parameter== "distance"){
		if (bound == LOWER) 
		    iter->cuts.distance.lower = value;
		else 
		    iter->cuts.distance.upper = value;
	    }

	    else if (parameter== "max1"){
		if (bound == LOWER) 
		    iter->cuts.max1.lower = value;
		else 
		    iter->cuts.max1.upper = value;
	    }

	    else if (parameter== "max2"){
		if (bound == LOWER) 
		    iter->cuts.max2.lower = value;
		else 
		    iter->cuts.max2.upper = value;
	    }

	    else if (parameter== "lensize"){
		if (bound == LOWER) 
		    iter->cuts.lensize.lower = value;
		else 
		    iter->cuts.lensize.upper = value;
	    }

	    else if (parameter== "dlength"){
		iter->cuts.dlength = value;
	    }

	    else if (parameter== "dwidth"){
		iter->cuts.dwidth = value;
	    }

	    else if (parameter== "radialcut"){
		iter->cuts.smoothing_radius = value;
	    }

	    else if (parameter== "alpha"){
		if (bound == LOWER) 
		    iter->cuts.alpha.lower = value;
		else 
		    iter->cuts.alpha.upper = value;
	    }

	    else {
		throw AnalysisException("unknown or 2D parameter: "+parameter);
	    }
		
	    cutter->process( *iter );
	    
	}
	
	stats = cutter->getTotalPairStatistics();
	outfile << setw(10) << value <<" "
		<< setw(10) << stats.significance <<" "
		<< setw(10) <<stats.excess <<" "
		<<endl;
	outfile.flush();

    }
	
    progress.printClear();
    delete cutter;
    outfile.close();

}


void 
optimize2DParameter( Config *conf, ofstream &outfile, string parameter, 
		     BoundType bound, int xpix, int ypix ) {


    RunInfo ri;
    double delta = (High-Low)/(double)Steps;
    ProgressBar progress(Steps,"Optimizing");

    deque<RunInfo> &runqueue = conf->getRunQueue();
    deque<RunInfo>::iterator iter;   

    double value;
    RunStatistics stats;

    Image2D excess2d(39,39);
    Image2D signif2d(39,39);

    Cutter *cutter = new Cutter;
    cutter->enableOutput(false);
    cutter->enable2D(true);
    cutter->enableDisplay(false);

    cout << "DEBUG: parameter="<<parameter<<endl;

    //--------------------------------------
    // loop over optimization parameter
    
    for (int i=0; i<=Steps; i++) {    

	value = Low + delta*i;
	
	// go through each run in the queue...
	
	cutter->clear();
	for (iter = runqueue.begin(); iter != runqueue.end(); iter++) {
	    
	    progress.print(i);
	    
	    if (parameter=="2d-asymmdist"){
		if (bound == LOWER) 
		    iter->cuts.asymmdist.lower = value;
		else 
		    iter->cuts.asymmdist.upper = value;
	    }
	    else if (parameter=="2d-radsmooth"){
		iter->cuts.smoothing_radius = value;
	    }
	    else if (parameter=="2d-elongation"){
		iter->cuts.elongation = value;
	    }
	    else {
		throw AnalysisException("unknown parameter: "+parameter);
	    }
		
	    cutter->process( *iter );
	    
	}

	stats = cutter->getTotalPairStatistics();	
	excess2d.load( "Totals/excess.im2d" );
	signif2d.load( "Totals/signif.im2d" );
	outfile << setw(10) << value <<" "
		<< setw(10) << stats.significance <<" "
		<< setw(10) << stats.excess <<" "
		<< setw(10) << signif2d.getPixel(xpix,ypix) <<" "
		<< setw(10) << excess2d.getPixel(xpix,ypix) <<" "
		<<endl;
	outfile.flush();

    }
	
    progress.printClear();
    delete cutter;
    outfile.close();

}


void 
getWUOptCommandLineOptions(int argc , char **argv) {

    char c=0;
    string bound;


    if (argc!=7) {
	cout << "DEBUG: argc = "<< argc << endl;
	showUsage();
	exit(1);
    }

    What_to_optimize = argv[1];
    bound = argv[2];
    Low = atof(argv[3]);
    High = atof(argv[4]);
    Steps = atoi(argv[5]);


    if (bound == "lower") Bound = LOWER;
    if (bound == "upper") Bound = UPPER;
    else {
	cout << "unknown bound '"<<bound<<"', using LOWER instead";
	Bound = LOWER;
    }



}

void showUsage() {

    cerr << "USAGE: wuopt <parameter> <bound> <lowvalue> <highvalue> <steps> <config file>"<<endl
	 <<endl
	 << "\t<parameter>    parameter to optimze:"<<endl
	 << "\t               length, width, size, distance, max1, max2 "<<endl
	 << "\t               lensize, dwidth,dlength,radialcut, alpha"<<endl
	 << "\t               2d-asymmdist 2d-radsmooth 2d-elongation"<<endl
	 << "\t<bound>        lower or upper"<<endl
	 << "\t<lowvalue>     floating point lower value"<<endl
	 << "\t<highvalue>    floating point high value"<<endl
	 << "\t<steps>        number of steps"<<endl
	 << endl << endl
	 << "All options are specified in the config file."<<endl
	 <<endl;

}

