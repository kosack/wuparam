// GaussianDeviate.cpp
// Karl Kosack <kosack@hbar.wustl.edu>

#include <iostream>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <iomanip>
#include <list>
#include <gsl/gsl_rng.h>       /* Random num generators */
#include <gsl/gsl_randist.h>   /* Random distributions */
#include <gsl/gsl_sf_trig.h>

#include "GaussianDeviate.h"
#include "ProgressBar.h"
#include "Histogram.h"

using namespace std;

GaussianDeviate::GaussianDeviate( string &dbdir ) 
    : _dbfilename(dbdir+"gaussdev.dat"), _current(0)
{

    ifstream dbfile( _dbfilename.c_str(), ios::binary);

    // Set up random number generator

    gsl_rng_env_setup();		// set up random number system
    _rngtype = gsl_rng_taus;	        // faster, simulation quality.
    _rng = gsl_rng_alloc( _rngtype );	// allocate random number generator
    _seed = (unsigned long) time(NULL);	// seed to the current time
    gsl_rng_set( _rng, _seed );

    // set current to a random starting point
    _current = (int)(gsl_rng_uniform(_rng)*DBSIZE);

    // Attempt to read the database, or generate it if needed.
    if (dbfile.fail()) {
	cout << "Generating Gaussian Deviate database (this "
	     << "only needs to be done once). "<<endl
	     << "Please stand by..."<< endl;
	generate();
    }
    else {
	dbfile.read( (char*)_deviates, int(DBSIZE*sizeof(double)) );
    }

      
}




/**
 * Return a gaussian deviate
 */
double
GaussianDeviate::getFastDeviate() {
    _current++;
    if (_current>=DBSIZE) _current=0;
    return _deviates[_current];
}

void
GaussianDeviate::generate() {

    // Calculate deviates...
    
    ProgressBar progress(DBSIZE);
    int count=0;

    for (int i=0; i<DBSIZE; i++) {
	_deviates[i] = gsl_ran_ugaussian(_rng);

	if (count == 10000){
	    progress.print(i);
	    count=0;
	}
	count++;
    }
    progress.printClear();

    // save to database...

    ofstream outfile( _dbfilename.c_str(),ios::binary|ios::out);
    
    if (outfile.fail()){
	cerr << "Couldn't save Gaussian Deviate database: "<< _dbfilename
	     << endl;
	return;
    }

    outfile.write( (char*)_deviates, DBSIZE*sizeof(double) );

    // output a histogram for diagnostics
    
    Histogram testdev(100,-5,5,"Gaussian Deviate Diagnostic");
    for(int i=0; i<DBSIZE; i++) {
	testdev.increment(_deviates[i]);
    }
    //    testdev.fitGaussian( 0,99 );
    testdev.save(_dbfilename);

}
