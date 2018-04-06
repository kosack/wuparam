// Histogram.cpp
// Karl Kosack <kosack@hbar.wustl.edu>
// Modified 030617 by JB to fix possible memory leak due to round-off error
// Modified 030716 by PFR: changed Histogram::save to change range format (%20g to %20.9g)

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstdio>
#include <assert.h>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_randist.h>

#include "Histogram.h"
#include "Exceptions.h"
#include "Config.h" // for tokenize()

using namespace std;


/**
 * Create a new Histogram.
 * @param numbins number of bins
 * @param min minimum value
 * @param max maximum value
 * @param name description of the histogram
 */
Histogram::Histogram(int nbins,double minx,double maxx, 
			   std::string name)
    : _name(name), _is_saved(false)
{

    if (nbins <1) 
	throw CriticalAnalysisException("Histogram "+name+" has < 1 bin!");

    if (minx>=maxx) {
	cout << "DEBUG: min = "<<minx<<" max="<<maxx<<endl;
	throw CriticalAnalysisException("Histogram "+name+" has min>=max");
    }

    _hist = gsl_histogram_alloc( nbins );

    if (_hist == NULL) {

	cout << "HISTOGRAM ALLOCATION FAILED: "
	     << "nbins="<<nbins
	     << "min="<<minx
	     << "max="<<maxx << endl;
	    
	throw CriticalAnalysisException("Couldn't allocate memory for "
					+_name+" histogram");
	

    }

    gsl_histogram_set_ranges_uniform( _hist, minx, maxx );

}


/**
 * Copy constructor
 * \bug: doesn't work yet
 */ 
Histogram::Histogram( const Histogram &h ) {

    _name = h._name;
    _underflow = h._underflow;
    _overflow = h._overflow;
    
    //    _hist = gsl_histogram_clone( h.getGSLHist() );

}

Histogram::~Histogram() {

//     double s = sum();

    gsl_histogram_free( _hist );

    // just a temporary note to check whether I'm saving all the histograms...
//     if (_is_saved == false && s > 0.0) {
// 	cout << "Histogram: info: '"
// 	     <<_name<<"' was destroyed without being saved"
// 	     << endl;
//     }

}



/**
 * Clear the histogram
 */
void
Histogram::
reset() {
    gsl_histogram_reset (_hist);
}

/**
 * finds the lower an upper range limits of the specified bin.
 */
void
Histogram::
getRange( int bin, double &lower, double &upper ) {

    if ( bin>=0 && bin<(int)gsl_histogram_bins( _hist )) {
	gsl_histogram_get_range( _hist, bin, &lower, &upper );
	
    }
    else { 

	cout << "DEBUG: bin="<<bin<<" nbins="<<gsl_histogram_bins(_hist)
	     << endl;
	
	throw CriticalAnalysisException("Histogram '"+_name
					+"' getRange() index out of range!");
    }

}


/**
 * Return index of maximum bin
 */
int
Histogram::maxBin() {
    return gsl_histogram_max_bin( _hist );
}

/**
 * Return index of minimum bin
 */
int
Histogram::minBin() {
    return gsl_histogram_min_bin( _hist );
}

/**
 * Return index of maximum value
 */
double
Histogram::maxValue() {
    return gsl_histogram_max( _hist );
}

/**
 * Return index of minimum value
 */
double
Histogram::minValue() {
    return gsl_histogram_min( _hist );
}




/**
 * Output a tab delimited datafile and corresponding GNUplot file.
 * @param filename base filename name of output text file (file
 * extension will be added.
 *
 * We follow the convention that the x-value of the histogram is
 * the value of the lower bound on the corresponding bin.  To calculate
 * the value at the center of the bin it is necessary to add half of
 * the bin width.
 */
void
Histogram::save(string filename) {

    string outfilename = filename + ".hist";
    FILE *fp;

    fp = fopen( outfilename.c_str(), "w");
    gsl_histogram_fprintf( fp, _hist, "%20.14g ", "%20g " );  //Changed %20g to %20.14g  (PFR)
    fclose( fp );

    saveStatistics( filename );
    
    _is_saved = true;

}


/**
 * Write out a text file with statistical info about the histogram.
 * This function is automatically called when Histogram::save() is
 * called, but is provided separately in case one wants to just write
 * statistics and not the whole histogram.  It would be cleaner to
 * append this to the top of the histogram file itself, but that would
 * make it harder for external programs to read the histogram files,
 * so I decided to keep the stats file separate.
 *
 * \param basename base filename (.hist.stats is appended)
 */
void
Histogram::
saveStatistics( string basename ) {

    string outbasename = basename + ".hist.stats";
    ofstream outfile(outbasename.c_str());
    
    outfile << "NAME: " << _name << endl;
    outfile << "NBINS: "<< gsl_histogram_bins(_hist) << endl;
    outfile << "RANGE: "<< minValue() << ", "<< maxValue() << endl;
    outfile << "MEAN: " << mean() << endl;
    outfile << "SIGMA: " << sigma() << endl;
    outfile << "SUM: " << sum() << endl;
    outfile << "OVERFLOW: " << _overflow << endl;
    outfile << "UNDERFLOW: " << _underflow << endl;
    outfile << "MAX_VAL: " << maxValue() << endl;
    outfile << "MIN_VAL: " << minValue() << endl;
    outfile << "MAX_BIN: " << maxBin() << endl;
    outfile << "MIN_BIN: " << minBin() << endl;
    
    outfile.close();


}

/**
 * Read in a histogram saved with the "save" function.  The histogram
 * to be read must have the number of bins as the Histogram object. 
 *
 */
void   
Histogram::
load(std::string filename) {

    FILE *fp;
    ifstream statsfile;

    // Check for a corresponding .stats file, and load the settings
    // from it if it exists, otherwise just go with what the user
    // specified.
    
    statsfile.open( (filename+".stats").c_str() );
    if (! statsfile.fail()) {
	
	string line;
	int n=1;

	vector<string> tokens;
	vector<string> tok2;


	while (statsfile && !statsfile.eof()) {

	    tokens.clear();
		    
	    getline( statsfile, line, '\n' );
	    tokenize( line, tokens, ":");
	    
	    if (tokens.size() == 2) {
		
		if (tokens[0] == "NAME") {
		    _name = tokens[1];
		}
		else if (tokens[0] == "NBINS") {
		    int nbins = atoi(tokens[1].c_str());
		    if (nbins != this->numBins()) {
			cout << "Histogram '"<<_name<<"': resizing from "
			     << numBins() << " to "
			     << nbins << " entries (from .stats file)" <<endl;
			gsl_histogram_free(_hist);
			_hist = gsl_histogram_alloc(nbins);
		    }
		}
		else if (tokens[0] == "RANGE") {
		    double min,max;
		    tok2.clear();
		    tokenize( tokens[1], tok2, "," );
		    min = atof(tok2[0].c_str());
		    max = atof(tok2[1].c_str());
		    gsl_histogram_set_ranges_uniform( _hist, min,max);
		}
		else if (tokens[0] == "OVERFLOW" ) {
		    _overflow = atoi(tokens[0].c_str());
		}
		else if (tokens[0] == "UNDERFLOW" ) {
		    _underflow = atoi(tokens[1].c_str());
		}

	    }
	    
	    n++;
	}
	statsfile.close();
	
    }

    // Load the histogram data:

    fp = fopen( filename.c_str(), "r" );
    
    if (fp == NULL) {
	throw CriticalAnalysisException( "Couldn't open: "+filename );
    }
    
    if ( gsl_histogram_fscanf( fp, _hist ) ) {
	throw CriticalAnalysisException("Couldn't read histogram: "+filename );
    }
    
    fclose(fp);

}

/**
 * \returns the mean value in the histogram.
 */
double
Histogram::
mean() {
    return gsl_histogram_mean( _hist );
}


/**
 * \returns the variance in the histogrammed value. 
 */
double
Histogram::
sigma() {
    return gsl_histogram_sigma( _hist );
}


/**
 * \returns the total number of counts in the histogram.
 */
double
Histogram::
sum() {
    return gsl_histogram_sum( _hist );
}



/**
 * \returns index of bin with the specified value, -1 on failure.
 */
int    
Histogram::
findBinWithValue( double val ) {

    size_t bin;
    int ret;
    
    ret = gsl_histogram_find( _hist, val, &bin );
    
    if ( ret == GSL_SUCCESS) 
	return static_cast<int>(bin);
    else 
	return -1;

}


/**
 * Multiplies all the bins in the histogram by the specified value
 */
void
Histogram::
scale( double val ) {

    gsl_histogram_scale( _hist, val );

}


/**
 * Subtracts the specified histogram from the current one. 
 *
 * \param hist histogram containing values to subtract. Must have
 * identical bin ranges as the current histogram.
 */
void 
Histogram::
subtract( Histogram &hist ) {

    int status;
    gsl_histogram *otherhist = hist.getGSLHist();

    if (gsl_histogram_equal_bins_p( _hist, otherhist )) {

	status = gsl_histogram_sub( _hist, otherhist );
	
    }
    else {
	
	throw CriticalAnalysisException(string("Histogram::subtract(): ")
					+"histograms have differnent sizes, "
					+" cannot subtract!");

    }
    

}
