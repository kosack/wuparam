// Histogram Class
// Karl Kosack <kosack@hbar.wustl.edu>


#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_histogram.h>
#include <vector>
#include <string>

/**
 * General 1-D histogram class.
 *
 * \todo: Make this class work with lightcurve data where the first
 * and last bins may be partially filled (or smaller width)
 */
class Histogram {
    
 public:
    
    Histogram(int,double,double, std::string);
    Histogram( const Histogram &h );
    ~Histogram();
   
    void   increment(double);
    void   accumulate(double x, double weight);
    int    maxBin();
    int    minBin();
    double maxValue();
    double minValue();
    void   getRange( int bin, double &lower, double &upper );
    double mean();
    double sigma();
    double sum();
    void   save(std::string);
    void   load(std::string);
    int    findBinWithValue( double val );
    void   saveStatistics( std::string );
    int    numBins() { return gsl_histogram_bins( _hist ); }  
    void   scale( double val );
    void   subtract( Histogram &h );
    double operator[](int);
    double get(int bin);
    void   reset();

    gsl_histogram* getGSLHist() { return _hist; }

 private:


    std::string _name;
    long int _overflow;
    long int _underflow;
    gsl_histogram *_hist;
    bool _is_saved;


};


/**
 * Increment the bin corresponding to the specified value. Implemented
 * as an inline function for speed.
 */
inline void
Histogram::increment(double x) {
    
    if (x>maxValue())      _overflow++;
    else if (x<minValue()) _underflow++;
    else {
	
	gsl_histogram_increment( _hist, x );

    }
}

/**
 * similar to increment(), but allows one to increment by the
 * specified weight instead of by 1.
 */
inline void
Histogram::accumulate(double x, double weight) {
    if (x>maxValue())      _overflow++;
    else if (x<minValue()) _underflow++;
    else {
      gsl_histogram_accumulate( _hist, x, weight );
    }
}


/**
 * Return histogram value at index i.
 */
inline double
Histogram::operator[](int i) {
    return gsl_histogram_get( _hist, i );
}

/**
 * Return histogram value at index i.  Can also use the [] operator
 * (i.e. myhist[3] is equivalent to myhist.get(3))
 */
inline double
Histogram::get(int i) {
    return gsl_histogram_get( _hist, i );
}



#endif
