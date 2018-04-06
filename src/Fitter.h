// Fitter.h
// Henric Krawczynski <krawcz@wuphys.wustl.edu>
// Modifications by Karl Kosack <kosack@hbar.wustl.edu>
//

/**
 * \todo: Notes from KPK:
 *
 *  This needs to be generalized as follows:
 * 
 * - there should be a generic "model function" object (functor?) that
 *      contains a vector of parameters and a method which returns the
 *     value of the model given an energy. This gets used wherever the
 *     model function is needed.
 *
 * - Then, a function to minimize (chisqr) should be set up which is a
 *     gsl_multimin_f, and takes the vector of parameters as an input
 *     (the X vector as described in the GSL manual). Given a set of
 *     parameters, this function should fill the simulated histogram
 *     using the weighting factor (to which the model function object
 *     is passed), and compute the value of ChiSqr between the real
 *     and simulated histograms.
 *
 * - Finally, GSL's multimin solver gets called (using the simplex
 *     method, which doesn't require a derivative function) to find
 *     the best ChiSqr value.
 * 
 * - when the output is generated at the end, the modelfunction object
     is used again to scale the values. 
 *
 */

#ifndef FITTER_H
#define FITTER_H

//
// SPECTRAL FITTING PACKAGE
//

#include <vector>
#include "Histogram.h"
#include "Config.h"
using namespace std;


/**
 * For Each Monte Carlo Event That Passes the Selection Cuts, We Want to
 * Know the True Energy and the Reconstructed Energy 
 */
struct MCRecord {
    double trueEnergy_in_TeV;
    double estimatedEnergy;
};


/**
 * We store the data in an Array, We only want to know the number of
 * excess counts, and the variance of this number.
 */
struct Data{
    double onOff,onOffVariance;
};

/**
 * Here We Store The Results of the Chi-Square Fit
 */
struct Results {
    double chiSquare;
    int    dof;
    double n0,n0P,n0M;
    double gamma0,gammaP,gammaM;
    double e0;
};

struct MCBinInfo {

    double start;
    double end;
    double radius;
    int nummc;

};

/**
 * Object that takes care of all the Monte Carlo (MC) dependent aspects:
 *
 * - it stores the MC events
 *
 * - given model parameters, it can produce histograms with
 *   "artficial" energy spectra that can be compared to the data
 *
 * - given a pointer to the data, and a filled "artificial" energy spectrum,
 *   it will compute the chi-square value.
 *
 */
class MCSpectrum {

 public:

    MCSpectrum( RunInfo &ri );
    ~MCSpectrum();
    
    void   fill(double N0,  double Gamma0, double E0,  double Time);

    // uses class-histogram _spectrum and pointer to the experimental data
    // to determine chi-square value
    double chiSquare(std::vector<Data> &);
    
    inline int    dof();
    
    double get(int);
    void   getRange( int bin, double &lower, double &upper );
    
    static const int NUMMCBINS=9;

    EnergySpectrum* getSpectrum() { return _spectrum; }

 private:

    double weight( double energy,  double N0, 
		   double Gamma0, double E0, double time);
    
    void   info  ( double energy,
		  double &numShower, double &de, double &a, 
		  double &e1, double &e2);
    

    std::vector<MCRecord> _event;

    EnergySpectrum* _spectrum;     //  histogram that holds
				   // "artificial" energy spectrum

    
    std::string _dbfilename;
    
    // Class-variable holding the degrees of freedom of the fit
    // (depends on # of bins, and on event statistics)
    int _ndof;

    std::vector<MCBinInfo> _info;

};   

/**
 * returns degrees of freedom of the last chi-square fit performed,
 *  stored in class-variable _ndof.
 */
int 
MCSpectrum::dof() {
    return _ndof;
}

/**
 * Return histogram value at index i.
 */
inline double 
MCSpectrum::get(int i) {
    return _spectrum->get(i);
}

/**
 * finds the lower an upper range limits of the specified bin.
 */
inline void 
MCSpectrum::getRange( int bin, double &lower, double &upper ) {
  _spectrum->getRange(bin,lower,upper);
}

/**
 * Class that performs actual chi-square grid search.
 */
class Fitter {

 public:
    Fitter( RunInfo &ri ){

	_mcSpectrum = new MCSpectrum( ri );

    }
    ~Fitter(){;}
    
    //
    // Start Actual Fitting Process
    // Save Chi-Square Table in File "Totals/chi.text"
    //
    void minimize(EnergySpectrum* energy_on,
		  EnergySpectrum* energy_off, 
		  double &N0, double &Gamma0, double &e0, double Time, 
		  double n0Range, int n0Step, 
		  double gRange, int gStep,double onoffratio);
    
    //
    // Print Out Results and Save Reconstructed Energy Spectra in
    // "Totals/espectrum.text"
    //
    void   save();
    
 private:
    // Object that takes care of the MC data
    MCSpectrum *_mcSpectrum;
    
    // Array with experimental data
    std::vector<Data>  _data;
    
    // Results of chi-square fit
    Results    _results;
};

#endif
    
