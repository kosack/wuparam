//
// Fitter.cpp
// Henric Krawczynski <krawcz@wuphys.wustl.edu>
//

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <list>
#include <cmath>
#include <unistd.h>
#include <cstring>
#include <sstream>

#include "SuperCutter.h"
#include "ZCutter.h"
#include "SpectralCutter.h"
#include "Config.h"
#include "Histogram.h"
#include "EnergySpectrum.h"
#include "FitterNew.h"
#include "Exceptions.h"
#include "ProgressBar.h"
#include "ImageAnalyzer.h"
#include "DataRecords.h"
#include "DataWriter.h"
#include "ParamDataReader.h"

//#define ES_DEBUG 0

using namespace std;
using namespace Cuts;

// ++++++++++++++++++
// + MCdata Methods +
// ++++++++++++++++++

/**
 * Creator, gets all the MC information
 *
 * \todo: make the specific energy estimator function user selectable
 * (LZA vs small zenith, etc.)
 */
MCSpectrum::MCSpectrum( RunInfo &ri ) {
   
    _spectrum = new EnergySpectrum(ri,"MC-Spectrum");

    // ======================================================================
    // open the  MC Events database
    

    ifstream simfile;
    simfile.open(string(ri.mcdatabase+".simhdr").c_str());
  
    if (simfile.fail())  {
	cout << "Can Not Read MC File" << endl;
	throw CriticalAnalysisException("Couldn't open simulation database: '"
					+ ri.mcdatabase + "'" );
    }
    else
	cout << "Reading Monte Carlo Data Set: " << ri.mcdatabase<<endl;
  
    double x;
    MCRecord temprecord;

    cout << "with the following cuts:" << endl;
    cout << ri.cuts << endl << endl;
    
    HillasParameterization param,unscaled_param;
    int numsims=0;
    ProgressBar prog(0,"Reading Sims");

    // ======================================================================
    // Load the simulated events

    // first load the header from the file

    string line;
    MCBinInfo tmpbin;

    cout  << "   E_start     E_end    Radius    N_sims" << endl;

    while (getline(simfile,line) && line != "#===END===") {

	if ( line.substr(0,2) == "#>" ) {

	    vector<string> token;	    
	    tokenize( line, token );
	    
	    if (token.size() != 5) {
		cout << token.size() <<endl;
		cout << "UNKNOWN INFO IN HEADER: "<<line <<endl;
	    }
	    else {
		tmpbin.start = atof(token[1].c_str());
		tmpbin.end = atof(token[2].c_str());
		tmpbin.radius = atof(token[3].c_str());
		tmpbin.nummc = atoi(token[4].c_str());
		_info.push_back( tmpbin );
		cout << setw(10) << tmpbin.start 
		     << setw(10) << tmpbin.end 
		     << setw(10) << tmpbin.radius
		     << setw(10) << tmpbin.nummc
		     << endl;
	    }

	}

    }
    
    if (_info.size() != NUMMCBINS) {
	throw CriticalAnalysisException("Bad header in "
					+ri.mcdatabase);
    }


    cout << "Reading event data from " << ri.mcdatabase << endl;


    // Read the event data

    ParamDataReader reader( ri.mcdatabase );

    numsims=0;

    EnergyEstimatorFactory *efactory = EnergyEstimatorFactory::instance();
    EnergyEstimator *energy = efactory->getEstimator(ri.energyestimator);

    try {
	while (1) {

	    if(numsims++ % 100==0) prog.print(1);

	    reader.getNextEventRecord( param );
	    x = energy->getEstimate(param.size,param.distance,param.zenith);
	    temprecord.estimatedEnergy = x;
	    temprecord.trueEnergy_in_TeV = param.sim.primary_energy;
	    _event.push_back( temprecord );


	}
    }
    catch (EOFException &e) {
	cout << "DEBUG: end of file: "<<e.what() << endl;
    } 
       
    prog.printClear();

    simfile.close();

    cout << "Succesfully Read MC File, " << _event.size() <<" of "
	 << numsims << " events passed cuts" <<endl;


}


MCSpectrum::~MCSpectrum() {
    delete _spectrum;
}




/**
 * Given model parameters, fills "artificial" energy spectrum into
 * histogram _spectrum
 */
void 
MCSpectrum::
fill(double N0, double Gamma0, double E0, double Time) {
    
    double w;
    
    // Reset Histogram  
    _spectrum->reset();
    
    double n0 = pow(10.,N0);
    
    // Loop over MC events, and fill "artifical" energy spectrum
    vector<MCRecord>::iterator evt;
    for (evt = _event.begin(); evt != _event.end(); evt++){
	
	w = weight( evt->trueEnergy_in_TeV,n0,Gamma0,E0,Time );
	_spectrum->accumulate( evt->estimatedEnergy, w );
	
    }

}


/**
 * Computes weights for MC events, uses method "info"
 */
double MCSpectrum::weight(double energy, double N0, double Gamma0, 
			  double E0, double Time) {
    
    // that's the number of events per time area and energy that we
    // want to have
    double soll;
    if (E0>0){ 
	soll = N0*pow(energy,Gamma0)*exp(-1*energy/E0)*Time;
    }else{
	soll = N0*pow(energy,Gamma0)*Time;
    }
    double numShower,de,a,e1,e2;
    
    info(energy,numShower,de,a,e1,e2);
    
    if (numShower>0) {

	// here we compute the number of events per time area and
	// energy that we have simulated some calculus is needed, to
	// take into account that we simulated with differential
	// spectral index of 1.5.
	
	double gamma_sim,gamma_red,n0,ist;
	
	gamma_sim = -1.5;  // should be differential index
	// note some versions of Kascade incorrectly call the input
	// "integral" when it is actually differnetial!
	gamma_red =  gamma_sim+1.;

	// n0 is normalization constant of the simulated differential
	// energy spectrum
	n0        = -gamma_red * numShower/(pow(e1,gamma_red)-
					    pow(e2,gamma_red));

	ist       = n0 * pow(energy,gamma_sim) / a ;
	
       // ratio of what we want to have and what we do have is the
       // weighting factor
	
	return soll/ist;
    }
    else
	return 0.;

}
 
/**
 * Contains information abouyt # of simulated events in energy range
 */
void MCSpectrum::info(double energy,double &numShower, 
		      double &de, double &a, double &e1, double &e2){
    
    if ((energy<_info[0].start)||(energy > _info[NUMMCBINS-1].end)) {
	numShower = -1;
	return ;
    }

    int i;
    for (i=0;i<NUMMCBINS;i++)
	if (_info[i].end >= energy) break;

    numShower = _info[i].nummc;

    // width of energy interval
    de        = _info[i].end - _info[i].start;
    
    // get boundaries of simulation interval
    e1        = _info[i].start;
    e2        = _info[i].end;
    
    // area
    a         = M_PI * _info[i].radius * _info[i].radius;

 }

/**
 * uses histogram _spectrum and pointer to the experimental data
 * to determine chi-square value
 */

double MCSpectrum::chiSquare( std::vector<Data> &data ) {
    
    int i;
    double chiSquare = 0.;
    double diff;
    _ndof = -2;
    
#ifdef ES_DEBUG
    cout << "New Chi^2"<<endl;
#endif
    
    // loop over bins
    for (i=0;i<_spectrum->numBins();i++) {

	// Only use if sum of on and off counts exceeds MIN_COUNTS (5),
	// that gurantees nice chi-square statistics


	if (data[i].onOffVariance>EnergySpectrum::MIN_COUNTS) {
	    
	    //	 kpk: also only keep points where the value is >
	    //	 sqrt(variance), that way really sick points aren't counted.
	    //	    	   && data[i].onOff > sqrt(data[i].onOffVariance))  {

	   diff       = data[i].onOff - _spectrum->get(i);
	   chiSquare += diff*diff/(data[i].onOffVariance);
	   
#ifdef ES_DEBUG 
	   cout << setw(10) << i << " " 
		<< setw(10)<< data[i].onOff << " " 
		<< setw(10)<< data[i].onOffVariance << " " 
		<< setw(10)<< _spectrum->get(i) << " "  
		<< chiSquare<<endl; 
#endif
	   
	   _ndof++;
       }
    }
#ifdef ES_DEBUG
    cout << endl <<endl;
#endif
    return chiSquare;
}

// ++++++++++++++++++ 
// + Fitter Methods +
// ++++++++++++++++++

/**
 * Start Actual Fitting Process
 * Save Chi-Square Table in File "Totals/chi.text"
 *
 * \todo: put "alpha" or tracking ratio into the onOff and
 * onOffVariance: onoff = on-a*off, variance=a*(on+off) or whatever
 *
 */
void
Fitter::
minimize( EnergySpectrum* energy_on,EnergySpectrum* energy_off, 
	  double &N0, double &Gamma0, double &E0, double Time, 
	  double n0Range, int n0Step, 
	  double gRange, int gStep,
	  double e0Range, double e0Step, double onoffratio) {
    
    int flagE0;
    
    if(E0>0.0){
	flagE0 = 1;
    }else{ flagE0 = 0;}

    // check that the histograms sizes are correct:

    if (energy_on->numBins() != energy_off->numBins() ) {
	throw CriticalAnalysisException(string("Wrong energy spectrum sizes: ")
					+"rerun wucut" );
    }
    if (energy_on->numBins() != _mcSpectrum->getSpectrum()->numBins() ) {
	throw CriticalAnalysisException(string("energy spectrum and ") 
					+"monte carlo spectrum have "
					+ "different sizes: rerun wucut");
    }

    // 
    // Get information out of histograms
    //

    if (_data.size() != energy_on->numBins()) {
	_data.resize(energy_on->numBins());
    }

    int i;
    for (i=0;i<energy_on->numBins();i++){
	_data[i].onOff         = (energy_on->get(i)-
				  onoffratio*energy_off->get(i));
	_data[i].onOffVariance = onoffratio*(energy_on->get(i)
					     +energy_off->get(i));
    }
    
    // Now we loop 3 times over the data.
    // (1) loop over n0 and g to minimized the chi-square value
    // (2) loop over n0 to get error on n0
    // (3) loop over g to get error on g
    
    double n0,n1,n2,dn;
    n0 = log10(N0);
    n1 = n0-n0Range;
    n2 = n0+n0Range;
    dn = (n2-n1)/n0Step;
    
    double g,g1,g2,dg;
    g1 = Gamma0-gRange;
    g2 = Gamma0+gRange;
    dg = (g2-g1)/gStep;

    double e0,e1,e2,de;
    if (flagE0 == 1){
	    e0 = log10(E0);
	    e1 = e0-e0Range;
	    e2 = e0+e0Range;
	    de = (e2-e1)/e0Step;
    }

    double chiSquare;
    _results.chiSquare = 1E+32;
    _results.n0        = 0.;
    _results.gamma0    = 0;
    if (flagE0 == 1){
	_results.e0        = 0.;
    }else { _results.e0 = -1.;}
    
    ofstream chiFile("Totals/chi.text");
    // TODO: put in a progress bar here!
    ProgressBar prog(n0Step,"Fitting");
    int tmp=0;
    cout << "E0 " << "N0 " << "Gamma "  << "ChiSqr "  << endl;
    if(flagE0 == 1){
	for (e0=e1;e0<e2;e0+=de){ //loop over E0 to minimize the chisqr value
	    //cout  <<"E = " << pow(10,e0) ;
	    //prog.print(tmp);
	    tmp++;
	    // (1) loop over n0 and g to minimize the chi-square value
	    for (n0=n1;n0<n2;n0+=dn){
		
		// Loop over Gamma0
		for (g=g1;g<g2;g+=dg) {
		    
		    _mcSpectrum->fill(n0,g,e0,Time);
		    chiSquare = _mcSpectrum->chiSquare(_data);
		    //chiFile << pow(10,n0) << " " << g << " " << pow(10,e0) << " " << chiSquare << endl;
		    chiFile << pow(10,e0) << " " << chiSquare << endl;

		    if (chiSquare<_results.chiSquare) {
			_results.chiSquare = chiSquare;
			_results.n0        = n0;
			_results.gamma0    = g;
			_results.e0        = e0;
			cout << pow(10,e0) << " " << pow(10,n0) << " " << g << " " <<chiSquare << endl;
			//chiFile << pow(10,e0) << " " << chiSquare << endl;
#ifdef ES_DEBUG
			cerr << "BEST ChiSqr:  "<<chiSquare
			     << " at n0="<<n0<<", gamma0="<<g<<" E0= "<<e0<< "   " <<endl;
#endif
		    }
		}
		
	    }
	    
	}
    }else{ 
	// (1) loop over n0 and g to minimize the chi-square value
	for (n0=n1;n0<n2;n0+=dn){
	    
	    prog.print(tmp);
	    tmp++;
	    
	    // Loop over Gamma0
	    for (g=g1;g<g2;g+=dg) {
		
		_mcSpectrum->fill(n0,g,e0,Time);
		chiSquare = _mcSpectrum->chiSquare(_data);
		chiFile << n0 << " " << g << " " << e0 << " " << chiSquare << endl;
		if (chiSquare<_results.chiSquare) {
		    _results.chiSquare = chiSquare;
		    _results.n0        = n0;
		    _results.gamma0    = g;
		    _results.e0        = e0;
#ifdef ES_DEBUG
		    cerr << "BEST ChiSqr:  "<<chiSquare
			 << " at n0="<<n0<<", gamma0="<<g<<"  " <<endl;
#endif
		}
	    }
		
	}
    }



    chiFile.close();
    prog.printClear();
    cout << "Done" << endl;
    
    ProgressBar prog2( 0,"Error on n0");
    cout << "Getting Errors..." << endl;
    // (2) loop over n0 to get error on n0
    _results.n0P = _results.n0;
    _results.n0M = _results.n0;
    
    for (n0=n1;n0<n2;n0+=dn/10.) {

	prog2.print(1);

	_mcSpectrum->fill(n0,_results.gamma0,_results.e0,Time);
	chiSquare = _mcSpectrum->chiSquare(_data);
	if (chiSquare<(_results.chiSquare+1)) {
	    _results.n0P      = GSL_MAX_DBL( _results.n0P, n0 );
	    _results.n0M      = GSL_MIN_DBL( _results.n0M, n0 );
	    if (chiSquare<_results.chiSquare) 
		_results.n0 = n0;
	}
    }
    
    // (3) loop over g to get error on g
    _results.gammaP = _results.gamma0;
    _results.gammaM = _results.gamma0;
    
    prog2.setName("Error on gamma");
    for (g=g1;g<g2;g+=dg/10.) {
	prog2.print(1);
	_mcSpectrum->fill(_results.n0,g,_results.e0,Time);
	chiSquare = _mcSpectrum->chiSquare(_data);
	
	if (chiSquare<_results.chiSquare+1) {
	    _results.gammaP   = GSL_MAX_DBL( _results.gammaP, g );
	    _results.gammaM   = GSL_MIN_DBL( _results.gammaM, g );
	    if (chiSquare<_results.chiSquare) 
		_results.gamma0 = g;
	}

    }

    // (2) loop over e0 to get error on e0
    if (flagE0 == 1)
	{
	    _results.e0P = _results.e0;
	    _results.e0M = _results.e0;
	    
	    for (e0=e1;e0<e2;e0+=de/10.) {
		
		prog2.print(1);
		
		_mcSpectrum->fill(_results.n0,_results.gamma0,e0,Time);
		chiSquare = _mcSpectrum->chiSquare(_data);
		if (chiSquare<(_results.chiSquare+1)) {
		    _results.e0P      = GSL_MAX_DBL( _results.e0P, e0 );
		    _results.e0M      = GSL_MIN_DBL( _results.e0M, e0 );
		    if (chiSquare<_results.chiSquare) 
			_results.e0 = e0;
		}
	    }
	}

    prog2.printClear();
    cout << "Done" << endl;

    //
    // Fill the histogram with best fit energy
    //
    _mcSpectrum->fill(_results.n0,_results.gamma0,_results.e0,Time);

    // TODO: make a histogram of simulated alpha after cuts filled by
    // weighting according to best fit parameters similar to the
    // fill() function.


    //
    // Fill all results into variables
    // 
    _results.n0       = pow(10,_results.n0 );
    _results.n0P      = pow(10,_results.n0P);
    _results.n0M      = pow(10,_results.n0M);
    _results.dof = _mcSpectrum->dof();
    
    N0    = _results.n0;
    Gamma0 = _results.gamma0;

    if (flagE0 == 1){
	    _results.e0       = pow(10,_results.e0 );
	    _results.e0P      = pow(10,_results.e0P);
	    _results.e0M      = pow(10,_results.e0M);
    }else{
	    _results.e0       = -1.;
	    _results.e0P      = -1.;
	    _results.e0M      = -1.;
    }
}

/**
 * Print Out Results and Save Reconstructed Energy Spectra in
 * "Totals/espectrum.text"
 */
void Fitter::save() {
    
    // output _results
    // first give fit parameters:
    double lower, upper;
    cout << endl << endl << "Fit Results" << endl;
    
    cout << "(chi_s/dof) : ("
	 << setprecision(3) << _results.chiSquare << " / " 
	 << _results.dof    << ")"   << endl ;
    
    lower = _results.n0 - _results.n0M;
    upper = _results.n0P - _results.n0;
    cout << "N0          :  "<< setprecision(3)<< _results.n0        
	 << " - " << setprecision(3) << lower 
	 << " + " << setprecision(3) << upper << endl;
    
    lower = _results.gamma0 - _results.gammaM;
    upper = _results.gammaP - _results.gamma0;
    cout << "Gamma0      :  " << setprecision(3) << _results.gamma0  
	 << " - " << setprecision(3) << lower 
	 << " + " << setprecision(3) << upper << endl << endl;

    if (_results.e0>0)
	{
	    lower = _results.e0 - _results.e0M;
	    upper = _results.e0P - _results.e0;
	    cout << "E0          :  "<< setprecision(3)<< _results.e0        
		 << " - " << setprecision(3) << lower 
		 << " + " << setprecision(3) << upper << endl;
	}
    cout << endl << endl;
    cout << " N0 low_error  high_error gamma0 low_error  high_error E0 low_error  high_error "  << endl;
	cout << _results.n0 << "  " << _results.n0M << "  "  << _results.n0P << "  "  
	 << _results.gamma0 << "  " << _results.gammaM << "  " << _results.gammaP  << "  " 
	 << _results.e0  << "  " << _results.e0M << "  " << _results.e0P << endl;
    
    // Now give fluxes
    cout << "Writing Energy Spectrum to file." << endl;
        
    int i;
    double logEnergy,energy,lowerError,upperError,model;
    double observed,scaleFactor;
    double fluxModel,fluxObserved,error;
    ofstream fluxesFile("Totals/espectrum.text");
    
    fluxesFile << "# energy lowerError upperError data variance "
	       << "model flux error" << endl;

    for (i=0;i<_mcSpectrum->getSpectrum()->numBins();i++) {
	// column output
	// (1)-(3) info on energy bin,
	// (4) data counts
	// (5) error on (4)
	// (6) model counts
	
	_mcSpectrum->getRange(i,lower,upper );
	logEnergy   = (lower+upper)/2.;
	energy      = pow(10.,logEnergy);
	lowerError  = energy-pow(10.,lower);
	upperError  = pow(10.,upper)-energy;
	model         = _mcSpectrum->get(i);
	
	fluxesFile << setprecision(6)<< energy << " " 
		   << setprecision(6)<< lowerError << " " 
		   << setprecision(6)<< upperError << " "
		   << setw(5)         << _data[i].onOff << " " 
		   << setprecision(6)<< sqrt(_data[i].onOffVariance) << " " 
		   << setprecision(6)<< model << " ";
	
	if (model>0) {
	    
	    // strategy: compute model energy spectrum and scale it by
	    // ratio of predicted counts to observed counts
	    
	    observed      = _data[i].onOff;
	    scaleFactor   = observed/model;
	    fluxModel     = _results.n0*pow(energy,_results.gamma0)*exp(-energy/_results.e0);
	    fluxObserved  = scaleFactor * fluxModel;
	    
	    if (_data[i].onOffVariance>0)
		error = sqrt(_data[i].onOffVariance)/_data[i].onOff * fluxObserved;
	    else
		error = 0.0;
	}
	else {
	    fluxObserved =0.;
	    error=0.;
	}
	
	// column output 
	// (7) estimated flux for ith bin
	// (8) error on (7)
	fluxesFile << setprecision(6)<< fluxObserved << " " 
		   << setprecision(6)<< error << endl ;
    }

    fluxesFile.close();

}
