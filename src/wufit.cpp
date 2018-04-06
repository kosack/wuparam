///////////////////////////////////////////////////////////////////////////
//  w  u  f i t 
//
//  Washington University Spectral Fitting Package for Whipple & VERITAS
//  Gamma-Ray telescope data.
//
//  by Henric Krawczynski (2003-07-30)
///////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <exception>
#include <stdexcept>
#include <iomanip>


#include "Histogram.h"
#include "EnergySpectrum.h"
#include "Fitter.h"
#include "Exceptions.h"
#include "Config.h"

using namespace std;

int
main( int argc , char **argv ) {
    
    cout << "wufit - spectral fitting "<< VERSION
	 << " <krawcz@wuphys.wustl.edu> " << endl
	 << "Washington University Physics Department, St. Louis M0"
	 << endl
	 << "============================================================" 
	 << endl << endl;
    
    // Get Input Parameters
    // The Data Histograms Have to Be in the Directory 
    // "Totals/energy-total-on.hist" and
    // Totals/energy-total-off.hist
    
    if (argc<10)  {
	cout << "Usage  : wufit <configfile> "
	     << "N0 Gamma Time(minutes) dN0 Steps dGamma Steps onoff-ratio E_CUTOFF(TeV)" << endl;
	cout << "Example: wufit sample.conf "
	     << "3.3E-7 -2.49 10000. 1. 20 1 20 1.0 4.33"<<endl;
	exit (-1);
    }

    // NOTE: Time should really be the total livetime.

    char *conffile = argv[1];
    double N0     = atof(argv[2]);
    double Gamma  = atof(argv[3]);
    double Time   = atof(argv[4]) *  60.0 ;
    
    double dN0    = atof(argv[5]);
    int    n0Steps= atoi(argv[6]);    
    double dGamma = atof(argv[7]);
    int    gSteps = atoi(argv[8]);
    double onoffratio = atof(argv[9]); 
    double E0 = atof(argv[10]);

    // Load the configuration information...
    // just need to look at one run to get the global information.
    RunInfo ri;
    Config conf( conffile );
    if (!conf.isDone()) {
	ri = conf.getNextRun();
    }
    else {
	cout << "Please specify at least one run in the configuration file"
	     << endl << "otherwise, the setup and cuts cannot be determined"
	     << endl << endl;

	exit(1);
	
    }

    cout << "N0    Starting Value: " << N0      <<endl
	 << "Gamma Starting Value: " << Gamma   <<endl
	 << "Time           Value: " << Time    <<endl
	 << "Delta N0            : " << dN0     << endl
	 << "N0 Steps            : " << n0Steps << endl
	 << "Delta Gamma         : " << dGamma  << endl
	 << "Gamma Steps         : " << gSteps  << endl 
	 << "ON/OFF ratio        : " << onoffratio << endl 
	 << "Energy Cutoff       : " << E0 << endl
	 << endl;

    try {
    
	//
	// Get Data Histograms
	// 
	
	EnergySpectrum* energy_on = new EnergySpectrum(ri,"EnergyEstimator-on");
	EnergySpectrum* energy_off= new EnergySpectrum(ri,"EnergyEstimator-off");
	energy_on->load("Totals/energy-total-on.hist");
	energy_off->load("Totals/energy-total-off.hist");
	
	int i;
	double lower,upper;
	cout << "Reading Data Histograms"<<endl
	     << "  Bin    E_min     E_max         On       Off"<<endl;
	for (i=0;i<energy_on->numBins();i++) {
	    energy_on->getRange( i, lower,upper );
	    cout << setw(5) << i 
		 << setw(10) << pow(10,lower)
		 << setw(10) << pow(10,upper)
		 << setw(10) << energy_on->get(i) 
		 << setw(10) << energy_off->get(i) << endl;
	}
	cout << endl << endl;
	
	//
	// Define Fitter, That Loads Monte Carlo Data
	// From File "/data/Whipple/Simulations/mall.dat"
	//
	Fitter* fitter = new Fitter( ri );
	
	//
	// Start Actual Fitting Process
	// Save Chi-Square Table in File "Totals/chi.text"
	//
	fitter->minimize(energy_on,energy_off,
			 N0,Gamma,E0,Time,
			 dN0,n0Steps,dGamma,gSteps, onoffratio);
	
	//
	// Print Out Results and Save Reconstructed Energy Spectra in
	// "Totals/espectrum.text"
	//
	fitter->save();
	
	delete energy_on;
	delete energy_off;
	delete fitter;
    }
    catch (CriticalAnalysisException &e) {
	cout << "Critical: "<< e.what() << endl;
    }
    catch (AnalysisException &e) {
	cout << "Exception: "<< e.what() << endl;
    }
    
	
}
