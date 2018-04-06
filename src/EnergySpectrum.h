// Energy Spectrum Histogram Class
// Henric Krawczynski <krawcz@wuphys.wustl.edu>

//
// Defines Standard Histogram for Determining Energy Spectra
//

#ifndef ENERGYSPECTRUM_H
#define ENERGYSPECTRUM_H
#include <cmath>
#include "Histogram.h"
#include "Config.h"


/**
 * Energy Spectrum 1-D histogram class.
 *
 */
class EnergySpectrum : public Histogram {

 public:
    

    // During chi-square fitting we use a bin only is the sum of 
    // on and off counts exceeds the canonical number 5.
    static const int MIN_COUNTS = 5;
    
    // the size of the EnergySpectrum comes from the config file RunInfo
    EnergySpectrum(RunInfo& ri,std::string name): 
	Histogram(ri.num_e_bins,log10(ri.emin_tev),log10(ri.emax_tev), 
		  name) {; }
    ~EnergySpectrum() {;}

}; 

struct EnergyEstimator {
    
    EnergyEstimator(){;}
    virtual double getEstimate(const double size, const double distance, 
			       const double zenith)=0;

};

struct LZAEnergyEstimator : EnergyEstimator {    
    double getEstimate(const double size, const double distance, 
		       const double zenith=1.04719755);
};

struct SZAEnergyEstimator : EnergyEstimator {    
    double getEstimate(const double size, const double distance,
		       const double zenith=0.366519);
};

struct Z40EnergyEstimator : EnergyEstimator {    
    double getEstimate(const double size, const double distance,
		       const double zenith=0.698131);
};

struct Z50EnergyEstimator : EnergyEstimator {    
    double getEstimate(const double size, const double distance,
		       const double zenith=0.872664);
};


class EnergyEstimatorFactory {

 public: 
    static EnergyEstimatorFactory* instance();
    EnergyEstimator* getEstimator( std::string type );

 protected:
    EnergyEstimatorFactory();
    

 private:
    static EnergyEstimatorFactory *pinstance;

};


//double energyEstimator2(const double siz,const double dis );
//double LZAEnergyEstimator(const double size, const double dist );
//double SZAEnergyEstimator(const double size, const double dist );
#endif
