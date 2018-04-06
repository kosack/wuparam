// PedestalFinder
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef _PEDESTALFINDER_H
#define _PEDESTALFINDER_H

#include <string>
#include <vector>
#include "Config.h"
#include "Types.h"

/**
 * Contains information about the pedestal for a single tube. 
 */
struct Pedestal {

    enum PEDESTAL_TYPES { GOOD, STAR, TUBEOFF };

    double pedestal;    //!< pedestal level in d.c.
    double dispersion;  //!< pedestal dispersion (sigma)
    int type;           //!< Code to mark special cases

};


/**
 *  Class to generate pedestal values for a series of runs.
 */
class PedestalFinder {

 public:
    
    PedestalFinder() : _badpedthresh(75.0) { setTubeOffThresholds(0.6,1.5); }
    ~PedestalFinder() {}
    void getOnSourcePeds( RunInfo &ri, std::vector<Pedestal> &ped, 
			  int telescope_id=0);
    void getOffSourcePeds( RunInfo &ri, std::vector<Pedestal> &ped, 
			   int telescope_id=0);
    void getPeds( RunInfo &ri, const std::string &id,
		  std::vector<Pedestal> &, int telescope_id=0 );
    void setTubeOffThresholds( double lower, double upper);
    void setBadPedThresh( double value ) {_badpedthresh=value;}
    
 private:

    double _lowthresh, _upthresh;
    double _badpedthresh;
    void readPeds( std::string &, std::vector<Pedestal> & );
    void writePeds( std::string &, std::vector<Pedestal> & );

};



#endif
 

