// GainFinder.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef _GAINFINDER_H
#define _GAINFINDER_H

#include <string>
#include "Config.h"
#include "Types.h"

/**
 * Retrieve or calculate nitrogen gains, given a nitrogen run ID. 
 */
class GainFinder {

 public:

    GainFinder() {setThresholds(50.0,1000.0,0.50);}
    ~GainFinder() {}
    void getGains( RunInfo &ri, Array_t &gains, int telescope_id=0 );
    void setThresholds( double,double,double );

 private:    

    void readGains( std::string &, Array_t & );
    void writeGains( std::string &, Array_t & );

    double _minadc, _maxadc, _minpercent;

};


#endif
