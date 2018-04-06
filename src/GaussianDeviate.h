// GaussianDeviate.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef GAUSSIANDEVIATE_H
#define GAUSSIANDEVIATE_H

#include <string>
#include <gsl/gsl_rng.h>       /* Random num generators */
#include <gsl/gsl_randist.h>   /* Random distributions */
#include <gsl/gsl_sf_trig.h>

class GaussianDeviate {

 public:
    
    GaussianDeviate( std::string &dbdir );
    double getDeviate();
    double getFastDeviate();

    static const int DBSIZE = 1000000;

 private:

    void generate();
    
    gsl_rng *_rng;
    const gsl_rng_type *_rngtype;
    unsigned long _seed;

    std::string  _dbfilename;
    double  _deviates[DBSIZE];
    int _current;

};

inline double
GaussianDeviate::getDeviate() {

    return gsl_ran_ugaussian( _rng );

}


#endif
