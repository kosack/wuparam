// ParamDataReader.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef PARAMDATAREADER_H
#define PARAMDATAREADER_H


#include <string>
#include <gsl/gsl_ntuple.h>
#include "ImageAnalyzer.h"
#include "DataReader.h"

/**
 * Reads parameterized datafiles
 */
class ParamDataReader : public DataReader {

 public:

    ParamDataReader( const std::string & );
    ~ParamDataReader();

    void getHeaderRecord( HeaderRecord &hdr );
    void getNextEventRecord( HillasParameterization &p );
    int  size(){return _datasize;}
    std::string getTypeString(){ return std::string("Binary Ntuple"); }
    bool isDone();

 private:

    int _datasize;
    int _currow;
    gsl_ntuple *_ntuple;
    HillasParameterization _ntuplerow;

    
};


#endif
