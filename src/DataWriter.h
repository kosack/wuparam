// DataWriter.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef _DATAWRITER_H
#define _DATAWRITER_H

#include <string>
#include <fstream>
#include <gsl/gsl_ntuple.h>
#include "DataRecords.h"
#include "ImageAnalyzer.h"

#define OUTPUT_FILE_VERSION "WUPARAM-FILE-1.8.0"

/**
 * Base class for all data writers (parameterized data and cut
 * data). Specifies a standard interface.
 */
class DataWriter {

 public:
    DataWriter(const std::string &filename,bool=false) :_filename(filename){;}
    virtual ~DataWriter(){;}
	
    virtual std::string getTypeString() = 0;
    std::string getFilename() { return _filename; } 

 private:
    std::string _filename;


};



/**
 * Base class for writing HillasParameterizations to disk. Subclasses
 * of this should implement the various file formats that can be
 * written.
 */
class ParamDataWriter : public DataWriter {

 public:

    ParamDataWriter(const std::string &filename, bool ow=false);
    virtual ~ParamDataWriter() {;}

    void writeHeader( HeaderRecord &hdr );
    virtual void writeParameterization(HillasParameterization &param) = 0;
    virtual void writeParameterization(SimShowerRecord &p1,
				       HillasParameterization &p2) {
	writeParameterization( p2 );
    }
    virtual int  size() =0;
    virtual std::string getTypeString() =0;

};


/**
 * Outputs parameterized events to a text file.
 */
class ParamDataWriterText : public ParamDataWriter {

 public:
    ParamDataWriterText( const std::string &, bool ow=false );
    ~ParamDataWriterText();

    void writeParameterization( HillasParameterization & );
    void writeParameterization(SimShowerRecord &,
			       HillasParameterization &);

    int  size() { return 0;}

    std::string getTypeString() { return std::string("tab delimited text"); }

 private:

    std::ofstream _outfile;

};

/**
 * Outputs parameterized events to an N-Tuple (for use in for example,
 * PAW).
 */
class ParamDataWriterNtuple : public ParamDataWriter {

 public:
    ParamDataWriterNtuple( const std::string &, bool ow=false );
    ~ParamDataWriterNtuple();

    void writeParameterization( HillasParameterization & );
    int  size() { return _numrows; }
    std::string getTypeString() { return std::string("binary N-tuple"); }

 private:

    gsl_ntuple *_ntuple;
    HillasParameterization _ntuplerow;
    int _numrows;
    
};

std::ostream &operator<<( std::ostream &stream,const SimShowerRecord &p);
//std::istream &operator>>( std::istream &stream, SimShowerRecord &p);

#endif
