// DataReader.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef _DATAREADER_H
#define _DATAREADER_H

#include <iostream>
#include <fstream>
#include <hdf.h>
#include <string>
#include "Types.h"
#include "Exceptions.h"
#include "DataRecords.h"
#include "Config.h"

class RawDataReaderFactory;

/**
 * Base class for all data readers (raw, parameterized, cut)
 */
class DataReader {

 public:
    DataReader( const std::string &filename ):_filename(filename){;}
    virtual ~DataReader(){;}
    
    std::string getFilename() { return _filename; }

    virtual bool isDone()=0;
    virtual int size() =0;
    virtual std::string getTypeString() =0;

 private:
    std::string _filename;

};



/**
 * Reads raw telescope data and produces coresponding records...
 */
class RawDataReader : public DataReader {

 public:

    RawDataReader( const std::string &filename );
    
    virtual void getHeaderRecord( RawHeaderRecord &header )=0;
    virtual void getNextEventRecord( RawEventRecord &event )=0;
    virtual int  size()=0;
    virtual std::string getTypeString() =0;

    virtual int getType()=0;
    enum RunType { OTHER,HDF4,SIM,VTEXT };

    virtual bool isDone()=0;
    virtual ~RawDataReader(){;}

    enum EventTypes { EVENT, PEDESTAL, UNKNOWN };

};


/**
 * Implementation of RawDataReader for HDF4 input files.
 */ 
class RawDataReaderHDF4 : public RawDataReader {

 public:

    RawDataReaderHDF4( const std::string &filename );
    ~RawDataReaderHDF4();

    void getHeaderRecord( RawHeaderRecord &header );
    void getNextEventRecord( RawEventRecord &event );
    int  size();
    bool isDone() { return _is_file_done; }
    std::string getTypeString();
    int getType(){ return RawDataReader::HDF4; }

 private:

    uint8* _infobuf;
    uint8* _eventbuf;
    int    _infobufsize;
    int    _eventbufsize;
    int    _hfileid;
    int32  _eventref;
    int32  _headerid;
    int32  _headerref;
    int32  _eventid;
    short* _eventadc;
    int    _nadc;
    bool   _is_header_initialized;
    bool   _is_event_initialized;
    bool   _is_file_done;

};

/**
 * Implementation of RawDataReader for simulation .tubes files
 */ 
class RawDataReaderSim : public RawDataReader {

 public:

    RawDataReaderSim( const std::string &filename );
    ~RawDataReaderSim();

    void getHeaderRecord( RawHeaderRecord &header );
    void getNextEventRecord( RawEventRecord &event );
    int  size();
    bool isDone() { return _is_file_done; }
    std::string getTypeString() { 
	return std::string("Simulation output file");
    }
    int getType(){ return RawDataReader::SIM; }

    SimShowerRecord getSimShowerRecord() { return _curshower; };

 private:


    SimShowerRecord _curshower;
    bool _is_file_done;
    std::ifstream _infile;
    int _count;
    int _nadc;
    double _faketime;

};

/**
 * Implementation of RawDataReader VERITAS text data as outputed by
 * Paul Rebillot's charge integration routines.
 */ 
class RawDataReaderVText : public RawDataReader {

 public:

    RawDataReaderVText( const std::string &filename );
    ~RawDataReaderVText();

    void getHeaderRecord( RawHeaderRecord &header );
    void getNextEventRecord( RawEventRecord &event );
    int  size();
    bool isDone() { return _is_file_done; }
    std::string getTypeString() { 
	return std::string("VERITAS Text Data");
    }
    int getType(){ return RawDataReader::VTEXT; }

 private:

    bool _is_file_done;
    std::ifstream _infile;
    int _count;
    int _nadc;
    int _nevents;

};



// Inserters
std::ostream &operator<<( std::ostream &stream, RawEventRecord &record );
std::ostream &operator<<( std::ostream &stream, RawHeaderRecord &record );


/**
 * Responsible for locating files and creating the correct data reader
 * object for the detected file type.  Implemented as a Singleton.
 */
class RawDataReaderFactory {

 public:
    
    static RawDataReaderFactory* instance();
    RawDataReader* getReader( const RunInfo &ri, string runid );
    std::string locateFile( const RunInfo &ri, 
			    const string basename,
			    const string ext );

    void clearCache( std::string );


 protected:

    RawDataReaderFactory();
    

 private:

    static RawDataReaderFactory *pinstance;


};


#endif


