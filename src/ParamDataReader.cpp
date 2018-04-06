// ParamDataReader.cpp
// Karl Kosack <kosack@hbar.wustl.edu>

#include <fstream>
#include <gsl/gsl_ntuple.h>
#include "DataWriter.h"
#include "ParamDataReader.h"
#include "Exceptions.h"
#include "Log.h"

using namespace std;

ParamDataReader::ParamDataReader( const string &filename )
    : DataReader(filename), _currow(0) {
    
    // check version information

    string versionfilename = filename+".ver";
    ifstream verfile(versionfilename.c_str());
    if (!verfile.fail()){
	string version;
	verfile >> version;
	if (version != OUTPUT_FILE_VERSION) {
	    cout << "NOTE: "<<filename<<" was parameterized using "
		 << "a different"<< endl
		 << "      version of WUPARAM ("<<version<<") "
		 << "and must be re-parameterized since the file format changed."<<endl;
	    throw CriticalAnalysisException("wrong file version for "
					    +filename);
	}
    }
    else {
	Logger *logger = Logger::instance();
	logger->printf("can't determine version of ntuple file, "
		       "assuming current");
    }

    // open the ntuple
    
    _ntuple = gsl_ntuple_open( (char*)getFilename().c_str(),
 			       &_ntuplerow, 
			       sizeof(_ntuplerow));

    if (_ntuple == NULL) {
	throw CriticalAnalysisException("Specified file didn't exist");
    }


    


}

ParamDataReader::
~ParamDataReader() {

    gsl_ntuple_close( _ntuple );

}

void
ParamDataReader::getNextEventRecord( HillasParameterization &param ) {

    if (gsl_ntuple_read( _ntuple )) {
	
	throw EOFException("end of ntuple file");
	
    }
    param = *((HillasParameterization *)_ntuple->ntuple_data);
    _currow++;

}

/**
 * Read in the Parameterized header record.  

 * \todo: should make a ParamHeaderRecord which is a subclass of
 * HeaderRecord which has all the extra parameterized info (like
 * whether zcuts was enabled, etc.)
 */
void 
ParamDataReader::getHeaderRecord( HeaderRecord &hdr ) {

    string hdrfile = getFilename() + ".hdr";
    ifstream infile( hdrfile.c_str() );
    int nadc;

    if (infile.fail())
	throw CriticalAnalysisException("Couldn't open the header file: "+hdrfile);

    infile >> _datasize 
	   >> hdr.ra 
	   >> hdr.dec
	   >> hdr.starttime
	   >> hdr.endtime
	   >> hdr.average_elevation
	   >> hdr.sourcename
	   >> hdr.num_telescopes;
    for (int i=0; i<hdr.num_telescopes; i++) {
	infile >> nadc;
	hdr.nadc.push_back(nadc);
    }
   
    infile.close();

}

bool
ParamDataReader::isDone() {

    if (_currow < _datasize-1) return false;
    else return true;

}
