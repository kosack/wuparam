//
// DataWriter.cpp
//
// DataWriter.cpp
// Karl Kosack <kosack@hbar.wustl.edu>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <gsl/gsl_ntuple.h>
#include <unistd.h>

#include "ImageAnalyzer.h"
#include "DataWriter.h"
#include "Exceptions.h"
#include "Types.h"

using namespace std;

ParamDataWriter::ParamDataWriter(const string &filename,bool overwrite)
    : DataWriter(filename,overwrite) {

    if (overwrite == false) {
	// test if header file exists, if it does, exit.
	string headerfilename = filename+".hdr";
	string versionfilename = filename+".ver";
	ifstream test(headerfilename.c_str());
	if (!test.fail()) {
	        
	    // Also check that the version of the file is current
	    // otherwise, must reparameterize
	    ifstream verfile(versionfilename.c_str());
	    if (!verfile.fail()){
		
		string version;
		verfile >> version;
		if (version != OUTPUT_FILE_VERSION) {
		    cout << "NOTE: "<<filename<<" was parameterized using "
			 << "a different"<< endl
			 << "      version of WUPARAM ("<<version<<") "
			 << "and must be re-parameterized"<<endl;
		}
 		else {
	 	    throw MildAnalysisException(filename+
						" is already parameterized.");
	        }
	    }
	}
	test.close();
    }
    else {
	cout << "Overwriting " << filename << " with new version"<< endl;
    }

}


//===============================================================================

ParamDataWriterText::ParamDataWriterText( const string &filename, 
					  bool ow ) 
    : ParamDataWriter(filename,ow) {
 
    _outfile.open( getFilename().c_str() );

    if (_outfile.fail()) 
	throw CriticalAnalysisException("Couldn't open output file");
    
}

ParamDataWriterText::~ParamDataWriterText() {
    
    _outfile.close();
    
}

void
ParamDataWriterText::
writeParameterization( HillasParameterization &param ) {

    _outfile << param << endl;
        
}

void
ParamDataWriterText::writeParameterization(SimShowerRecord &p1,
                                           HillasParameterization &p2){
    _outfile << p1 << p2 << endl;
}


//===============================================================================

ParamDataWriterNtuple::
ParamDataWriterNtuple( const string &filename,
		       bool overwrite) 
    : ParamDataWriter(filename,overwrite), _ntuple(NULL), _numrows(0) {
    
    _ntuple = gsl_ntuple_create( (char *)filename.c_str(), 
				 &_ntuplerow, 
				 sizeof(_ntuplerow) );
    
    if (_ntuple == NULL) {
	throw CriticalAnalysisException("Couldn't create N-tuple");
    }

    if (overwrite) {
	// remove the old header file if it exists...
	string headerfilename = filename+".hdr";
	unlink( headerfilename.c_str() );
    }
}

ParamDataWriterNtuple::~ParamDataWriterNtuple() {

    gsl_ntuple_close( _ntuple );

}


void
ParamDataWriterNtuple::
writeParameterization(  HillasParameterization &param ) {
    
    _ntuplerow = param;

    gsl_ntuple_write( _ntuple );

    _numrows++;

}


/**
 * Write out header.  Should be done after all NTuple rows have been
 * written, so the number of rows is correct.
 */ 
void 
ParamDataWriter::
writeHeader( HeaderRecord &hdr ) {

    string hdrname = getFilename() + ".hdr";
    ofstream outfile( hdrname.c_str() );
    
    string tsourcename = hdr.sourcename;

    // replace spaces in source name with underscore, so it's all one
    // token
    for (int i=0; i<tsourcename.length(); i++)
	if (tsourcename[i] == ' ')
	    tsourcename[i] = '_';
	

    outfile.precision(20);
    outfile << size() << " "
	    << hdr.ra << " "
	    << hdr.dec << " "
	    << hdr.starttime << " "
	    << hdr.endtime << " "
	    << hdr.average_elevation << " "
	    << tsourcename << " "
	    << hdr.num_telescopes << " ";

    for (int i=0; i<hdr.num_telescopes; i++) {
	outfile << hdr.nadc[i] << " ";
    }
    outfile << endl;

    outfile.close();

    // also write out version file:
    string vername = getFilename() + ".ver";
    ofstream voutfile( vername.c_str() );
    voutfile << OUTPUT_FILE_VERSION << endl;
    voutfile.close();

}

/**
 * Inserter for outputting SimShowerRecord in text format
 */
ostream& operator<<(ostream &stream,const SimShowerRecord &p){
 
    stream << (int) p.event_number << " "
           << (int) p.primary_type << " "
           << setprecision(14) << p.primary_energy << " "
           << setprecision(14) << p.impact_parameter.x << " "
           << setprecision(14) << p.impact_parameter.y << " "
           << setprecision(14) << p.direction_cos.x << " "
           << setprecision(14) << p.direction_cos.y << " " ;
 
    return stream;
}



/**
 * Extractor for SimShowerRecords
 */
// istream &operator>>( istream &stream, SimShowerRecord &p) {

//     stream >> p.event_number
//            >> p.primary_type
//            >> p.primary_energy 
//            >> p.impact_parameter.x
//            >> p.impact_parameter.y
//            >> p.direction_cos.x 
//            >> p.direction_cos.y;
    
//     return stream;
    
// }
