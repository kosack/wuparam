// DataReader.cpp
// Karl Kosack <kosack@hbar.wustl.edu>

#include <iostream>
#include <string>
#include <iomanip>
#include <glob.h>
#include <sstream>
#include <sys/types.h>
#include <unistd.h>
#include <hdf.h>

#include "ProgressBar.h"
#include "DataReader.h"
#include "Exceptions.h"
#include "Config.h"

using namespace std;


RawDataReader::
RawDataReader( const string &filename ) : DataReader(filename) {
    

}

RawDataReaderHDF4::
RawDataReaderHDF4( const string &filename )  : RawDataReader( filename ) {

    int status;

    // =================================================
    // Allocate memory and set defaults... 
    // =================================================

    _infobufsize = 
	2*sizeof(int32) +
	3*sizeof(float64) +
	80*sizeof(char8);

    _infobuf = new uint8[ _infobufsize ];
    _eventbufsize = -1;
    _is_header_initialized = false;
    _is_event_initialized = false;
    _is_file_done = false;

    // =================================================
    // open the HDF file
    // =================================================    

    _hfileid = Hopen( getFilename().c_str(), DFACC_READ, 0 );

    if (_hfileid == -1) {
	throw CriticalAnalysisException( "HDF Open Failed" );
    }

    // ==================================================
    // Attach to the required VDatas
    // ==================================================

    status = Vstart( _hfileid );

    _eventref = VSfind( _hfileid, "10M Event" );
    if( _eventref < 0 ) {
	throw CriticalAnalysisException("Couldn't locate 10M Event Vdata in "
					+getFilename());
    }
    
    _eventid = VSattach( _hfileid, _eventref, "r");
    if( _eventid < 0 ) {
	throw CriticalAnalysisException("Couldn't find event in "
					+getFilename());
    }
        

    _headerref = VSfind( _hfileid, "10M Run Information" );
    _headerid  = VSattach( _hfileid, _headerref, "r");
    if( _headerid < 0 ) {
	throw CriticalAnalysisException("Couldn't find run information in "
					+getFilename());
    }
        

    
}

RawDataReaderHDF4::~RawDataReaderHDF4() {

    VSdetach( _headerid );
    VSdetach( _eventid );
    Hclose( _hfileid );

    delete [] _infobuf;
    delete [] _eventbuf;
    delete [] _eventadc;

}

void
RawDataReaderHDF4::
getHeaderRecord( RawHeaderRecord &header ) {

    static const int NFIELDS = 6;
    static const char FIELDNAMES[128] = 
	"VERSION,NADC,UTC_START,SOURCE_NAME,SOURCE_RA,SOURCE_DEC";
    int ret;
    VOIDP field[NFIELDS];
    char sourcename[128];
    int version;

    header.num_telescopes = 1; // whipple data always has 1 telescope
    header.nadc.resize(1);
    int nadc;
        
    // =================================================
    // Read the header information
    // =================================================

    ret = VSsetfields( _headerid, FIELDNAMES );
    if (ret == -1) {
	throw CriticalAnalysisException("VSsetfields() failed");
    }

    ret = VSread( _headerid, _infobuf, 1, FULL_INTERLACE );
    if(ret != 1){
	throw CriticalAnalysisException("Couldn't read the run header record");
    }

    // =================================================
    // Unpack the fields into the correct variables
    // =================================================

    field[0] = &version;
    field[1] = &(nadc);
    field[2] = &(header.starttime);
    field[3] = &(sourcename[0]);
    field[4] = &(header.ra);
    field[5] = &(header.dec);

    VSfpack( _headerid, 
	     _HDF_VSUNPACK, 
	     FIELDNAMES, 
	     _infobuf,
	     _infobufsize,
	     1, NULL,
	     (VOIDP *)field);

    sourcename[80] = '\0';
    header.sourcename = sourcename;
    header.nadc[0] = nadc;
    header.windowsize = 0;

   
    // =================================================
    // Allocate memory for the event record if needed...
    // =================================================
    if (_is_header_initialized == false) {
	try {
	    _eventbufsize = (nadc+3)*sizeof(int16)+3*sizeof(float64);
	    _eventbuf = new uint8[_eventbufsize];
	    _eventadc = new short[ nadc ];
	    _nadc = nadc;
	}
	catch(bad_alloc e) {
	    throw CriticalAnalysisException( "Memory allocation failed: event buffer");
	}

	_is_header_initialized = true;

    }

}

/**
 * Return number of records in data file
 */
int
RawDataReaderHDF4::
size() {

    return VSelts( _eventid );

}


void
RawDataReaderHDF4::
getNextEventRecord( RawEventRecord &event ) {

    static const char FIELDNAMES[128] =
	"CODE,GPSTIME,OSCTIME,LIVETIME,ADC,ELEVATION,AZIMUTH";
    int ret=0;
    short code;
    VOIDP field[7];

    if (_is_file_done == true ) {
	throw EOFException( "no more events" );
    }

    // =================================================
    // If not initialized, fetch the header to determine
    // buffer sizes (and allocate memory).  Then 
    // set up to read from the event vdata...
    // =================================================
    
    if ( _is_header_initialized == false  ) {
	
	RawHeaderRecord header;
	getHeaderRecord( header );
	
    }
    
    if ( _is_event_initialized == false ) {


	ret = VSsetfields( _eventid, FIELDNAMES );
	if (ret == -1) {
	    throw CriticalAnalysisException("VSsetfields() failed");
	}

	_is_event_initialized = true;
	
    }


    // =================================================
    // Read an event record from the file...
    // =================================================

    ret = VSread( _eventid, _eventbuf, 1, FULL_INTERLACE );
    if(ret != 1){
	_is_file_done = true;
    }

    // =================================================
    // Unpack the data and put it into the record object
    // =================================================
    
    field[0] = &(code);
    field[1] = &(event.gpstime);
    field[2] = &(event.osctime);
    field[3] = &(event.livetime);
    field[4] = &(_eventadc[0]);
    field[5] = &(event.elevation);
    field[6] = &(event.azimuth);

    VSfpack( _eventid, 
	     _HDF_VSUNPACK, 
	     "CODE,GPSTIME,OSCTIME,LIVETIME,ADC,ELEVATION,AZIMUTH", 
	     _eventbuf,
	     _eventbufsize,
	     1, NULL,
	     (VOIDP *)field);

    // convert the adc counts to doubles: eventually gains will be
    // applied, so we want floating point
    if ((int)event.adc.size() != _nadc)  event.adc.resize( _nadc );
    for( int i=0; i<_nadc; i++ ) {
	event.adc[i] =  (double) _eventadc[i];
    }


    // Check this! 
    if (code & 0x4) {
	event.type = RawDataReader::EVENT;
    }
    else if (code & 0x1) {
	event.type = RawDataReader::PEDESTAL;
    }
    else {
	event.type = RawDataReader::UNKNOWN;
    }

    event.telescope_id = 0;


}


string
RawDataReaderHDF4::
getTypeString() { 
    
    uint32 maj,min,rel;
    char str[100];

    Hgetfileversion( _hfileid, &maj,&min,&rel,str );

    string s((const char *)str);
    return s;
}

/**
 * Open .rec input file
 * \bug: if file doesn't exist, doesn't always fail!
 */ 
RawDataReaderSim::RawDataReaderSim( const string &filename ) 
    : RawDataReader( filename ), _is_file_done(false), 
      _count(0), _nadc(490), _faketime(0)
{

    _infile.open( filename.c_str() );

    if (!_infile || _infile.fail()) {
	throw CriticalAnalysisException("Open failed for "+filename);
    }
    
}

RawDataReaderSim::~RawDataReaderSim() {
    
    _infile.close();
    
}


/**
 * Return number of records in the file
 */ 
int
RawDataReaderSim::size() {
    
    ifstream infile( getFilename().c_str() );
    string line;
    int count=0;

    while ( !infile.eof() && getline( infile, line, '\n') ){
	if (line[0] == 'Z' || line[0] == 'Q')
	    count++;
    }

    infile.close();

    return count;

}

/**
 * Get header record. 
 * \todo: fill in the correct values
 */
void
RawDataReaderSim::
getHeaderRecord( RawHeaderRecord &header ) {

    header.nadc.resize(1);
    header.nadc[0] = 490;
    header.sourcename = "Simulation Data";
    header.starttime = 0;
    header.ra = 0.0;
    header.dec = 0.0;
    

}

/**
 * Returns the next RawEventRecord from the simulation. 
 */ 
void
RawDataReaderSim::
getNextEventRecord( RawEventRecord &event ) {

    string line, temp;

    if (event.adc.size() != _nadc);
	event.adc.resize(_nadc);

    while (!_infile.eof() && getline( _infile, line, '\n')) {

	if ( line[0] == 'S' ) {
	    // Process a shower record Only shower records which
	    // precede a 'Q' record are ones that actually triggered
	    // the telescope, so others can be ignored.
	    
	    istringstream ist( line );
	    ist >> temp  // the 'S'
		>> _curshower.event_number
		>> _curshower.primary_type
		>> _curshower.primary_energy
		>> _curshower.impact_parameter.x
		>> _curshower.impact_parameter.y
		>> _curshower.direction_cos.x
		>> _curshower.direction_cos.y;

	}

	if ( line[0] == 'Z' || line[0] =='Q' ) {
	    // Process a pedestal (Z) or actual event (Q)
	    
	    if (line[0] == 'Z') event.type = RawDataReader::PEDESTAL;
	    else event.type = RawDataReader::EVENT;
	    
	    istringstream ist( line );
	    ist >> temp // record type
		>> temp // placeholder
		>> event.telescope_id; // telescope number 
	    
	    for (int i=0; i<_nadc; i++) {
		if (!(ist >> event.adc[i])) 
		    throw CriticalAnalysisException("Not enough PMTS in data");
	    }

	    event.gpstime = _faketime;
	    event.osctime = _faketime;
	    event.livetime =_faketime;
	    _faketime += 1e-6;
	    event.elevation = 0;
	    event.azimuth = 0;

	    // found event, so stop looking
	    return;

	}

    }

    // if we got here, then the file is done
    _is_file_done = true;
    return;

}




//========================================================================


/**
 * Return number of records in the file
 */ 
int
RawDataReaderVText::size() {


    RawHeaderRecord header;
    getHeaderRecord( header );

    if (_nevents != -1) {

	return _nevents;
    }
    
    cout << "RawDataReaderVText: size of file not found in header, calculating it"<<endl
	 << "                    with brute force method! (please wait)"<< endl;

    ProgressBar progress( 0, "CalcSize" );

    ifstream infile( getFilename().c_str() );
    string line;
    int count=0;
    int pcount=0;

    while ( !infile.eof() && getline( infile, line, '\n') ){
	if (line[0] == 'v' || line[0] == 'p')
	    count++;
	if (pcount++ >= 2048) {pcount=0; progress.print(1);}
	
    }

    infile.close();
    progress.printClear();

    _nevents = count;

    cout << "RawDataReaderVText: size = "<< _nevents << endl;

    return count;

}


/**
 * Open .vtext input file
 */ 
RawDataReaderVText::RawDataReaderVText( const string &filename ) 
    : RawDataReader( filename ), _is_file_done(false), _count(0), _nadc(500),_nevents(-1)
{
    _infile.open( filename.c_str() );

    if ( !_infile || _infile.fail()) {
	throw CriticalAnalysisException("Open failed for "+filename);
    }

}

RawDataReaderVText::~RawDataReaderVText() {
    _infile.close();
}

/**
 * Get header record. 
 * \todo: fill in the correct values
 */
void
RawDataReaderVText::
getHeaderRecord( RawHeaderRecord &header ) {
    
    string line;
    vector<string> tokens;
    bool foundheader=false;
    int count=0;

    ifstream infile( getFilename().c_str() );

    if (!infile) {
      	throw CriticalAnalysisException("Open failed for "+getFilename());
    }


    //cout << "DEBUG: RawDataReaderVText: searching for header record..."
    //	 << endl;

    header.num_telescopes = 1;
    header.nadc.resize(header.num_telescopes);

    
    while (!infile.eof() && getline( infile, line, '\n')) {

	count++;

	if  (count > 100) break; // bail out early - header shuld be near top!

	tokenize( line, tokens, "\t " );

	if (tokens[0] == "h") {

	    if (tokens.size() <= 5) {
		infile.close();
		throw CriticalAnalysisException("Bad header "+getFilename());
	    }

	    // found the header
	    // TODO: make this work for multi-telescope data
	    header.nadc[0] = atoi(tokens[1].c_str());
	    _nadc = header.nadc[0]; // TODO: make this the maximum of all nadc
	    header.starttime = atof(tokens[2].c_str());	    
	    header.sourcename = tokens[3];
	    header.ra = atof(tokens[4].c_str());	    
	    header.dec = atof(tokens[5].c_str());	    

	    if (tokens.size() >= 7) {
		_nevents = atoi(tokens[6].c_str());
	    }
	    foundheader=true;
	    header.windowsize = atoi(tokens[7].c_str());
	    break;

	}

    }

    infile.close();

    if (foundheader==false) {
	throw CriticalAnalysisException("Couldn't find the header 'h'"
					"record in "
					+getFilename());
    }

}

/**
 * Returns the next RawEventRecord from the simulation. 
 *
 * format is: type event_num gpstime osctime livetime elevation
 * azimuth telescope_id adcvalues[499]
 */ 
void
RawDataReaderVText::
getNextEventRecord( RawEventRecord &event ) {

    string line, temp;
    int nsamp;

    if (_is_file_done == true) 
	throw MildAnalysisException("There are no more events!");


    while (!_infile.eof() && getline( _infile, line, '\n')) {


	if (line.size() == 0 ) continue;

	// --------------------------------------
	// Process veritas event
	// --------------------------------------
	if ( line[0] == 'v' || line[0] == 'p') {

	    _count++;
	    
	    if (line[0] == 'v') event.type=RawDataReader::EVENT;
	    else if (line[0] == 'p') event.type=RawDataReader::PEDESTAL;

	    istringstream ist(line);


	    ist >> temp; // type    
	    ist >> temp; // event number
	    ist >> event.gpstime
		>> event.osctime
		>> event.livetime 
		>> temp //event.elevation is a short
		>> temp //event.azimuth  is a short
		>> event.telescope_id;
	    ist >> nsamp; // num samples;
		
	    if (event.adc.size() != _nadc);
	    event.adc.resize(_nadc);

	    for (int i=0; i< _nadc; i++) {
		ist >> event.adc[i]; 
	    }
	    	 
	    // found event, so stop looking
	    return;
   
	}
	else {
	    cout << "Skipping unknown event at line " <<_count
		 << " of "<< getFilename()<<endl;
	}

    }

    // if we got here, then the file is done
    _is_file_done = true;
    return;

}



//========================================================================

/**
 * Inserter for printing RawHeaderRecords.
 */
ostream &operator<<( ostream &stream, RawHeaderRecord &record ) {

    stream <<"NTEL: "<< setw(5) << record.num_telescopes << " ";
    for (int i=0; i<record.num_telescopes; i++) {
	cout   <<"N"<<i<<":   "<< setw(5) << record.nadc[i] << " ";
    }
    cout    <<"RA:  "<< setw(12) << record.ra << " "
	    <<"DEC: "<< setw(12) << record.dec << " "
	    <<"STT: "<< setw(12) << record.starttime << " "
	    <<"SRC: "<< record.sourcename << " ";

    return stream;

}

/**
 * Inserter for printing RawEventRecords.
 */
ostream &operator<<( ostream &stream, RawEventRecord &record ) {

    stream <<"typ: " << setw(5)  << record.type << " "
	   <<"gps: " << setw(12) << record.gpstime << " "
	   <<"osc: " << setw(12) << record.osctime << " "
	   <<"liv: " << setw(12) << record.livetime<< " "
	   <<"ele: " << setw(6)  << record.elevation<< " "
	   <<"azi: " << setw(6)  << record.azimuth << " "
	   <<"tel: " << setw(6)  << (int)record.telescope_id<<" ";
       
    return stream;

}


RawDataReaderFactory* RawDataReaderFactory::pinstance = NULL;

/**
 * Returns pointer to the global RawDataReaderFactory instance
 */
RawDataReaderFactory* 
RawDataReaderFactory::
instance() {

    if (pinstance == NULL)
	pinstance = new RawDataReaderFactory();
    
    return pinstance;

}

RawDataReaderFactory::
RawDataReaderFactory() {


}


/**
 * Returns a pointer to a data reader for the given id
 */
RawDataReader*
RawDataReaderFactory::
getReader( const RunInfo &ri, string runid ) {

    string filename;
    int type=RawDataReader::OTHER;
    RawDataReader *rawdatareader;

    // first look for an HDF4 File:

    filename = locateFile( ri, runid, "hdf" );

    if (filename == "") {

	// Try to find a vtext file:
	filename = locateFile( ri, runid, "vtext");
	if (filename == "") {
	
	    // Try to find a simulation file
	    filename = locateFile( ri, runid, "rec" );
	    if (filename =="") {
		
		// bail out!
		throw CriticalAnalysisException("No file with id "+runid
						+ " could be found under "
						+ ri.datadir );
	    }
	    else {
		type = RawDataReader::SIM;
	    }
	}
	else {
	    type = RawDataReader::VTEXT;
	}
    }
    else {
	type = RawDataReader::HDF4;
    }

    
    switch (type) {

    case RawDataReader::HDF4:
	rawdatareader = new RawDataReaderHDF4(filename);
	break;
	
    case RawDataReader::SIM:
	rawdatareader = new RawDataReaderSim(filename);
	break;

    case RawDataReader::VTEXT:
	rawdatareader = new RawDataReaderVText(filename);
	break;

    default:
	throw(CriticalAnalysisException("Unknown run data type for "+runid));
	break;
    }
    
    return rawdatareader;

}


/**
 * Locate a file in a directory structure, uncompressing if needed.
 *
 * \param startdir directory under which to look
 * \param basename base name of the file (without extensions)
 * \param file extension (i.e. hdf)
 *
 * \returns full pathname to the located file.
 */
string 
RawDataReaderFactory::
locateFile( const RunInfo &ri, const string basename, const string ext ) {
    
    glob_t globbuf;
    string pattern;

    // Check in cache directory first:
    pattern = ri.cachedir + basename +"."+ext;
    glob(pattern.c_str(), 0, NULL, &globbuf);

    if (globbuf.gl_pathc > 0)  
	return string(globbuf.gl_pathv[0]);

    // Next, look in DataDir for the original, uncompressed file...

    for ( int i=0; i<5; i++ ) {
	pattern = ri.datadir;
	for (int j=0; j<i; j++) {
	    pattern += "*/";
	}
	pattern += basename + "." + ext;
	glob(pattern.c_str(), GLOB_APPEND, NULL, &globbuf);

	if (globbuf.gl_pathc > 0) {
	    return string(globbuf.gl_pathv[0]);
	}
    }

    // If we've got this far, then no file has been located. Check for
    // a compressed version...

    for ( int i=0; i<5; i++ ) {
	pattern = ri.datadir;
	for (int j=0; j<i; j++) {
	    pattern += "*/";
	}
	pattern += basename + "." + ext+".bz2";
	if (i==0) 
	    glob(pattern.c_str(), 0, NULL, &globbuf);
	else
	    glob(pattern.c_str(), GLOB_APPEND, NULL, &globbuf);

	if (globbuf.gl_pathc > 0) {

	    // Found file, needs to be uncompressed...

	    string found = globbuf.gl_pathv[0];
	    string uncompressed = ri.cachedir+basename+"."+ext;
	    string command = "bunzip2 -v -c "+found+" > "+uncompressed+"~";
	    string command2= "mv "+uncompressed+"~ "+uncompressed;
	    cout << "Uncompressing " << basename
		 <<", this may take some time... "<< endl;
	    system(command.c_str());
	    system(command2.c_str());

	    // Record the fact that this file has been cached...
	    
	    pid_t pid = getpid();
	    char pidstr[16];
	    sprintf( pidstr, "%d", pid );
	    string listname = ri.cachedir+string(pidstr)+"-cached_files.list";
	    ofstream filelist(listname.c_str(), ios::app);
	    filelist << basename+"."+ext << endl;
	    filelist.close();

	    return uncompressed;

	}

    }

    // File wasn't found anywhere!
    return string("");

}

/**
 * Purge the cache of stored files.
 */
void
RawDataReaderFactory::
clearCache(string cachedir) {
    
    pid_t pid = getpid();
    char pidstr[16];
    sprintf( pidstr, "%d", pid );
    string listname = cachedir+string(pidstr)+"-cached_files.list";
    ifstream filelist(listname.c_str());
    string file;

    if (!filelist.fail()) {

	cout << "Clearing cached uncompressed files in "
	     << cachedir <<"..."<<endl;

	while (filelist >> file) {
	    string fullfile = cachedir+file;
	    unlink( fullfile.c_str() );
	}

	filelist.close();

    }

    unlink(listname.c_str());

}




