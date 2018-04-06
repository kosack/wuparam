
#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include "Exceptions.h"
#include "Log.h"
using namespace std;


Logger* Logger::pinstance = NULL;

Logger* 
Logger::
instance() {
    if (pinstance == NULL)
	pinstance = new Logger();
    
    return pinstance;
}


/**
 * writes to the log, but doesn't increment the warning counter
 */
void 
Logger::info( char *format, ...) {


    va_list ptr;
    
    va_start( ptr, format );
    vfprintf( fp, format, ptr );
    va_end( ptr );
    fprintf( fp, "\n");

    fflush(fp);


}

/**
 * Writes warning string to the log and to stdout
 */
void 
Logger::printf( char *format, ...) {

    const char ESC = 0x1b;

    va_list ptr;
    
    fprintf( stdout, "%c[01;31m *** %c[0m",ESC,ESC); // escape codes for bright red

    va_start( ptr, format );
    vfprintf( fp, format, ptr );
    va_end( ptr );

    va_start( ptr, format );
    vfprintf( stdout, format, ptr );
    va_end( ptr );

    fprintf( fp, "\n");
    fprintf( stdout, "\n");
    
    fflush(fp);

    _errorcount++;

}


int
Logger::
getErrorCount() {
    return _errorcount;
}


void
Logger::
printLog() {

    std::ifstream logfile( "warning-log.txt");
    char buffer[128];

    if (logfile.fail()) {
	cout << "COULDN'T OPEN LOGFILE FOR READING!"<< endl;
	return;
    }

    cout << "===================================" << endl;
    cout << "THE FOLLOWING WARNINGS WERE LOGGED:" << endl;
    cout << "===================================" << endl;
    while (logfile) {
	logfile.getline( buffer, 128);
	cout << buffer << endl;
    }
    logfile.close();

}

string 
Logger::
getFileName() {
    return _filename;
}

Logger::Logger()  {

    char *thetime;
    time_t lt;

    lt = time( NULL );
    thetime = ctime( &lt );
    thetime[strlen(thetime)-1] = '\0';

    _filename = "warning-log.txt";
	
    fp = fopen( _filename.c_str(), "a");

    if (fp == NULL) {
	throw(AnalysisException("Couldn't open log file: "+_filename));
    }

    _errorcount=0;
    
    fprintf( fp,"#----------------- LOG STARTED: %s -----------------#\n",
	     thetime ); 

}


Logger::~Logger() {
    
    fclose( fp );
    
} 
