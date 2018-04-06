// ProgressBar
// Karl Kosack <kosack@hbar.wustl.edu>

#include <iostream>
#include <iomanip>
#include <ctime>
#include <sys/time.h>

#include "ProgressBar.h"

using std::cout;
using std::cerr;
using std::endl;

#define ANSITERM 1

ProgressBar::ProgressBar( int total, std::string name ) 
    : _total(total),_last(0),_spin(0),_name(name),_avgpersec(0) {

    _spinner[0] = '|';
    _spinner[1] = '/';
    _spinner[2] = '-';
    _spinner[3] = '\\';

}

/**
 * Displays the progress bar, at position n
 */
void
ProgressBar::print( int n ) {

    int i;
    float percent;
    int secleft, minleft;
    double evtpersec;
    int prog;
    double delta = n - _last;
    const char ESC = 0x1b;

    if (n<=0) return;

    _spin++; if (_spin>=4) _spin=0;
    percent = (float)n/_total;
    prog = (int)(percent * 30); 
    cerr << _name<<": "<< _spinner[_spin] << " [" ;
    if (_total == 0) {
	cerr << "please wait...]                       ";
	cerr <<"  \r" << std::flush;
	return;
    }

    if (ANSITERM) {

	for (i=0; i<=prog; i++)   cerr<<ESC<<"[01;47m "<<ESC<<"[0m";//"#";
	for (i=prog+1; i<30; i++) cerr<<ESC<<"[01;44m "<<ESC<<"[0m"; //cerr << "_";
    }
    else {
	for (i=0; i<=prog; i++)   cerr<<"#";
	for (i=prog+1; i<30; i++) cerr<<"_";
    }

    _now = getTime();

    if (delta > 0) {
	evtpersec = (_now-_before > 0)? (delta/((_now-_before))):-1;
	_before = _now;
	_last = n;
    }


    _avgpersec = (evtpersec + _avgpersec)/2.0;
    secleft = (_avgpersec>0)?(int)((_total-n)/_avgpersec):0;
    minleft = secleft/60;
    secleft = secleft % 60;
    cerr << ESC << "[0m"
	 << "] " 
	 << std::setw(2) << (int)(percent*100) << "% ";


    if ( evtpersec > 1e8 || evtpersec < 0) 
	cerr << "??? per sec, ";
    else 
	cerr << std::setw(6) <<std::setprecision(4)<<evtpersec<<" per sec, ";
	
    if ( minleft < 60 && secleft > 0) {
	cerr << std::setw(2) << minleft <<":" 
	     << std::setw(2) << std::setfill('0')<< secleft << " left"
	     << std::setfill(' ');
    }
    else {
	cerr << "           ";
    }
	

    cerr <<"  \r" << std::flush;

    
}


/**
 * Returns the current system time in microseconds
 */
double
ProgressBar::getTime() {

    struct timeval tv;
    struct timezone tz;

    gettimeofday( &tv,&tz);

    return (tv.tv_sec + (double)tv.tv_usec/1.0e6);
    
}

/**
 * Print a string of spaces over the progress bar area.
 */
void
ProgressBar::printClear() {
    
    cerr << "                                              "
	"                                \r"
	 <<std::flush;
}
