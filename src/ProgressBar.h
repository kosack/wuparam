// ProgressBar
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <string>

/**
 * Displays a text-based progress bar with estimated time to finish
 * for a particular task.
 */ 
class ProgressBar {

 public:
    
    ProgressBar( int total, std::string name="Progress" );
    void setName( std::string name ){ _name=name; }
    void print( int );
    void printClear();
    double getTime();


 private:
    
    std::string _name;
    int _last;
    char _spinner[4];
    int _total;
    int _spin;
    double _now;
    double _before;
    double _avgpersec;

};



#endif
