// Exceptions.h
// Karl Kosack <kosack@hbar.wustl.edu>

#ifndef _EXCEPTIONS_H
#define _EXCEPTIONS_H

//
//  Define a set of exceptions for use by the analysis programs
//

#include <string>
#include <iostream>
#include <exception>
#include <stdexcept>

/**
 * Standard exception thrown by analysis routines
 * \todo: AnalysisException should be a subclass of the STL exception
 */

class AnalysisException : public std::runtime_error {
 public:
    AnalysisException( std::string str ): std::runtime_error(str) {;}
};

/**
 * Thrown on errors in the analysis where the program cannot recover
 */
class CriticalAnalysisException : public AnalysisException {
    
 public:
    CriticalAnalysisException(const std::string &str)
	: AnalysisException(str) {;}

};


/**
 * Thrown for warnings and recoverable errors
 */
class MildAnalysisException : public AnalysisException {
    
 public:
    MildAnalysisException(const std::string &str)
	: AnalysisException(str) {;}

};

/**
 * Thrown by datareaders when the file is over
 */ 
class EOFException : public MildAnalysisException {

 public:
    EOFException(const std::string &str)
	: MildAnalysisException(str) {;}
    
};


/**
 * Thrown by datareaders when the file is over
 */ 
class RangeException : public MildAnalysisException {

 public:
    RangeException(const std::string &str, int code=0)
	: MildAnalysisException(str) {;}
    
};

#endif
