#ifndef LOG_H
#define LOG_H

#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <ctime>
#include <string>
#include <cstdio>

/**
 * Error log singleton.  Calling the printf() function will write to
 * the log.   
 */
class Logger {

 public:

    static Logger* instance(); 
    void printf( char *format, ...);
    void info( char *format, ...);
    int getErrorCount();
    void printLog();
    std::string getFileName();

    ~Logger();

 protected:

    Logger();

 private:
   
    FILE *fp;
    int _errorcount;
    std::string _filename;
    static Logger *pinstance;
    fpos_t *_start_pos;

};



#endif
