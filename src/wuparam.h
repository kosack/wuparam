///////////////////////////////////////////////////////////////////////////
//  w  u  p  a  r  a  m 
//
//  Washington University data parameterization program for Whipple
//  Gamma-Ray telescope data.
//
//  by Karl Kosack (2000-08-21)
///////////////////////////////////////////////////////////////////////////

//======================================================================
// PROTOTYPES:
void  getCommandLineOptions(int , char**);
void showUsage();



//======================================================================
// Options and parameters:

struct Option_s{

    bool verbose;
    bool interactive;
    bool display;
    bool overwrite;
    bool clearcache;

} Options;






