#include <iostream>
#include <fstream>
#include <list>
#include <vector>
#include <string>
#include <gsl/gsl_math.h>
#include <iomanip>

#include "StarCatalog.h"
#include "Exceptions.h"
#include "Image2D.h"

using namespace std;


bool operator<(const Star &s1, const Star &s2 ) {
    return s1.name < s2.name;
}
bool operator>(const Star &s1, const Star &s2 ) {
    return s1.name > s2.name;
}
bool operator==(const Star &s1, const Star &s2 ) {
    return s1.name == s2.name;
}
bool operator!=(const Star &s1, const Star &s2 ) {
    return s1.name != s2.name;
}


/**
 * Load a database of stars from the specified .edb file (such as
 * those included with XEphem.
 */
void
StarCatalog::
load( string stardbfilename ) {

    ifstream stardb( stardbfilename.c_str() );
    string line;
    vector<string> tokens;
    vector<string> dhms;
    Star tmpstar;
    double d,h, m,s;

    tmpstar.catalog_id = _cur_catalog_id;

    if (stardb.fail()){
	throw MildAnalysisException("Couldn't load star catalog");
    }

    while ( getline( stardb, line, '\n') ) {

	if (line[0]=='#') continue;
	tokens.clear();
	tokenize( line, tokens, "," );

	if (tokens.size() == 6) {
	    
	    tmpstar.name = tokens[0];
	    replaceChar( ' ','_', tmpstar.name );
	    tmpstar.spectclass = tokens[1];
	    
	    // read and convert ra from HMS -> radians

	    dhms.clear();
	    tokenize( tokens[2], dhms, ":" );
	    if (dhms.size() == 3) {
		h = atof(dhms[0].c_str());
		m = atof(dhms[1].c_str());
		s = atof(dhms[2].c_str());
		tmpstar.ra = (M_PI/12.0) * ( h + m/60.0 + s/3600.0 );
		tmpstar.ora = tmpstar.ra;
	    }
	    else {
		cerr << "Skipping bad line in star catalog..."<< endl;
	    }

	    // read and convert dec from DMS -> radians

	    dhms.clear();
	    tokenize( tokens[3], dhms, ":" );
	    if (dhms.size()==3){
		d = atof(dhms[0].c_str());
		m = atof(dhms[1].c_str());
		s = atof(dhms[2].c_str());
	    }
	    else if (dhms.size() == 2) {
		d = atof(dhms[0].c_str());
		m = atof(dhms[1].c_str());
		s = 0;
	    }
	    else {
		cerr << "Skipping bad line in star catalog..."<< endl;
		continue;
	    }
	    tmpstar.dec = (M_PI/180.0) * ( d + m/60.0 + s/3600.0 );
	    tmpstar.odec = tmpstar.dec;

	    // read in magnitude and epoch

	    tmpstar.magnitude = atof(tokens[4].c_str()); 
	    tmpstar.epoch     = atoi(tokens[5].c_str()); 
	    tmpstar.tx   = 0; // tangential coords unknown at this point
	    tmpstar.ty   = 0; // unknown at this point
	    
	    _starcatalog.push_back(tmpstar);

	}

    } 
    
    _cur_catalog_id++;


}

/**
 * Produces a vector of stars that fall within a specified radius of a
 * specified ra+dec coordinate.  Also figures out the position of each
 * star in tangential cartesian coordinates, relative to the center
 * position for easy plotting.
 *
 * \param ra right ascension where the camera is pointing
 * \param dec declination where the camera is pointing
 * \param r radius (in degrees) around which to look for stars (should be 
 *     FOV of camera)
 * \param nearbystars empty vector of Stars which will be filled up
 * with the results
 *
 * \todo: check that this ra and dec from the tracking computer are
 * already precessed, or from j2000 or something.  If not already
 * precessed, they need to be.
 */
void 
StarCatalog::
findNearbyStars( double ra, double dec, double r,list<Star> &nearbystars) {

    list<Star>::iterator star;        
    double rrad = r*M_PI/180.0;
    double denom;

    for (star=_starcatalog.begin(); star != _starcatalog.end(); star++) {

	// is the star within radius*2 of the declination
	if ( fabs( star->dec - dec ) < rrad*2 ) {

	    // calculate the tangential coordinates of the star
	    denom = (sin(dec)*sin(star->dec) + 
		     cos(dec)*cos(star->dec)*cos(star->ra - ra));

	    if (denom>0) {

		// Not sure why minus sign is needed for tx, but
		// seems to be wrong otherwise

		star->tx = -(180.0/M_PI)*(cos(star->dec)*
					 sin(star->ra-ra))/denom;
		
		star->ty = (180.0/M_PI)*(cos(dec)*sin(star->dec)-
					 sin(dec)*cos(star->dec)*
					 cos(star->ra-ra))/denom;
		
		
		// if it's in range:

		if (sqrt(pow(star->tx,2)+pow(star->ty,2)) < r){

		    nearbystars.push_back( *star );

		}
	    }
	}
    }

}


/**
 * Precess star catalog to the current epoch.
 *
 * \param utdate for example '020506' for May 6, 2002
 *
 * \todo: Should check this - it looks like stars are moving too much
 * over the course of a year or so, so it's possible this is not
 * working exactly right.  However, it's pretty close.
 */
void
StarCatalog:: 
precessToDate( int utdate ) {

    // expansion parameters:
    const double ara  = 1.2812323;
    const double bra  = 0.0003879;
    const double cra  = 0.0000101;
    const double adec = 0.5567530;
    const double bdec = -0.0001185;
    const double cdec = -0.0000116;

    double yy = (int)(utdate/10000.0);
    double nn = (int)((utdate-yy*10000)/100.0);
    double epoch = 2000.0 + yy + nn /12.0;

    double alpham, deltam, alpha,delta;
    double d2r = 180.0/M_PI;
    double tp=0,np=0,mp=0;
    double lastepoch=0;

    list<Star>::iterator star;        

    // precess each star

    for (star=_starcatalog.begin(); star!=_starcatalog.end(); star++) {

	if (lastepoch != star->epoch) {
	    tp = (epoch - star->epoch)/100.0;
	    mp = ara*tp  + bra*tp*tp  + cra*tp*tp*tp;
	    np = adec*tp + bdec*tp*tp + cdec*tp*tp*tp;
	}

	// precess to mid-time point

	alpham = star->ora*d2r  + 0.5*(mp+np*sin(star->ora)*tan(star->odec));
	deltam = star->odec*d2r + 0.5*np*cos(star->ora);

	// precess to epoch

	alpha = star->ora*d2r + mp+np*sin(alpham)*tan(deltam);
	delta = star->odec*d2r + np*cos(alpham);
	
	star->ra  = alpha/d2r;
	star->dec = delta/d2r;

	lastepoch = star->epoch;

	
    }
    
}

/**
 * Function to split a string into tokens
 */ 
void 
StarCatalog::
tokenize(const string& str, vector<string>& tokens,
	 const string& delimiters)
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}

void
StarCatalog::
replaceChar( char search, char replace, string &str) {

    for (int i=0; i<str.length(); i++) {
	if (str[i] == search)
	    str[i] = replace;
    }

}


ostream& operator<<( ostream &stream, Star &star ) {

    stream << setw(10) << star.ora  << ' '
	   << setw(10) << star.odec  << ' '
	   << setw(10) << star.epoch  << ' '
	   << setw(10) << star.tx   << ' '
	   << setw(10) << star.ty   << ' '
	   << setw(10) << star.ra   << ' '
	   << setw(10) << star.dec  << ' ' 
	   << setw(10) << star.magnitude << ' '
	   << setw(10) << star.catalog_id << ' '
	   << setw(10) << star.name << ' '
	   << setw(10) << star.spectclass << ' '
	   << endl;
    
    return stream;

}

istream& operator>>( istream &stream, Star &star ) {

    stream >> star.ora
	   >> star.odec
	   >> star.epoch
	   >> star.tx 
	   >> star.ty 
	   >> star.ra 
	   >> star.dec 
	   >> star.magnitude
	   >> star.catalog_id 
	   >> star.name 
	   >> star.spectclass;

    return stream;
}
