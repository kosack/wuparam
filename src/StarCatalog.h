#ifndef STARCATALOG_H
#define STARCATALOG_H

#include <string>
#include <list>
#include <vector>

/**
 * representation of a star
 */
class Star {

 public:
    std::string name;		//!< Name of star
    std::string spectclass;	//!< Spectral class of star
    double ra;			//!< precessed ra 
    double dec;			//!< precessed dec
    double magnitude;		//!< magnitude of the star
    int epoch;			//!< base epoch for the given ora/odec
    double ora, odec;		//!< base ra and dec 
    double tx, ty;		//!< tangential cartesian coordinates
    int catalog_id;		//!< Index of star catalog

    // operator definitions for sorting
    friend bool operator<(const Star &s1, const Star &s2 );
    friend bool operator>(const Star &s1, const Star &s2 );
    friend bool operator==(const Star &s1, const Star &s2 );
    friend bool operator!=(const Star &s1, const Star &s2 );

    friend std::ostream& operator<<( std::ostream &stream, Star &s );
    friend std::istream& operator>>( std::istream &stream, Star &s );

    
};




/**
 * A searchable star database
 */
class StarCatalog {


 public:
    StarCatalog():_cur_catalog_id(0)
	{_starcatalog.clear();} //!< Create a star catalog object

    void load( std::string stardbfilename );
    void findNearbyStars( double , double , double ,std::list<Star> &);
    void precessToDate( int utdate );   

 private:

    void tokenize(const std::string& str, std::vector<std::string>& tokens,
		  const std::string& delimiters = " ");
    
    void replaceChar( char search, char replace, std::string &);

    std::list<Star> _starcatalog;
    int _cur_catalog_id;

};



#endif
