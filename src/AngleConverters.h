#ifndef ANGLECONVERTERS_H
#define  ANGLECONVERTERS_H

inline void 
radToDMS( double rad, double &d, double &m, double &s ) {

    double decimal = rad *180.0/M_PI;
    double tmp;

    d   = (int)trunc(decimal);

    tmp = decimal - d;
    m = (int)trunc(tmp*60);

    tmp = (tmp*60) - m;
    s = tmp*60;
    
}


inline void 
radToHMS( double rad, double &h, double &m, double &s ) {

    double decimal = rad *12.0/M_PI;
    double tmp;

    h   = (int)trunc(decimal);

    tmp = (decimal - h)*60;
    m = (int)trunc(tmp);

    tmp = (tmp - m)*60;
    s = tmp;
    
}


inline std::string 
getHMSString( double rad ) {

    double h,m,s;
    char tmp[50];

    radToHMS( rad, h, m, s );

    sprintf( tmp, "%dh %dm %ds (%f rad) (%f deg)",
	     (int)h, (int)m, (int)s, rad, rad*12.0/M_PI);

    return std::string(tmp);

}

inline std::string 
getDMSString( double rad ) {

    double d,m,s;
    char tmp[50];

    radToDMS( rad, d, m, s );

    sprintf( tmp, "%dd %d' %d\" (%f rad) (%f deg)",
	     (int)d, (int)m, (int)s, rad, rad*180.0/M_PI);

    return std::string(tmp);

}

#endif
