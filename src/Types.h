//  types.h - various typedefs for the analysis programs
//

#ifndef _TYPES_H
#define _TYPES_H

#include <valarray>

#define DEBUG 0

typedef std::valarray<double> Array_t;
typedef unsigned char Mask_t;


/**
 * A 2-d coordinate.
 */
struct Coordinate_t {
    double x;  //!< X coordinate
    double y; //!< Y coordinate
};




#endif
