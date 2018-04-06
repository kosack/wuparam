#ifndef HIRES_H
#define HIRES_H

#include <valarray>
#include <vector>

void filter_init( int npix, 
		  std::valarray<double> &xc, std::valarray<double> &yc,
		  std::valarray<double> &newxc,
		  std::valarray<double> &newyc, 
		  int msub=3, double afiltr=0.25, int nc=2, double spas=0 );

double filter_image(  std::valarray<double> &im, 
		    std::valarray<double> &fim  );


#endif
