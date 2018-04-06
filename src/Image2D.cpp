//
// Image2D.cpp
// Karl Kosack <kosack@hbar.wustl.edu>
//

#include <iostream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <assert.h>

#include "Image2D.h"
#include "Exceptions.h"
#include "ImageAnalyzer.h"
#include "Camera.h"
#include "ProgressBar.h"

using namespace std;


Image2D::Image2D( int nxbins, int nybins ) 
    : _nxbins(nxbins), _nybins(nybins) 
{

    // allocate memory for the grid and initialize to 0:

    _grid.resize( nxbins );
    for (int i=0; i<_grid.size(); i++) {
	_grid[i].resize( nybins );
    }

    _xcoord.resize( nxbins );
    _ycoord.resize( nybins );

    for (int i=0; i<_nxbins; i++) {
	for (int j=0; j<_nybins; j++) {
	    _grid[i][j] = 0;
	}
    }

    setCoordinateBox( 0,0,  nxbins,nybins );

}

/**
 * Set up the coordinate values of the image.  Specify the boundaries
 * of the image, not the positions of the center of the boundary
 * pixels!
 */
void
Image2D::setCoordinateBox( double minx, double miny, double maxx, double maxy ) {

    // _xcoord[] and _ycoord[] are the coordinates of the CENTER of the bin!

    double dx = (maxx-minx)/(double)(_nxbins);
    double dy = (maxy-miny)/(double)(_nybins);

    for (int i=0; i<_nxbins; i++) {
	_xcoord[i] = (double)i*dx + minx + dx/2.0;
    }

    for (int j=0; j<_nybins; j++) {
	_ycoord[j] = (double)j*dy + miny + dy/2.0;
    }

    _minx=minx;
    _maxx=maxx;
    _miny=miny;
    _maxy=maxy;

}

/**
 * Add another image to this image. Sums each pixel.  Images must be
 * of the same dimensions.
 */ 
void
Image2D::addImage( Image2D &im ) {

    if (im.getXDim() > _nxbins || im.getYDim() > _nybins) {
	throw MildAnalysisException("Image2D::addImage: image dimensions too large!");
    }
    
    for (int i=0; i<im.getXDim(); i++) {
	for (int j=0; j<im.getYDim(); j++) {

	    _grid[i][j] += im.getPixel(i,j);

	}
    }

}


/**
 * Add the specified value to the pixel which contains the coordinate
 * (x,y).  This is useful for histogramming in 2-D.
 */
void   
Image2D::
addHist( double x, double y, double val ) {

    int xbin,ybin;

    xbin = getXHistBinOf( x );
    ybin = getYHistBinOf( y );

    if (_nxbins>xbin && xbin>=0 && _nybins>ybin && ybin>=0) {
    
	_grid[xbin][ybin] += val;

    }
    else {
	_outside++;
    }



}


/**
 * Adds the value val to all pixels in the image that are within
 * radius of the point (x,y).  This provides another method to
 * radially smooth an image while filling it. It's slower than calling
 * applyRadialSmoothing() after an image has been created, but it's
 * more accurate.
 *
 * \todo: can speed this up considerably with geometric considerations
 * (by only iterating over range of pixels that fall withing radius*2
 * for example)
 */
void 
Image2D::
addHistRadially( double x, double y, double val, float radius) {
    
    register int i,j;
    float dist;

    // the following is a speed optimization: calculate the region
    // over which we should look (instead of iterating over the whole
    // grid)

    int istart = this->getXHistBinOf( x-2.0*radius );
    int iend = this->getXHistBinOf( x+2.0*radius );
    int jstart = this->getYHistBinOf( y-2.0*radius );
    int jend = this->getYHistBinOf( y+2.0*radius );
    
    if (istart<0) istart =0;
    if (jstart<0) jstart =0;
    if (iend>=_nxbins) iend =_nxbins-1;
    if (jend>=_nybins) jend =_nybins-1;

    // now look at each pixel in the region, and see if it's within
    // radius of the specified pixel:

    for (i=istart; i<=iend; i++) {
	for (j=jstart; j<=jend; j++) {
	    dist = hypot(x - _xcoord[i],y - _ycoord[j]);
	    if ( dist < radius ) {
		_grid[i][j]+=val;
	    }
	}
    }

}

/**
 * Write the image to disk.
 */
void
Image2D::
save( string filename ) {

    string fullfilename = filename + ".im2d";
    ofstream outfile(fullfilename.c_str());

    if (outfile.fail()) {
	throw MildAnalysisException("Couldn't write image: "+filename);
    }

//     outfile << _nxbins  << " "
// 	    << _nybins  << " "
// 	    << _minx <<" " 
// 	    << _miny <<" " 
// 	    << _maxy <<" " 
// 	    << _maxy <<" " << endl;
    

    for (int i=0; i<_nxbins; i++) {
	for (int j=0; j<_nybins; j++) {
	    
	    outfile << setw(15) << _xcoord[i] << " "
		    << setw(15) << _ycoord[j] << " "
		    << setw(15) << _grid[i][j] << " "
		    << endl;

	}
    }

    outfile.close();

}

/**
 * Write the image to disk.
 */
void
Image2D::
saveGrid( string filename ) {

    string fullfilename = filename + ".grid2d";
    ofstream outfile(fullfilename.c_str());

    if (outfile.fail()) {
	throw MildAnalysisException("Couldn't write image: "+filename);
    }

    for (int i=0; i<_nxbins; i++) {
	for (int j=0; j<_nybins; j++) {
	    outfile << setw(15) << _grid[i][j] << " ";
	}
	outfile << endl;
    }

    outfile.close();

}

/**
 * Load a saved image.
 */
void
Image2D::
load( string filename ) {

    ifstream infile( filename.c_str() );
    double nxbins, nybins, minx=1.0e9, miny=1.0e9, maxx=-1.0e9, maxy=-1.0e9;

    if (infile.fail()) {
	throw MildAnalysisException("Couldn't read image: "+filename);
    }


//     infile >> nxbins >> nybins >> minx >>  miny >> maxx >> maxy;

//     if ((nxbins != _nxbins) || (nybins != _nybins)) {
// 	throw CriticalAnalysisException(string("attempted to load image of ")+
// 					"incorrect size! "+ filename);
//     }

    for (int i=0; i<_nxbins; i++) {
	for (int j=0; j<_nybins; j++) {
	    infile >> _xcoord[i] >> _ycoord[j] >> _grid[i][j];
	    if (_xcoord[i] < minx) minx = _xcoord[i];
	    if (_xcoord[i] >= maxx) maxx = _xcoord[i];
	    if (_ycoord[i] < miny) miny = _ycoord[i];
	    if (_ycoord[i] >= maxy) maxy = _ycoord[i];
    	}
    }
    
    infile.close();


    double binsizex = _xcoord[1] - _xcoord[0];
    double binsizey = _ycoord[1] - _ycoord[0];

    // this requires edges of the box, whereas each pixel coordinate
    // is defined as the "center" of the box
    setCoordinateBox( minx-binsizex/2.0, miny-binsizey/2.0, 
		      maxx+binsizex/2.0, maxy+binsizex/2.0 );


}


/**
 * Set all image elements to 0
 */
void
Image2D::
clear() {

    for ( int i=0; i<_nxbins; i++) {
	for(int j=0; j<_nybins; j++) {
	    _grid[i][j] = 0.0;
	}
    }

}

/**
 * Write out a PGM: this is just a hack right now, should really make
 * a better 2D histogram class that does all this junk.
 */
void 
Image2D::
savePGM( string filename ) {

    int i,j;
    double min=1e10,max=-1e10;
    double scalefactor;
    unsigned char value;
    string fullfilename = filename + ".pgm";
    ofstream outfile(fullfilename.c_str());

    // determine min/max values for scaling...

    for ( i=0; i<_nxbins; i++) {
	for ( j=0; j<_nybins; j++) {
	    if (_grid[i][j] < min) min = _grid[i][j];
	    if (_grid[i][j] > max) max = _grid[i][j];
	}
    }

    scalefactor = 255.0/(max-min);
    
    // write PGM:

    outfile << "P5 "<<_nxbins<<" "<<_nybins<<" 255\n";

    for ( j=_nxbins-1; j>=0; j--) {
	for ( i=0; i<_nybins; i++) {

	    value = (unsigned char) (scalefactor*(_grid[i][j]-min));
	    outfile << value;

	}

    }

    outfile.close();
    
}



/**
 * Smooth the image to the specified radius.
 */
void 
Image2D::
applyRadialSmoothing( double r ) {
    
    int i,j,k,l;
    double xt,yt, xb,yb, dist2;
    double r2 = r*r;

    // Store the current image for use later...
    // probably need to define a copy constructor!

    Image2D oldimage = *(this);
    
    // Clear the current image

    for(i=0; i<_nxbins; i++)
	for(j=0; j<_nybins; j++)
	    _grid[i][j] = 0;

    // Generate a smoothed image based on the old stored image:

    for (i=0; i<_nxbins; i++) {

	for (j=0; j<_nybins; j++) {
	    
	    xt = _xcoord[i];
	    yt = _ycoord[j];

	    for (k=0; k<_nxbins; k++) {
                for (l=0; l<_nybins; l++) {
                    
                    xb = _xcoord[k];
                    yb = _ycoord[l];
		    
                    dist2  = pow(xb-xt,2) + pow(yb-yt,2);
		    
                    if (dist2 < r2) {
                        
			_grid[i][j] += oldimage.getPixel( k,l );
                        
                    }
                }
            }
	}
    }


}



double 
Image2D::
minValue() {
    
    double min = 1e10;

    for (int i=0; i< _nxbins; i++) {
	for (int j=0; j< _nybins; j++) {    
	    
	    if (_grid[i][j] < min)
		min = _grid[i][j];
	    
	}
    }

    return min;

}

double 
Image2D::
maxValue() {

    double max = -1e10;

    for (int i=0; i< _nxbins; i++) {
	for (int j=0; j< _nybins; j++) {    
	    
	    if (_grid[i][j] > max)
		max = _grid[i][j];

	}
    }

    return max;

}


/**
 * Expand the image by a factor of 2, and interpolate the results
 */
void
Image2D::
expand() {

    Image2D newimage( _nxbins*2, _nybins*2 );
    int j,k;
    double newpix;
    double a,b, x,y;
    double t1,t2,t3,t4;

    // Stuff the new image with pixels from the old one, using
    // trilinear filtering

    for (int n=1; n<newimage.getYDim(); n++) {
	for (int m=1; m<newimage.getXDim(); m++) {

	    // trilinear filter  

	    x = m/2.0;
	    y = n/2.0;

	    j = static_cast<int>(  x  );
	    k = static_cast<int>(  y  );

	    if ((j<1) || (j>(getXDim()-3))  || (k<1) || (k>(getYDim()-3)) ) {
		newpix = _grid[j][k];
	    }
	    else {
		
		a = x-j;
		b = y-k;
		
		t1 = a*(1-a)*(1-a)*_grid[j-1][k-1] + 
		    (1 - 2*a*a+a*a*a)*_grid[j][k-1] +
		    a*(1 + a-a*a)*_grid[j+1][k-1] - 
		    a*a*(1-a)*_grid[j+2][k-1];
		
		t2 = -a*(1-a)*(1-a)*_grid[j-1][k]+
		    (1 - 2*a*a + a*a*a)*_grid[j][k]+
		    a*(1+a-a*a)*_grid[j+1][k]-
		    a*a*(1-a)*_grid[j+2][k];
		
		t3 = -a*(1-a)*(1-a)*_grid[j-1][k+1]+
		    (1 - 2*a*a + a*a*a)*_grid[j][k+1]+
		    a*(1+a-a*a)*_grid[j+1][k+1]-
		    a*a*(1-a)*_grid[j+2][k+1];
		
		t4 = -a*(1-a)*(1-a)*_grid[j-1][k+2]+
		    (1 - 2*a*a + a*a*a)*_grid[j][k+2]+
		    a*(1+a-a*a)*_grid[j+1][k+2]-
		    a*a*(1-a)*_grid[j+2][k+2];
		
		
		newpix = -b*(1-b)*(1-b)*t1 + (1-2*b*b+b*b*b)*t2 +
		    b*(1+b-b*b)*t3 + b*b*(b-1)*t4;
		
	    }

	    newimage.setPixel( m,n, newpix );
	    
	}
    }
    

    // Now, copy the new image to the current one...
    
    _nxbins *= 2;
    _nybins *= 2;
    
    _grid.resize( _nxbins );
    for (int i=0; i<_grid.size(); i++) {
	_grid[i].resize( _nybins );
    }
    
    _xcoord.resize( _nxbins );
    _ycoord.resize( _nybins );
    
    for (int k=0; k<_nxbins; k++) {
	for (int l=0; l<_nybins; l++) {
	    _grid[k][l] = newimage.getPixel( k,l );
	}
    }    
    
    // resize the coordinate boxes as well..
    setCoordinateBox( _minx,_miny,_maxx,_maxy );
    
}



/**
 * Compute the cross-corrolation with the specified image.  Returns a
 * new Image2D of the correlation coefficients.  Remember to delete
 * the new image when you're done with it.
 *
 * \f[ C_{ij} = \sum_{k,l} I_1(k,l) * I_2((k+i) mod N,(l+j) mod M) \f]
 */
Image2D* 
Image2D::
crossCorrelateWith( Image2D &image ) {
  
    register int i,j,k,l,m,n;
    double sum=0;
    
    
    if (image.getXDim() != _nxbins || image.getYDim() != _nybins ) {
	cout << "Image2d::crossCorrolateWith(): image dimensions don't match"
	     << endl;
	cout <<"source: "<< _nxbins<<"x"<<_nybins
	     <<" dest: "<< image.getXDim()<<"x"<<image.getYDim()<<endl;

	throw AnalysisException( "Image dims don't match!" );
    }


    Image2D *corr = new Image2D( _nxbins, _nybins );
    corr->setCoordinateBox( _minx,_miny, _maxx,_maxy );
    
    for (i=0; i<_nxbins; i++) {
	for (j=0; j<_nybins; j++) {

	    // Calculate the C[i][j] correlation coefficient:
	    
	    sum=0;
	    for (k=0; k<_nxbins; k++) {
		for (l=0; l<_nybins; l++) {
		    sum += 
			(getPixel( k,l )) *
			(image.getPixel((k+i)%_nxbins,(l+j)%_nybins));
		}
	    }

	    m = (i+_nxbins/2)%_nxbins;
	    n = (j+_nybins/2)%_nybins;

	    corr->setPixel(m,n,
			   sum/((double)(_nxbins*_nybins)));
	}
    }

    return corr;

}


/**
 * Sets i, j to the indices of the max pixel
 */
void
Image2D::
getIndexOfMax( int &ii, int &jj ) {
    
    double max =-1e10;

    for (int i=0; i< _nxbins; i++) {
	for (int j=0; j< _nybins; j++) {    
	    
	    if (_grid[i][j] >= max) {
		max = _grid[i][j];
		ii=i;
		jj=j;
	    }

	}
    }

    assert(ii>=0 && ii<_nxbins);
    assert(jj>=0 && jj<_nybins);

    // now, (ii,jj) is the index of the max value
    
}
