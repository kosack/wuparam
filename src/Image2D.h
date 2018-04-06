
#ifndef IMAGE2D_H
#define IMAGE2D_H

#include "Types.h"
#include <cstdio>
#include <list>
#include <vector>
#include <string>

/**
 * Representation of a 2-D raster image (for 2D analysis)
 */
class Image2D {

 public:


    Image2D( int nxbins=39, int nybins=39 );


    void   setCoordinateBox(double minx,double miny, double maxx,double maxy);
    
    double getXCoord( int i ) { return _xcoord[i];}
    double getYCoord( int j ) { return _ycoord[j];}
    double getMinX(){ return _minx; }
    double getMinY(){ return _miny; }
    double getMaxX(){ return _maxx; }
    double getMaxY(){ return _maxy; }

    void   load( std::string filename );
    void   save( std::string filename );
    void   savePGM(std::string filename);
    void   saveGrid(std::string filename);
    void   clear();

    void   setPixel( int i, int j, double val ) {_grid[i][j] = val; }
    double getPixel( int i, int j ) {return _grid[i][j]; }
    void   addToPixel( int i, int j, double val) {_grid[i][j] += val;}
    void   addHist( double x, double y, double val );
    void   addHistRadially( double x, double y, double val, float radius);
    int    getOutsideHits(){ return _outside; }

    void   addImage( Image2D &im );
    void   applyRadialSmoothing( double r );
    
    int    getXDim() { return _nxbins; }
    int    getYDim() { return _nybins; }

    double minValue();
    double maxValue();
    void   getIndexOfMax( int &i, int &j );

    void   expand(); 

    Image2D* crossCorrelateWith( Image2D & );

    int getXHistBinOf( double x ) {
	return (int)((x - _minx)*(double)_nxbins / (_maxx - _minx));  
    }
    int getYHistBinOf( double y ) {
	return (int)((y - _miny)*(double)_nybins / (_maxy - _miny));  
    }

 private:

    int _nxbins;
    int _nybins;
    double _minx, _miny, _maxx, _maxy;
    int _outside;   		
    std::vector< std::vector<double> > _grid;// 2-d array of the image
    std::vector<double> _xcoord; // array of x-coordinates of bins
    std::vector<double> _ycoord; // array of y-coordinates of bins

};






#endif
