#ifndef PLOTMAKER_H
#define PLOTMAKER_H

#include <plot.h>
#include "Image2D.h"
#include "Exceptions.h"
#include "Camera.h"
#include "Types.h"
#include <stack>

struct Pedestal;
struct HillasParameterization;
struct MuonParameterization;
class Star;

/**
 * Plots stuff to any PlotUtils graphics context
 */
class PlotMaker {

 public:

    PlotMaker( std::string type, std::string filename="", 
	       double minx=-3.1,double miny=-3.1, 
	       double maxx=3.1, double maxy=3.1);
    ~PlotMaker() {
	pl_closepl_r( plotter );
	pl_deletepl_r( plotter );
	fclose(_fp);
    }

    enum ColorMap {
	COLORMAP_BLUES,
	COLORMAP_REDBLUE,
	COLORMAP_YELLOWGREEN,
	COLORMAP_INVERSEGREY,
	COLORMAP_GREY,
	COLORMAP_RAINBOW,
	COLORMAP_HOT,
	COLORMAP_PRINTABLE,
    };

    enum ContourLineStyle {
	LINE_SOLID,
	LINE_DOTTED,
	LINE_DASHED
    };

    enum CameraPlotType {
	CAMERA_IMAGE,
	CAMERA_PEDS,
	CAMERA_PEDDISPS,
	CAMERA_TUBENUMS
    };

    void setAxisBox( double minx, double miny, double maxx, double maxy ) {
	_minx=minx; 
	_maxx=maxx; 
	_miny=miny; 
	_maxy=maxy;
    }

    void setAxisBox( Image2D &im2d ) {
	setAxisBox( im2d.getMinX(),
			im2d.getMinY(),
			im2d.getMaxX(),
			im2d.getMaxY() );
    }

    void setAxisBox( Camera &cam ) {
	setAxisBox( cam.getMinX(),
		    cam.getMinY(),
		    cam.getMaxX(),
		    cam.getMaxY() );
    }

    void setSpace( double minx, double miny, double maxx, double maxy ){
	_pageminx = minx;
	_pageminy = miny;
	_pagemaxx = maxx;
	_pagemaxy = maxy;
	pl_fspace_r( plotter, _pageminx, _pageminy, _pagemaxx, _pagemaxy); 
    }

    void flush() { pl_flushpl_r( plotter ); }
    void newPage(){ pl_erase_r( plotter ); }
    void pushState();
    void popState();

    void rotate( double theta );
    void plotAxes();
    void plotTitle( std::string title );
    void plotSubtitle( std::string title );
    void plotXTicAt( double x );
    void plotYTicAt( double y );
    void plotXMinorTicAt( double x );
    void plotYMinorTicAt( double y );
    void contourPlot( Image2D &im2d, double zmin, double zmax, double dz );
    void plotRADecGrid( Image2D &ragrid, Image2D &decgrid );

    // Higher-level plot commands:
    void plot( Image2D &im2d );
    void plot( Camera &cam, const Array_t &image, const vector<Pedestal> &peds,
	       vector<int> &cleanpixels, CameraPlotType plottype=CAMERA_IMAGE);
    void plot( HillasParameterization &p );
    void plot( MuonParameterization &p );
    void plot( Star &s, int symbol=16 );
    void plot( vector<double> &xpoints, vector<double> &ypoints );
    void plot( string text, double x, double y  );
    void plot( TelescopeArray &array );
    
    // Lower-level plot commands
    void plotEllipse( double x, double y, double r1, double r2,double theta=0);
    void plotMarker( double x, double y );

    void setColorMap( ColorMap c ) {_state.top().colormap = c; }
    void setColorMap( string name ) {
	if      (name == "blues")_state.top().colormap = COLORMAP_BLUES;
	else if (name == "redblue")_state.top().colormap = COLORMAP_REDBLUE;
	else if (name == "yellowgreen")_state.top().colormap = COLORMAP_YELLOWGREEN;
	else if (name == "inversegrey")_state.top().colormap = COLORMAP_INVERSEGREY;
	else if (name == "grey")_state.top().colormap = COLORMAP_GREY;
	else if (name == "rainbow")_state.top().colormap = COLORMAP_RAINBOW;
	else if (name == "hot")_state.top().colormap = COLORMAP_HOT;
	else if (name == "printable")_state.top().colormap = COLORMAP_PRINTABLE;
	else throw MildAnalysisException("unknown colormap "+name);
    }

    void getMappedColor(double value, int &red, int &green, int &blue);

    void setPixelScale( double pixscale ) {_pixelscale = pixscale;}
    void setScaleFactor( double factor ) {_scalefactor = factor; }
    void setTicLevel( int l ) { _ticlevel = l; }

    void setColorScale( double lower, double upper) {
	if (lower<=upper) {
	    _state.top().colorscale_lower=lower;
	    _state.top().colorscale_upper=upper;
	}
    }
    void enableAutoScale( bool val=true) { _state.top().autoscale=val;  }
    void setContourColorName( string name ) {
	_state.top().contourcolor = name;
    }
    void setTicMarkColorName( string name ) {
	_state.top().ticcolor = name;
    }
    void setContourLineStyle( ContourLineStyle i ){
	_state.top().contourlinestyle = i;
    }
    void setContourLineStyle( string style ){
	if (style == "solid") setContourLineStyle(LINE_SOLID);
	else if (style == "dotted") setContourLineStyle(LINE_DOTTED);
       	else if (style == "dashed") setContourLineStyle(LINE_DASHED);
	else throw MildAnalysisException("unknown line style "+style);
    }
    
    void setContourLineWidth( double width ) {_state.top().contourwidth=width;}
    void setLabelColorName( string name ) {_state.top().labelcolor = name; }
    void setLabelAngle( double theta ){_state.top().labelangle=theta;}
    void setLabelSize( double size ){_state.top().labelsize=size; }
    void setMarkerType( int type ) {_state.top().markertype = type; }
    void setMarkerSize( double size ) {_state.top().markersize = size; }
    void setMarkerColorName( string name ) {_state.top().markercolor = name; }
    


 private:

    void setLineMod();
    void niceRaDecContours(double ra,  double dec, 
		  double xcont[],double ycont[], 
		  double ra_nice[], double dec_nice[], 
		  int drawx[], int drawy[], 
		  char ra_str[17][128],
			   char dec_str[17][128]);

    static const int NCONTSTEPS =17;
    
    plPlotter *plotter;
    plPlotterParams *plotter_params;
    FILE *_fp;

    double _pageminx,_pageminy, _pagemaxx, _pagemaxy;
    double _minx,_miny,_maxx,_maxy;
    double  _pixelscale;
    double  _maxmagnitude;
    int     _numkeycolors;
    double  _keyoffset;
    double  _keywidth;
    double  _labeloffset;
    double  _titleoffset;
    double  _subtitleoffset;
    int     _ticlevel;
    double  _ticmarksize;
    double  _scalefactor;

    struct PlotState {

	int     contourlinestyle;
	int     markertype;
	double  markersize;
	string  markercolor;
	bool    use_inward_tics;
	ColorMap colormap;
	string  contourcolor;
	double  contourwidth;
	string  ticcolor;
	string  labelcolor;
	double  labelangle;
	double  labelsize;
	bool    autoscale;
	double  colorscale_lower;
	double  colorscale_upper;

    };

    std::stack<PlotState> _state;

    double xIntercept(double , double , double , double ,double , double );
    double yIntercept(double , double , double , double ,double , double );

};


struct ContourSegment {
    Coordinate_t start;
    Coordinate_t end;
};

#endif
