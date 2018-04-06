/* -*- c-basic-offset: 4 -*- */
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>
#include <list>
#include <cmath>
#include <gsl/gsl_math.h>
#include <plot.h>

#include  "PlotMaker.h"
#include "ImageAnalyzer.h"
#include "Camera.h"
#include "PedestalFinder.h"
#include "MuonImageAnalyzer.h"
#include "StarCatalog.h"
#include "Types.h"
#include "AngleConverters.h"


using namespace std;

PlotMaker::PlotMaker( string type, string filename,
		      double minx,double miny, double maxx, double maxy ) {
    
    PlotState s;
    
    _minx=_maxx=_miny=_maxy=0;
    _pageminx = minx;
    _pagemaxx = maxx;
    _pageminy = miny; 
    _pagemaxy = maxy;
    _keyoffset = 0.2;
    _keywidth = 0.1;
    _numkeycolors = 20;
    _pixelscale = -1;
    _maxmagnitude = 10.0;
    _ticlevel = 0;
    _labeloffset=0.1;
    _titleoffset=0.12;
    _subtitleoffset=0.03;
    _ticmarksize=0.06;
    _scalefactor = 1.0; 

    s.use_inward_tics = true;
    s.contourcolor = "LightGray";
    s.contourlinestyle = LINE_SOLID;
    s.ticcolor = "black";
    s.labelcolor="white";
    s.contourwidth = 0.009;
    s.markertype = 1;
    s.markersize = 0.2;
    s.markercolor = "white";
    s.labelangle = 0.0;
    s.labelsize = 0.1;
    s.autoscale = true;

    _state.push(s);
    
    setColorMap( COLORMAP_BLUES );

    plotter_params = pl_newplparams();
    pl_setplparam( plotter_params, "PAGESIZE", (char*)"letter");
    pl_setplparam( plotter_params, "BITMAPSIZE", (char*)"640x640");
    pl_setplparam( plotter_params, "USE_DOUBLE_BUFFERING", (void*)"yes");    

    if (type=="X") 
	_fp = stdout;
    else
	_fp = fopen( filename.c_str(), "w" );
    
    if ((plotter = pl_newpl_r( type.c_str(), stdin, _fp, stderr,
			       plotter_params)) == NULL) {
	throw CriticalAnalysisException("Couldn't create plotter");
    }
    
    pl_openpl_r( plotter );
    pl_fspace_r( plotter, _pageminx, _pageminy, _pagemaxx, _pagemaxy); 
    pl_flinewidth_r( plotter, 0.25); // line thickness in user coordinates
    pl_pencolorname_r(plotter, "black"); // line color defaults
    pl_erase_r( plotter );	// erase Plotter's graphics display
    pl_fmove_r(plotter ,0.0, 0.0); // position the graphics cursor
    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, _scalefactor*0.15);
    

}


/**
 * Plot a set of axes around the AxisBox region.  
 *
 * \todo: don't plot minor tics where a major tic is
 */ 
void 
PlotMaker::
plotAxes() {

    double niceminx, niceminy, nicemaxx, nicemaxy;
    int niceticlevel;
    double p;
    double x,y;

    pl_savestate_r( plotter );

    pl_pentype_r( plotter, 1);
    pl_filltype_r( plotter, 0);
    pl_flinewidth_r( plotter, 0.02*_scalefactor );

    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, _scalefactor*0.12);
    
    // plot some nice major tic marks

    //    niceticlevel = (int) log10(_maxx-_minx);
    //    _ticlevel = 

    p = pow( 10.0, -_ticlevel );
    niceminx =  ceil(_minx * p)/p;
    niceminy =  ceil(_miny * p)/p;
    nicemaxy = floor(_maxy * p)/p;
    nicemaxx = floor(_maxx * p)/p; 

    for ( x = niceminx; x<=nicemaxx; x+=1.0/p ) {
	plotXTicAt( x );
    }

    for ( y = niceminy; y<=nicemaxy; y+=1.0/p ) {
	plotYTicAt( y );
    }

    // plot the minor tic marks (at one higher ticlevel magnitude)

    p = pow( 10.0, -(_ticlevel-1) );
    niceminx =  ceil(_minx * p)/p;
    niceminy =  ceil(_miny * p)/p;
    nicemaxy = floor(_maxy * p)/p;
    nicemaxx = floor(_maxx * p)/p; 
   
    for ( x = niceminx; x<=nicemaxx; x+=1.0/p ) {
	plotXMinorTicAt( x );
    }

    for ( y = niceminy; y<=nicemaxy; y+=1.0/p ) {
	plotYMinorTicAt( y );
    }


    // plot the bounding box
    pl_pentype_r( plotter, 1);
    pl_filltype_r( plotter, 0);
    pl_flinewidth_r( plotter, 0.02 );
    pl_fbox_r( plotter, _minx, _miny, _maxx, _maxy );



    pl_restorestate_r( plotter );

}


void 
PlotMaker::
plotXTicAt( double x ) {

    char text[64];
    double ticoffs;

    pl_flinewidth_r( plotter, 0.02*_scalefactor );

    if ( _state.top().use_inward_tics == true )
	ticoffs = -_ticmarksize*_scalefactor;
    else 
	ticoffs =  _ticmarksize*_scalefactor;
        
    pl_pencolorname_r( plotter,_state.top().ticcolor.c_str() );
    pl_fline_r( plotter, x ,_maxy,  x,_maxy + ticoffs ) ;
    pl_fline_r( plotter, x ,_miny,  x,_miny - ticoffs ) ;
    pl_pencolorname_r( plotter, "black" );

    pl_fmove_r( plotter, x, _miny - _labeloffset );
    sprintf( text, "%.1f", x );
    pl_alabel_r( plotter, 'c','t',  text );

}

void 
PlotMaker::
plotYTicAt( double y ){

    char text[64];
    double ticoffs;

    pl_flinewidth_r( plotter, 0.02*_scalefactor );

    if (_state.top().use_inward_tics == true )
	ticoffs = -_ticmarksize*_scalefactor;
    else 
	ticoffs =  _ticmarksize*_scalefactor;

    pl_pencolorname_r( plotter,_state.top().ticcolor.c_str() );
    pl_fline_r( plotter, _maxx , y,  _maxx + ticoffs, y ) ;
    pl_fline_r( plotter, _minx , y,  _minx - ticoffs, y ) ;
    pl_pencolorname_r( plotter, "black" );

    pl_fmove_r( plotter, _minx - _labeloffset, y );
    sprintf( text, "%.1f", y );
    pl_alabel_r( plotter,'r','c',  text );
   

}

void 
PlotMaker::
plotXMinorTicAt( double x ) {

    double ticoffs;

    pl_flinewidth_r( plotter, 0.01 );

    if ( _state.top().use_inward_tics == true )
	ticoffs = - _ticmarksize/1.5 *_scalefactor;
    else 
	ticoffs = _ticmarksize/1.5 *_scalefactor;
        
    pl_pencolorname_r( plotter,_state.top().ticcolor.c_str() );
    pl_fline_r( plotter, x ,_maxy,  x,_maxy + ticoffs ) ;
    pl_fline_r( plotter, x ,_miny,  x,_miny - ticoffs ) ;
    pl_pencolorname_r( plotter, "black" );
}

void 
PlotMaker::
plotYMinorTicAt( double y ){

    double ticoffs;

    pl_flinewidth_r( plotter, 0.01 );

    if (_state.top().use_inward_tics == true )
	ticoffs = - _ticmarksize/1.5 *_scalefactor;
    else 
	ticoffs = _ticmarksize/1.5 *_scalefactor;

    pl_pencolorname_r( plotter,_state.top().ticcolor.c_str() );
    pl_fline_r( plotter, _maxx , y,  _maxx + ticoffs, y ) ;
    pl_fline_r( plotter, _minx , y,  _minx - ticoffs, y ) ;
    pl_pencolorname_r( plotter, "black" );

}


void
PlotMaker:: 
plotTitle( std::string title ) {

    pl_savestate_r( plotter );

    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, _scalefactor*0.15);    

    pl_fmove_r( plotter, (_minx+_maxx)/2.0, _maxy+_titleoffset );
    pl_alabel_r( plotter, 'c','b', title.c_str() );

    pl_restorestate_r(  plotter );

}

void
PlotMaker:: 
plotSubtitle( std::string title ) {

    pl_savestate_r( plotter );

    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, _scalefactor*0.09);    

    pl_fmove_r( plotter, (_minx+_maxx)/2.0, _maxy+_subtitleoffset );
    pl_alabel_r( plotter, 'c','b', title.c_str() );

    pl_restorestate_r(  plotter );

}

void 
PlotMaker::
plot( Image2D &im2d ) {

    double minx = im2d.getMinX();
    double miny = im2d.getMinY();
    double maxx = im2d.getMaxX();
    double maxy = im2d.getMaxY();
    double x1, y1, x2, y2;
    //     double boxsizex = (maxx-minx)/(double)im2d.getXDim();
    //     double boxsizey = (maxy-miny)/(double)im2d.getYDim();
    double boxsizex = im2d.getXCoord(1)- im2d.getXCoord(0);
    double boxsizey = im2d.getYCoord(1)- im2d.getYCoord(0);
    double min; 
    double max; 
    double colorlevel;
    int  red, green, blue;

    if (_state.top().autoscale == true) {
	min = im2d.minValue();
	max = im2d.maxValue();
    }
    else {
	min = _state.top().colorscale_lower;
	max = _state.top().colorscale_upper;
    }

    pl_savestate_r( plotter );

    pl_filltype_r( plotter, 1 );
    pl_pentype_r(  plotter, 0 );

    // plot the image as a grid of boxes:

    for (int i=0; i< im2d.getXDim(); i++) {
	for (int j=0; j< im2d.getYDim(); j++) {

	    x1 = im2d.getXCoord(i) - boxsizex/2.0;
	    y1 = im2d.getYCoord(j) - boxsizey/2.0;
	    x2 = im2d.getXCoord(i) + boxsizex/2.0;
	    y2 = im2d.getYCoord(j) + boxsizey/2.0;

	    colorlevel = ((im2d.getPixel(i,j)-min) / fabs(max-min));

	    getMappedColor( colorlevel, red, green, blue );
  
	    pl_fillcolor_r(  plotter, red, green, blue );
	    pl_fbox_r(  plotter, x1,y1,x2,y2 );

	}
    }

    // now plot a colormap key

    double h = (maxy-miny)/(double)_numkeycolors;
    x1 = maxx + _keyoffset;
    x2 = x1 + _keywidth;
    for (int i=0; i<_numkeycolors; i++) {

	colorlevel = (double)i/(double)_numkeycolors;
	
	getMappedColor( colorlevel, red, green, blue );
	pl_fillcolor_r( plotter, red, green, blue );

	y1 = miny + i*h;
	y2 = miny + (i+1)*h;

	pl_fbox_r( plotter, x1,y1 ,x2,y2 );

    }

    // plot the labels for the colormap

    char text[64];

    pl_filltype_r( plotter, 0 );
    pl_pentype_r(  plotter, 1 );
    pl_flinewidth_r( plotter,0.02 ); 
 
    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, _scalefactor*0.1 );        

    sprintf(text, "%.1f", min);
    pl_fmove_r( plotter, x2+_ticmarksize, miny );
    pl_alabel_r( plotter, 'l','b', text );

    sprintf(text, "%.1f", max);
    pl_fmove_r( plotter, x2+_ticmarksize, maxy );
    pl_alabel_r( plotter, 'l','t', text );



    if (max>0 && min<0) {

	y1 = miny +  fabs(min/(max-min)) * (maxy-miny);
	pl_pencolorname_r( plotter, "black" );
	pl_fline_r( plotter, x1 -_ticmarksize, y1, x2+_ticmarksize, y1 );
	pl_fmove_r( plotter, x2+_ticmarksize, y1 );
	pl_alabel_r( plotter, 'l','c', " 0" );
	    
    }

    pl_fbox_r( plotter, x1, miny ,x2, maxy );

    pl_restorestate_r(  plotter );

}


void 
PlotMaker::
plot( Camera &cam, const Array_t &image, const vector<Pedestal> &peds,
      vector<int> &cleanpixels, CameraPlotType plottype ) {

    double qr;
    char text[128];
    Array_t &xc = cam.xCoords();
    Array_t &yc = cam.yCoords();
    Array_t &rad  = cam.radii();
    double max;
    double value;
    register int i;

    switch (plottype) {
    case CAMERA_IMAGE:
	max = image.max();
	break;

    case CAMERA_PEDS:
	max = -1e10;
	for (i=0; i<peds.size(); i++) {
	    if (peds[i].pedestal > max)
		max = peds[i].pedestal;
	}
	break;

    case CAMERA_PEDDISPS:
	max = -1e10;
	for (i=0; i<peds.size(); i++) {
	    if (peds[i].dispersion > max)
		max = peds[i].dispersion;
	}
	break;

    default:
	max = image.max();
	break;
    }

    pl_savestate_r( plotter );
    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r(plotter, _scalefactor*0.07);

    // ===============================================================
    // Set View parameters

    pl_flinewidth_r( plotter, 0.005);  // line thickness in user coords 
    pl_bgcolorname_r(plotter, "white");	// background color for the window
    pl_pencolorname_r( plotter, "grey60"); 

    // ===============================================================
    // draw camera

    pl_filltype_r(plotter, 0);			    // don't fill
    for(int i=0; i<cam.getNumPixels(); i++) {
	pl_fellipse_r(plotter, xc[i], yc[i], rad[i], rad[i],0);
    }

    // ===============================================================
    // draw signal

    pl_filltype_r( plotter, 1 );
    for(int i=0; i<cam.getNumPixels(); i++) {

	pl_fmove_r( plotter, 0.0, _maxy - 0.05);
	
	// --------------------------------------------------
	// choose value based on type:

	switch (plottype) {
	case CAMERA_IMAGE:
	    value = image[i];
	    break;
	case CAMERA_PEDS:
	    value = peds[i].pedestal;
	    break;
	case CAMERA_PEDDISPS:
	    value = peds[i].dispersion;
	    break;
	case CAMERA_TUBENUMS:
	    value = 0;
	default:
	    value = image[i];
	    break;
	    
	}

	// --------------------------------------------------
	// choose a color for the plot

	if (value >= 0) {
	    pl_pencolorname_r(plotter, "LightSkyBlue"); 
	    pl_fillcolorname_r(plotter, "LightSkyBlue");
	}
	else { 
	    pl_pencolorname_r( plotter, "goldenrod"); 
	    pl_fillcolorname_r( plotter, "goldenrod");
	}
	

	// --------------------------------------------------
	// apply scaling: 

	if (_pixelscale ==-1) {
	    // Autoscale to max
	    qr =  rad[i]*fabs(value)/max;
	}
	else {
	    // Fixed scale to _pixelscale
	    qr = rad[i]*fabs(value)/_pixelscale;
	    if (qr > rad[i]) qr = rad[i]; // saturated
	}


	// --------------------------------------------------
	// deal with tubes which are off

	if (peds[i].type != Pedestal::GOOD) {
	    if (peds[i].type == Pedestal::TUBEOFF) // off
		pl_fillcolorname_r( plotter, "gray60");
	    else if (peds[i].type == Pedestal::STAR)
		pl_fillcolorname_r( plotter, "gray90");
	    else
		pl_fillcolorname_r(plotter, "white");
	    pl_pencolorname_r(plotter, "gray88"); 
	    qr = rad[i];
	}
	

	// plot the value
	if (value>0)
	    pl_fellipse_r( plotter, xc[i],yc[i],qr,qr,0); 
	    
	
    }

    // ===============================================================
    // mark picture tubes

    if (plottype == CAMERA_IMAGE ) {
	for (vector<int>::iterator iter=cleanpixels.begin();
	     iter!=cleanpixels.end(); iter++){
	    pl_pencolorname_r( plotter, "darkblue");
	    pl_fillcolorname_r( plotter, "darkblue");
	    if (_pixelscale==-1) qr = fabs(rad[*iter]*image[*iter]/max);
	    else qr = rad[*iter]*image[*iter]/_pixelscale; 
	    if (qr > rad[*iter]) qr = rad[*iter];
	    pl_fellipse_r( plotter, xc[*iter],yc[*iter],qr,qr,0);  
	}

	// display the signal value if over a threshold
	
	for (int i=0; i<cam.getNumPixels(); i++) {
	    if(image[i]>max*0.75 && rad[i]) {
		sprintf(text,"%d", (int) image[i]);
		if (pl_flabelwidth_r( plotter,text ) < rad[i]*1.2){
		    pl_pencolorname_r( plotter, "white"); 
		    pl_fmove_r( plotter,  xc[i],yc[i] );	    
		    pl_alabel_r(plotter, 'c','c', text );
		}
	    }
	}
    }

    // ===============================================================
    // Display the tube numbers if requested

    pl_fontname_r( plotter, "HersheySans");
    pl_ffontsize_r( plotter, _scalefactor*0.05 );

    if (plottype == CAMERA_TUBENUMS) {

	pl_pencolorname_r( plotter, "black"); 
	for (int i=0; i<cam.getNumPixels(); i++) {
	    sprintf(text,"%d", i);
	    pl_fmove_r( plotter,  xc[i],yc[i] );	    
	    pl_alabel_r(plotter, 'c','c', text );
	}
    }
    
    pl_restorestate_r( plotter );

}

/**
 * Plot a curve specified by a set of x,y points
 */
void 
PlotMaker::
plot( vector<double> &xpoints, vector<double> &ypoints ) {

    if (xpoints.size() != ypoints.size()){
	throw MildAnalysisException( "Couldn't plot x and y "
				     "vectors with different sizes");
    }

    if (xpoints.size() ==0 || ypoints.size()==0){
	throw MildAnalysisException( "Couldn't plot curve with no points ");
    }

    pl_savestate_r( plotter );
    pl_flinewidth_r( plotter, 0.009 );

    setLineMod();
    
    pl_fmove_r( plotter, xpoints[0], ypoints[0] );

    for (int i=1; i< (int)xpoints.size(); i++) {
	
	pl_fcont_r( plotter, xpoints[i], ypoints[i] );

    }

    pl_endpath_r( plotter );
    pl_restorestate_r( plotter );

}


/**
 * plot a text label 
 */
void 
PlotMaker::
plot( string text, double x, double y ){
    
    pl_savestate_r( plotter );

    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, _state.top().labelsize );
    pl_ftextangle_r( plotter, _state.top().labelangle );
    pl_pencolorname_r( plotter,_state.top().labelcolor.c_str() );
    pl_fmove_r( plotter, x,y );
    pl_alabel_r( plotter, 'c','c', text.c_str() );
    pl_restorestate_r( plotter );

}


/**
 * Plot an array image
 */
void
PlotMaker::
plot( TelescopeArray &array ) {

    double x,y,r;

    pl_savestate_r( plotter );

    pl_pencolorname_r( plotter,"black" );

    for (int i=0; i<array.getNumTelescopes(); i++){
	r = 0.5*(array.getCamera(i)->getMaxX() 
		 - array.getCamera(i)->getMinX());

	x = array.getLocationX(i);
	y = array.getLocationY(i);
	
	pl_fellipse_r( plotter, x,y, r,r, 0 );
	
    }


    pl_restorestate_r( plotter );

}

void
PlotMaker::
pushState() {
    pl_savestate_r( plotter );
    _state.push( _state.top() );
}

void
PlotMaker::
popState() {
    pl_restorestate_r( plotter );
    _state.pop();
};

void 
PlotMaker::
rotate( double theta ) {

    pl_frotate_r( plotter, theta );

}


void 
PlotMaker::
plot( HillasParameterization &p ) {

    const double step = 0.11;
    char text[128];

    pl_savestate_r( plotter );

    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r(plotter, 0.07);

    // draw the centeroid

    pl_pencolorname_r( plotter, "red");
    pl_fillcolorname_r( plotter,"red"); 
    pl_fmarker_r(  plotter,p.centroid.x, p.centroid.y, 9,0.1);

    // draw the points or origin

    pl_pencolorname_r(  plotter,"DarkGreen");
    pl_fmarker_r(plotter,p.point_of_origin_a.x, p.point_of_origin_a.y, 10,0.4);
    pl_fmarker_r(plotter,p.point_of_origin_b.x, p.point_of_origin_b.y, 10,0.3);
    
    // draw the hillas ellipse
    pl_pencolorname_r( plotter,"red");
    pl_flinewidth_r( plotter,0.02); // line thickness in user coordinates 
    pl_filltype_r( plotter,0);  
    pl_pencolorname_r( plotter,"red"); 
    pl_fellipse_r( plotter,p.centroid.x, p.centroid.y, 
		   p.length, 
		   p.width, (p.psi * 180.0/M_PI) );
    
    // draw the labels

    pl_fmove_r(  plotter,_minx+0.1,_maxy-step*1 );
    sprintf(text, "max1,2,3: %5.2f %5.2f %5.2f", p.max[0], p.max[1], p.max[2]);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*2 );
    sprintf(text, "alpha: %3.3g", p.alpha *180.0/M_PI);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*3 );
    sprintf(text, "miss: %3.3g", p.miss);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*4 );
    sprintf(text, "dist: %3.3g", p.distance);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*5 );
    sprintf(text, "length: %3.3g", p.length);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*6 );
    sprintf(text, "width: %3.3g", p.width);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*7 );
    sprintf(text, "size: %.2f", p.size);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*8 );
    sprintf(text, "psi: %3.3g deg", p.psi*180.0/M_PI);
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1,_maxy-step*9 );
    sprintf(text, "asym: %3.3g", p.asymmetry);
    pl_label_r( plotter, text );

    if (p.invalid != 0){
	pl_fmove_r( plotter, _minx+0.1,_maxy-step*10 );
	sprintf(text, "INVALID: %d ",(int) p.invalid);
	pl_label_r( plotter, text );
    }

    pl_ffontsize_r( plotter,0.09);
    pl_pencolorname_r(  plotter,"darkgreen" );
    pl_fmove_r( plotter, _maxx-0.1, _maxy-step );
    sprintf(text,"%d", p.event_number);
    pl_alabel_r(  plotter,'r','c',text);

    pl_pencolorname_r(  plotter,"DarkKhaki" );
    pl_fmove_r( plotter, _maxx-0.1, _maxy-2*step );
    sprintf(text,"T%2d", p.telescope_id);
    pl_alabel_r(  plotter,'r','c',text);

    // plot center marker
    pl_flinewidth_r( plotter,0.01); // line thickness in user coordinates 
    pl_fmarker_r(  plotter,0.0, 0.0,9,0.2);

    pl_restorestate_r( plotter );    

}

void 
PlotMaker::
plot( MuonParameterization &p ) {

    const double step = 0.11;
    char text[128];

    pl_savestate_r( plotter );

    // draw the other ring centers
    
    pl_pencolorname_r( plotter,"goldenrod1");
    for (int i=0; i<p.nplotpoints; i++) {
	pl_fmarker_r(  plotter,p.plot[i].x, p.plot[i].y, 10, 0.08 );
    }

    // set color

    pl_pencolorname_r(  plotter,"orange");
    pl_fillcolorname_r(  plotter,"orange");     

    // draw the muon center position

    pl_fmarker_r(  plotter,p.ringcenter.x, p.ringcenter.y, 2, 0.4 );
    
    // draw the ring itself
   
    pl_linemod_r( plotter, "shortdashed" );
    pl_filltype_r( plotter,0);		// don't fill
    pl_flinewidth_r( plotter,0.02);	// line thickness in user coordinates 
    pl_fcircle_r(  plotter,p.ringcenter.x, p.ringcenter.y, p.radius );

    // draw the arc (philo to phihi) region:

    double x0 = p.radius*cos(p.philo) + p.ringcenter.x;
    double y0 = p.radius*sin(p.philo) + p.ringcenter.y;
    double x1 = p.radius*cos(p.phihi) + p.ringcenter.x;
    double y1 = p.radius*sin(p.phihi) + p.ringcenter.y;
    double xm = p.radius*cos(p.phimid) + p.ringcenter.x;
    double ym = p.radius*sin(p.phimid) + p.ringcenter.y;

    pl_fline_r( plotter, p.ringcenter.x, p.ringcenter.y, xm,ym );
    pl_flinewidth_r( plotter,0.05 );	
    pl_linemod_r( plotter, "solid" );
    pl_farc_r( plotter, p.ringcenter.x, p.ringcenter.y, x0,y0,x1,y1 ); 

    pl_pencolorname_r( plotter,"black");
    pl_fmarker_r( plotter, p.xcs, p.ycs, 6, 0.2 );
    pl_pencolorname_r( plotter,"orange");

    // display the relevant parameters

    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r(plotter, 0.07);    
    
    pl_fmove_r( plotter, _minx+0.1, _miny+step*1 );
    sprintf( text, "S/Arclen: %.2f", p.soal );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*2 );
    sprintf( text, "Radius: %f", p.radius );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*3 );
    sprintf( text, "ArcStrength: %f", p.arcstrength );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*4 );
    sprintf( text, "RingFrac: %.2f %%", p.ringfrac * 100 );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*5 );
    sprintf( text, "MuSkew: %f", p.muskew );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*6 );
    sprintf( text, "MuGain: %f", p.mugain );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*7 );
    sprintf( text, "gain: %f", p.gain );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*8 );
    sprintf( text, "muonness: %.2f", p.muonness );
    pl_label_r( plotter, text );

    pl_fmove_r( plotter, _minx+0.1, _miny+step*9 );
    sprintf( text, "arclength: %.2f", p.arclen );
    pl_label_r( plotter, text );


    pl_restorestate_r( plotter );

    

}


/**
 * Returns the rgb color for the specified grey value (in range 0..1)
 *
 * \param value grey value which should be mapped to a color (0..1)
 * \param red returned color value (0..65535)
 * \param green returned color value (0..65535)
 * \param blue returned color value (0..65535)
 */
void
PlotMaker::getMappedColor( double value, int &red, int &green, int &blue ) {

  
    switch ( _state.top().colormap) {
	
    case COLORMAP_BLUES:
	red = static_cast<int>( (value*value)*65535.0);
	green = static_cast<int>(value*65535.0);
	blue = static_cast<int>(sin(M_PI_2*value)*65535.0);
	break;

    case COLORMAP_YELLOWGREEN:
	red = static_cast<int>( value*(pow(sin(M_PI_2*value),2) + 
				       pow(cos(3*M_PI_2*value),2) )*65535.0);
	green = static_cast<int>(sin(M_PI_2*value)* 65535.0);
	blue = static_cast<int>( (pow(sin(3*M_PI_2*value),2)/2.0 + value/2.0)
				 *65535.0);	
	break;

    case COLORMAP_REDBLUE:
	red   = static_cast<int>((sqrt(value))*65535.0);
	green = static_cast<int>((pow(value,3))*65535.0);
	blue  = static_cast<int>((sin(2*M_PI*value))*65535.0);
	break;

    case COLORMAP_INVERSEGREY:
	red   = static_cast<int>((1.0-value)*65535.0);
	green = static_cast<int>((1.0-value)*65535.0);
	blue  = static_cast<int>((1.0-value)*65535.0);

	break;

    case COLORMAP_GREY:
	red   = static_cast<int>((value)*65535.0);
	green = static_cast<int>((value)*65535.0);
	blue  = static_cast<int>((value)*65535.0);
	break;

    case COLORMAP_RAINBOW:
	red   = static_cast<int>((abs(2*value-0.5))*65535.0);
	green = static_cast<int>((sin(M_PI*value))*65535.0);
	blue  = static_cast<int>((cos(M_PI/2.0*value))*65535.0);
	break;

    case COLORMAP_HOT:
	red   = static_cast<int>((3*value)*65535.0);
	green = static_cast<int>((3*value-1)*65535.0);
	blue  = static_cast<int>((3*value-2)*65535.0);
	break;

    case COLORMAP_PRINTABLE:
	red   = static_cast<int>((value/0.32-0.78125)*65535.0);
	green = static_cast<int>((2*value-0.84)*65535.0);
	blue  = static_cast<int>((4*value)*65535.0);
	break;

    default:
	throw MildAnalysisException("unknown color map value");

    };


    if (red > 65535)   red   = 65535;
    if (green > 65535) green = 65535;
    if (blue > 65535)  blue  = 65535;

    if (red < 0)    red   = 0;
    if (green < 0 ) green = 0;
    if (blue < 0 )  blue  = 0;



}


void 
PlotMaker::
plot( Star &star, int symbol ) {
    
    double size = (_maxmagnitude - star.magnitude)/40.0;
    if (size<=0) return;

    pl_savestate_r( plotter );
    pl_fontname_r( plotter, "HersheySerif");
    pl_ffontsize_r( plotter, 0.06 );    
    pl_pencolorname_r( plotter, "yellow" );
    pl_fmarker_r(plotter,star.tx,star.ty, symbol ,size );
    if ( star.catalog_id == 0 ) {
	pl_fmove_r( plotter, star.tx, star.ty-size );
	pl_alabel_r( plotter, 'c','t',star.name.c_str() );
    }
    pl_restorestate_r( plotter );
    
}


/**
 * Plot contours for the given Image2D
 */
void 
PlotMaker::
contourPlot( Image2D &im2d, double zmin, double zmax, double dz ) {

    int nlevels = static_cast<int>( (zmax-zmin)/dz );
    int nx = im2d.getXDim();
    int ny = im2d.getYDim();
    double x00,x01,x10,x11;
    double y00,y01,y10,y11;
    double z00,z01,z10,z11;
    double zc;  
    int k, red,green,blue;
    bool iny, inx;
    double x[5], y[5];
    double xin, yin; 		// x, y intercepts

    ContourSegment segment;

    if (fabs(zmin-zmax) < 1e-10) {
	cout << "contourPlot: too small a contour range! (" 
	     << zmin << ","<<zmax<<")"<<endl;
	return;
     }

    // Set up the colors for plotting 

    pl_savestate_r( plotter );
    //    pl_flinewidth_r( plotter, 0.009 );

    setLineMod();
        
    // First we need to figure out where all the contours go, then
    // we'll plot them
    
    // For each grid box
    for ( int yi=0; yi<nx-1; yi++ ) {
	for (int xi=0; xi<ny-1; xi++ ) {
	    
	    // get the coordinate centers of the current block of four
	    // pixels (just to make them easier to type)

	    x00 = im2d.getXCoord( xi   );
	    x10 = im2d.getXCoord( xi+1 );
	    x01 = im2d.getXCoord( xi   );
	    x11 = im2d.getXCoord( xi+1 );

	    y00 = im2d.getYCoord( yi );
	    y10 = im2d.getYCoord( yi );
	    y01 = im2d.getYCoord( yi+1 );
	    y11 = im2d.getYCoord( yi+1 );

	    z00 = im2d.getPixel( xi   , yi );
	    z10 = im2d.getPixel( xi+1 , yi );
	    z01 = im2d.getPixel( xi   , yi+1 );
	    z11 = im2d.getPixel( xi+1 , yi+1 );
	    

	    // For each contour j:
	    
	    zc = zmin;
	    for (int j=0; j<nlevels; j++) {
		
		k = 0;

		// Check if the contour passes through the 4 pixel block

		iny = ( ((z00<=zc)&&(zc<=z01))||((z00>=zc)&&(zc>=z01)) ) ;
		iny = iny ||( ((z10<=zc)&&(zc<=z11))||((z10>=zc)&&(zc>=z11)) );
		inx = ( ((z00<=zc)&&(zc<=z10))||((z00>=zc)&&(zc>=z10)) ) ;
		inx = inx ||( ((z01<=zc)&&(zc<=z11))||((z01>=zc)&&(zc>=z11)) );
		
		if (inx || iny) {
		    
		    xin = xIntercept(z00,z01,z10,z11,zc,0.0);
		    yin = yIntercept(z00,z01,z10,z11,zc,0.0);
		    if ( xin > 0 ){
			x[k] = x00 + xin*(x10-x00);
			y[k] = y00 + xin*(y10-y00);
			k++;
		    }
		    if( yin > 0.00) {
			x[k] = x00+yin*(x01-x00) ;
			y[k] = y00+yin*(y01-y00) ;
			k++;
		    }

		    xin=xIntercept(z00,z01,z10,z11,zc,1.0);
		    yin=yIntercept(z00,z01,z10,z11,zc,1.0);
		    if ( xin > 0.00 ){
			x[k] = x01+xin*(x11-x01) ;
			y[k] = y01+xin*(y11-y01) ;
			k++;
		    }
		    if( yin > 0.00) {
			x[k] = x10+yin*(x11-x10) ;
			y[k] = y10+yin*(y11-y10) ;
			k++;
		    }
		    
		    for (int i=0; i<k-1; i++ ) {
			
			// contour segment is: x[i],y[i] - x[i+1],y[i+1]

			segment.start.x = x[i];
			segment.start.y = y[i];
			segment.end.x   = x[i+1];
			segment.end.y   = y[i+1];

			// plot it!
			getMappedColor(((double)j/(double)nlevels), 
				       red,green,blue );
			if ( _state.top().contourcolor == "MAPPED")
			    pl_pencolor_r( plotter, red,green,blue  );

			pl_fline_r( plotter, 
				    segment.start.x,
				    segment.start.y,
				    segment.end.x,
				    segment.end.y );	
		    }
		    
		}

		// move to the next contour level
		zc += dz; 

	    }
	}
    }

    // Now, we know all the contours, so on to plotting...
    
    pl_restorestate_r( plotter );
    

}


/**
 * Plot an ellipse centered on x,y with radii r1 and r2 and angle theta
 */
void 
PlotMaker::
plotEllipse( double x, double y, double r1, double r2,double theta) {

    pl_savestate_r( plotter );
    pl_filltype_r( plotter, 0 );
    pl_flinewidth_r( plotter, 0.009 );
    setLineMod();
    pl_fellipse_r( plotter, x,y, r1,r2, theta );
    pl_restorestate_r( plotter );    

}


/**
 * Plot RA/DEC contours and labels
 */
void 
PlotMaker::
plotRADecGrid( Image2D &ragrid, Image2D &decgrid ) {
    
    // get RA/DEC of center pixel:
    
    int xcenterpix = (ragrid.getXDim()/2);
    int ycenterpix = (ragrid.getYDim()/2);

    double ra = ragrid.getPixel( xcenterpix, ycenterpix ) *M_PI/12.0;
    double dec = decgrid.getPixel( xcenterpix, ycenterpix )*M_PI/180.0;


    // get the nice contour values:

    double xcont[NCONTSTEPS];
    double ycont[NCONTSTEPS];
    double nice_dec[NCONTSTEPS];
    double nice_ra[NCONTSTEPS];
    int    drawx[NCONTSTEPS];
    int    drawy[NCONTSTEPS];
    char   ra_str[NCONTSTEPS][128];
    char   dec_str[NCONTSTEPS][128];

    niceRaDecContours( ra, dec, 
		       xcont, ycont,
		       nice_ra, nice_dec,
		       drawx, drawy,
		       ra_str, dec_str );

    for (int i=0; i<NCONTSTEPS; i++) {


	// Plot the RA contours and labels
	
	contourPlot( ragrid, 	// kind of a hack to plot just one contour
		     nice_ra[i]*12.0/M_PI, 
		     nice_ra[i]*12.0/M_PI + 4.0, 3.9 );
	
	if ( drawx[i] ){
	    plot( ra_str[i], xcont[i], _miny + _ticmarksize*2 );
	}


	// plot the dec contours and labels 

	contourPlot( decgrid, 
		     nice_dec[i]*180.0/M_PI, 
		     nice_dec[i]*180.0/M_PI + 4.0, 3.9 );

	if ( drawy[i] ){
	    pushState();
	    _state.top().labelangle = 90.0;
	    plot( dec_str[i], _minx+_ticmarksize*2, ycont[i] );
	    popState();
	}	
	
	
    }

}


/**
 * helper func to switch line style based on current value
 */
void 
PlotMaker::
setLineMod() {
    
    switch (_state.top().contourlinestyle) {

    case LINE_SOLID:
	pl_linemod_r( plotter, "solid" );
	break;

    case LINE_DOTTED:
	pl_linemod_r( plotter, "dotted" );
	break;	
	
    case  LINE_DASHED:
	pl_linemod_r( plotter, "shortdashed" );
	break;

    default:
	pl_linemod_r( plotter, "solid" );
	break;
    }

    if ( _state.top().contourcolor != "MAPPED" )
	pl_pencolorname_r( plotter,_state.top().contourcolor.c_str() );
    
    pl_flinewidth_r( plotter,_state.top().contourwidth );

}

void 
PlotMaker::
plotMarker( double x, double y ) {
    
    pl_savestate_r( plotter );    
    pl_pencolorname_r( plotter,_state.top().markercolor.c_str() );
    pl_fmarker_r( plotter, x,y,_state.top().markertype,_state.top().markersize );
    pl_restorestate_r( plotter);

}


double	
PlotMaker::
yIntercept(double z00, double z01, double z10, double z11, 
	   double zc, double xr) {
    double	yr ;
    double	t1,t2 ;
    
    t1 = zc-xr*(z10-z00)-z00 ;
    t2 = xr*(z11-z01-z10+z00)+z01-z00 ;
    
    yr=(double)t1/(double)t2 ;
    if((yr>=0.00) && (yr<=1.00)) 
	return(yr) ;
    else
	return(-1.00) ;
    
}

double
PlotMaker::	
xIntercept(double z00, double z01, double z10, double z11, 
	   double zc, double yr) {

    double xr ;
    double t1,t2 ;
    
    t1 = zc-yr*(z01-z00)-z00 ;
    t2 = yr*(z11-z01-z10+z00)+z10-z00 ;
    
    xr = (double)t1/(double)t2 ;
    if((xr>=0.00) && (xr<=1.00))
	return(xr) ;
    else
	return(-1.00) ;
    
}


/**
 * 	nice_contours()
 * 	J. Buckley
 *
 * 	Does calculations required to draw nice RA and DEC contours
 * 	around a source at the given ra and dec.
 *
 * 	Inputs:
 * 	ra	Right ascension in hhmmss.s
 * 	dec	Declination in ddmmss.s
 * 	xcont[]	Array of values for text position for RA labels
 * 	        along the x-axis in camera/tangential coordinates (deg).
 * 	ycont[] Array of values giving the text position for the DEC
 * 	        labels along the yaxis in camera coordinates (deg).
 * 	dec_nice[]	Nice declination values (in radians)
 * 	ra_nice[]	Nice RA values in radiatns
 * 	ra_str[] Array of strings to place at xcont[] positions along
 * 	 	the x-axis giving the formatted right ascension to the
 * 	 	nearest minute
 * 	dec_str[] Array of strings to place at ycont[] positions along
 * 	        the y-axis containing formatted declination valures to
 * 	        the nearest arcmin
 *
 *      drawx[] and draw[y] are arrays of flags for which contours to draw
 */
void 
PlotMaker::
niceRaDecContours(double ra,  double dec, 
		  double xcont[],double ycont[], 
		  double ra_nice[], double dec_nice[], 
		  int drawx[], int drawy[], 
		  char ra_str[17][128],
		  char dec_str[17][128]) {

    register int i, j, k ;
    double hh, mm, ss ;
    double fraz ;	/* RA in radians but truncated to the nearest min */
    double radeg ;	/* RA in decimal degrees */
    double dd, arcmm, archm, arcss ;
    double fdecz ;  /* DEC in radians but truncated to the nearest arcmin */
    double decdeg ; /* DEC in decimal degrees */
    double ratmp, dectmp,dectmp2,denom;
    static int	nbins = 6 ;
    static double xoff_dec[] = {-32.0, -16.0, -8.0,-4.0, -8.0, -16.0, -32.0} ;
    static double xstep_dec[] = {4.0, 2.0, 1.0, 0.5, 1.0, 2.0, 4.0} ;
    static double dec_bins[] = {-90.1, -73.0, 
				-67.0 ,-45, 
				45.0, 67.0, 
				73.0, 90.1};
    double xoff ;	/* Starting RA offset in degrees */
    double dx ;	/* RA step in degrees */
    char tempstra[100] ;
    double frahh[NCONTSTEPS+1] ;
    double framm[NCONTSTEPS+1] ;
    int irahh[NCONTSTEPS+1] ;
    int iramm[NCONTSTEPS+1] ;
    double fdecdd[NCONTSTEPS+1] ;
    double fdecmm[NCONTSTEPS+1] ;
    int idecdd[NCONTSTEPS+1] ;
    int	idecmm[NCONTSTEPS+1] ;
    
    // KPK: this routine (which was taken from the old quicklook kumac
    // files originally assumed that ra and dec were in "formatted"
    // format (ie, hhmmss.ss) instead of radians.  We use radians, so
    // we need to first convert to the hhmmss format!  (confusing, but
    // easier than modifying this code!)

    double ra_h, ra_m, ra_s;
    double dec_d, dec_m, dec_s;

    radToHMS( ra, ra_h,ra_m,ra_s );
    radToDMS( dec, dec_d,dec_m,dec_s );

    ra = (((int)ra_h)*10000 +
	  ((int)ra_m)*100 +
	  ra_s);
    
    dec = (((int)dec_d)*10000 +
	   ((int)dec_m)*100 +
	   dec_s);
	   
    // done with conversion!

    hh = int(ra/10000.) ;
    mm = int((ra-hh*10000)/100.) ;
    ss = ra - hh*10000 - mm*100. ;

    fraz = (hh+mm/60.)*15.0/(180.0/M_PI) ;
    radeg = (hh+mm/60.+ss/3600.)*15.0 ;

    dd = int(dec/10000.) ;
    arcmm = int((dec-dd*10000.)/100.) ;
    archm = int((dec-dd*10000.)/(100.*30.)) ;
    arcss = dec-dd*10000.-arcmm*100. ;
    fdecz = (dd+0.5*archm)/(180.0/M_PI) ;  /* DEC in rad truncated to min */
    decdeg = dd+arcmm/60.0+arcss/3600. ;  /* DEC in degrees */

    for(i=0; i<nbins; i++) {	/* Find good values for RA given DEC */
	if((decdeg > dec_bins[i]) && (decdeg < dec_bins[i+1])) {
	    xoff = xoff_dec[i] ;
	    dx = xstep_dec[i] ;
	}
    }
    ra = radeg/(180.0/M_PI) ;	
    dec = decdeg/(180.0/M_PI) ;
      
    for(i=0; i<NCONTSTEPS; i++) {
        ratmp = fraz*(180.0/M_PI)+xoff ;
        dectmp = fdecz*(180.0/M_PI)+xoff ;
        frahh[i] = floor(ratmp*24./360.+0.0001) ;
        framm[i] = floor((ratmp*24./360.-frahh[i])*60.+0.0001) ;
        irahh[i] = (int)floor(frahh[i]+0.0001) ;
        iramm[i] = (int)floor(framm[i]+0.0001) ;

	if (dectmp<0) dectmp2 = dectmp*-1;
	else dectmp2 = dectmp;

	fdecdd[i] = floor((dectmp2+0.0001)) ;
	fdecmm[i] = floor((dectmp2-fdecdd[i])*60.+0.0001) ;
	idecdd[i] = (int)floor(fdecdd[i]+0.001) ;
	idecmm[i] = (int)floor(fdecmm[i]+0.001) ;

	if (dectmp<0) {
	    fdecdd[i] *= -1;
	    idecdd[i] *= -1;
	}

        xoff += dx ;
	
        dec_nice[i] = dectmp /= (180.0/M_PI) ;
        ra_nice[i] = ratmp /= (180.0/M_PI) ;

        denom = sin(dec)*sin(dectmp)+cos(dec)*cos(dectmp)*cos(ratmp-ra) ;
        if(denom > 0.0) {
	    xcont[i]=-(180.0/M_PI)*cos(dectmp)*sin(ratmp-ra)/denom ;
	    ycont[i]=(180.0/M_PI)*(cos(dec)*sin(dectmp)-
				  sin(dec)*cos(dectmp)*cos(ratmp-ra))/denom ;
	    if((xcont[i] > _minx) && (xcont[i] < _maxx)) 
		drawx[i] = 1 ;
	    else 
		drawx[i] = 0 ;
	    if((ycont[i] > _miny) && (ycont[i] < _maxy)) 
		drawy[i] = 1 ;
	    else 
		drawy[i] = 0 ;
	} 
        else {
	    xcont[i] = 1000.0 ;
	    ycont[i] = 1000.0 ;
	    drawx[i] = 0 ;
	    drawy[i] = 0 ;
	}
	
	sprintf(ra_str[i], "%2d\\sph\\ep%2d\\spm\\ep", irahh[i], iramm[i]) ;
	sprintf(dec_str[i], "%3d\\de%02d'", idecdd[i], idecmm[i]) ;
    } 

} 
                           
    
