//
// MuonImageAnalyzer.cpp  - ImageAnalyzer for looking for muon rings
//
// Jim Buckley <buckley@wuphys.wustl.edu>
// Karl Kosack <kosack@hbar.wustl.edu>
//

#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <valarray>
#include <assert.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_math.h>

#include "Camera.h"
#include "Exceptions.h"
#include "Types.h"
#include "ImageAnalyzer.h"
#include "MuonImageAnalyzer.h"
#include "Image2D.h"
#include "PlotMaker.h"
#include "Histogram.h"

using namespace std;

/**
 * Create a new MuonImageAnalyzer
 */
MuonImageAnalyzer::MuonImageAnalyzer( Camera *cam )
    : ImageAnalyzer(cam), _plot_points_is_enabled(false), 
      _calibration_is_enabled(true)
{

    register int i,j;


    //    pm = new PlotMaker( "X","x.out" );


    // Set the sizes of all the lookup tables and grids

    _single_pmt_mask.resize(getNumPixels()) ;
    
//     cout << "DEBUG: resizing centers to "<< getNumPixels() << "^3 size="
// 	 << pow((double)getNumPixels(),3)*sizeof(Triplet) 
// 	 << " sizeof triplet is "<< sizeof(Triplet)
// 	 << endl;

//     _centers.resize(getNumPixels()) ;   // Allocate the 3d centers array 
//     for(i=0;i<getNumPixels();i++) {
// 	_centers[i].resize(getNumPixels()) ;
// 	for(j=0;j<getNumPixels();j++) {
//             _centers[i][j].resize(getNumPixels()) ;
// 	}
//     }
    

    _pmtmask.resize(NX) ;   // Now initalize all of the arrays that 
    _ntrip.resize(NX) ;     // cover the NX x NY 2-d grid covering 
    _xsum.resize(NX) ;      // the camera.                         
    _ysum.resize(NX) ;
    _rsum.resize(NX) ;
    
    for(i=0;i<NX;i++) {
	_pmtmask[i].resize(NY) ;
	_ntrip[i].resize(NY) ;
	_xsum[i].resize(NY) ;
	_ysum[i].resize(NY) ;
	_rsum[i].resize(NY) ;
	for(j=0;j<NY;j++) {
            _pmtmask[i][j].resize(_nwords = getNumPixels()/sizeof(int)) ;
	}
    }
        
    // Initialize a structure containing the bitmask for each PMT.
    // Since there are a large number of PMTs in the camera, divide
    // the mask up into a number of 32 bit words
     
    for(i=0; i<getNumPixels(); i++) {
	_single_pmt_mask[i].word = i/sizeof(int) ;
	_single_pmt_mask[i].bit = 0x00000001 << (i % sizeof(int)) ;
    }
    
    // Loop through the N*(N-1)*(N-2)/(3*2*1) triplets, where N=nadc
    // vals and calculate the (x,y) coordinate of the ring center for
    // each triplet, the index (ix,iy) of the ring-center grid and the
    // radius r of the ring defined by these three points.

//     int i1,i2,i3;
    
//     for(i1=0; i1<getNumPixels(); i1++) {
// 	for(i2=i1+1; i2<getNumPixels(); i2++) {
//             for(i3=i2+1; i3<getNumPixels(); i3++) {

// 		getCenter( i1,i2,i3, _centers[i1][i2][i3] );
		
//             } // end for i3 
// 	} // end for i2 
//     } // end for i1 
    

}

/**
 * get the perpendicular bisector center position.
 */ 
inline void 
MuonImageAnalyzer::
getCenter( int i1, int i2, int i3, Triplet & t ) {

    double denom, numx, numy;
    double x0,y0,x1,y1,x2,y2,x3,y3;
    double  x1sq,y1sq,x2sq,y2sq,x3sq,y3sq ;
    double rsq;

    x1 = _x[i1] ; x2 = _x[i2] ; x3 = _x[i3] ;
    y1 = _y[i1] ; y2 = _y[i2] ; y3 = _y[i3] ;
    
    x1sq = _x2[i1] ; x2sq = _x2[i2] ; x3sq = _x2[i3] ;
    y1sq = _y2[i1] ; y2sq = _y2[i2] ; y3sq = _y2[i3] ;
    
    denom = 2.0*((x2-x3)*(y1-y2)-(x1-x2)*(y2-y3)) ;
    
    t.ix = t.iy = -1;

    if(fabs(denom) > 1.0e-10) {
	numx = ((x2sq-x3sq+y2sq-y3sq)*(y1-y2)
		-(x1sq-x2sq+y1sq-y2sq)*(y2-y3)) ;
	x0 = numx/denom ;
	t.ix = static_cast<short>(NXD2 + floor(x0/DX));
	
	if((t.ix > 0)&&(t.ix <= NX)) {
	    numy = ((x1sq-x2sq+y1sq-y2sq)*(x2-x3)
		    -(x2sq-x3sq+y2sq-y3sq)*(x1-x2)) ;
	    y0 = numy/denom ;
	    t.iy = static_cast<short>(NYD2 + floor(y0/DY));
	    
	    if((t.iy > 0)&&(t.iy <= NY)) {
		rsq = (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0) ;
		if((rsq > 0.0)&&(rsq < RSQMAX)) {
		    t.x0 = x0 ;
		    t.y0 = y0 ;
		    t.r = sqrt(rsq) ;
		} // end if rsq 
		else { 
		    // Flag that value is out of range 
		    t.ix = -1 ;
		    t.iy = -1 ;
		}
	    } //end if iy 
	    else {
		t.iy = -1 ;
		t.ix = -1 ;
	    }
	} // end if ix 
	else {
	    t.ix = -1 ;
	    t.iy = -1;
	}
    } // end if denom 

    
}



/**
 * Find muon arc and calculate parameters.
 *
 * \param image the image
 * \param cleanpixels vector of clean pixel numbers from 
 */
void 
MuonImageAnalyzer::
parameterize( const Array_t &image, std::vector<int> &cleanpixels ){

    int i,j,k,i1,i2,i3,n,ix,iy;
    int nbox;
    int ixmax,iymax,nmax ;
    double dnmax;
    double phsum;
    int nph;
    vector<unsigned int> totpmtmask ;
    double dxa, dya, stot, sarc;
    double wsum, wphi, wphisq, wt, phi;
    double wt_sum;
    double ringfrac1, ringfrac2;
    double xcs, ycs, dphi;
    double size_clean, size_raw;
    int npix, narc;
    int npict;
    int word1, word2, word3, tubenum1, tubenum2, tubenum3;
    int tubenum, wordnum;
    Triplet center;

    string str;
    Histogram phihist( 15, 0.0, M_PI*2.0, "Muon Phi" );

    // KPK: added this to calculate area of pixel.  It assumes no
    // spacing between pixels, which is fair due to the light cones.
    // (just pi * r^2)
    double apix =  (pow(_x[1]-_x[0],2) + pow(_y[1]-_y[0], 2))/4.0 * M_PI;
    
    npict = cleanpixels.size() ; 
    _mparam.nplotpoints = 0 ;

    if(npict < 3) {
	
	_mparam.clear();
	
	return ;

    }

    totpmtmask.resize(_nwords) ;

    for(i=0; i < NX; i++) {
	for(j=0; j < NY; j++) {
            _ntrip[i][j] = 0 ;
            _xsum[i][j] = 0.0 ;
            _ysum[i][j] = 0.0 ;
            _rsum[i][j] = 0.0 ;

	    for (k=0; k<_nwords; k++) {
		_pmtmask[i][j][k] = 0;
	    }

	}
    } 
    n=0 ;
    
    // Loop through the N*(N-1)*(N-2)/(3*2*1) triplets, where N=npict
    // Note that cleanpixels[i] gives the pixel number of the ith
    // element of the list of picture pixels.

    int total_trips=0;

    for(i1=0; i1 < npict; i1++) {
	word1 = _single_pmt_mask[tubenum1=cleanpixels[i1]].word ;
	for(i2=i1+1; i2 < npict; i2++) {
            word2 = _single_pmt_mask[tubenum2=cleanpixels[i2]].word ;
            for(i3=i2+1; i3 < npict; i3++) {
		word3 = _single_pmt_mask[tubenum3=cleanpixels[i3]].word ;
		n++ ;

		getCenter( tubenum1,tubenum2,tubenum3, center );
    		ix = static_cast<int>(center.ix);
		iy = static_cast<int>(center.iy);
		assert( center.ix <= NX && center.iy <= NY );


		if((ix > 0) && (iy > 0)) {

		    _ntrip[ix][iy]++ ;

		    total_trips++;

		    // Accumulate averages at the appropriate bin
                
		    _rsum[ix][iy] += center.r ;
		    _xsum[ix][iy] += center.x0 ;
		    _ysum[ix][iy] += center.y0 ;
		
		    if(_plot_points_is_enabled) {     
			// Save some of the ring centers for display 
			if(_mparam.nplotpoints < MAX_PLOT_CENTERS) {
			    _mparam.plot[_mparam.nplotpoints].x = center.x0 ;
			    _mparam.plot[_mparam.nplotpoints].y = center.y0 ;
			    _mparam.nplotpoints++ ;
			} 
		    }
		
		    // Set the bits in pmtmask corresponding to the pmts
		    // indexed by i1,i2,i3
                
		    _pmtmask[ix][iy][word1] |= _single_pmt_mask[tubenum1].bit ;
		    _pmtmask[ix][iy][word2] |= _single_pmt_mask[tubenum2].bit ;
		    _pmtmask[ix][iy][word3] |= _single_pmt_mask[tubenum3].bit ;

		} // end if ix,iy in range 
            } // end for i3 
	} // end for i2 
    } // end for i1 

    // First find the index of the box-car smoothed ntrip array that
    // contains the maximum number of triplet ring centers.  Note: It
    // is not a mistake that the for loop runs as long as ix < NX-1 -
    // the -1 is needed because of the boxcar sum.

    nmax = 0 ;
    for(ix=0; ix < NX-1; ix++) {
	// find maximum smoothed ntrip 
	for(iy=0; iy < NY-1; iy++) { 
	    nbox = _ntrip[ix][iy]+_ntrip[ix+1][iy]+_ntrip[ix+1][iy+1]+
		_ntrip[ix][iy+1] ;
	    if(nbox > nmax) {
		nmax = nbox ;
		ixmax = ix ;
		iymax = iy ;
	    } 
	} 
    } 

    
    // TEST PLOT OF ntrip
//     Image2D ntripimg(NX,NY);
//     ntripimg.setCoordinateBox( -static_cast<double>(NX)/2.0*DX,
// 			       -static_cast<double>(NY)/2.0*DY,
// 			      static_cast<double>(NX)/2.0*DX,
// 			      static_cast<double>(NY)/2.0*DY  );
//     for(ix=0; ix < NX-1; ix++) {
// 	for(iy=0; iy < NY-1; iy++) { 
// 	    ntripimg.setPixel( ix,iy, _ntrip[ix][iy] );
// 	}
//     }

//     pm->setAxisBox( ntripimg );
//     pm->plot( ntripimg );
//     pm->plotAxes();
//     pm->flush();
//     pm->newPage();
    
    // calculate the centroid and width/length of the triplet-space image:
    
    double  momx=0;   
    double  momy=0 ;  
    double  momx2=0;  
    double  momy2=0 ; 
    double total=0;
    double twidth;

    for(ix=0; ix < NX-1; ix++) {
	for(iy=0; iy < NY-1; iy++) { 
	    momx    += _ntrip[ix][iy] * ix;
	    momy    += _ntrip[ix][iy] * iy;
	    momx2   += _ntrip[ix][iy] * ix*ix;
	    momy2   += _ntrip[ix][iy] * iy*iy;
	    total   += _ntrip[ix][iy];
	}
    }
    
    momx  /= total;   
    momy  /= total;   
    momx2 /= total;   
    momy2 /= total;   

    double vx2  = (momx2 - momx*momx);	     // <x^2> - <x>^2
    double vy2  = (momy2 - momy*momy);	     // <y^2> - <y>^2

    twidth = sqrt( (vx2 + vy2)/2.0  );
 
    // calculate the ring center
    if(nmax > 0) {
	dnmax = (double)nmax ;
	_mparam.ringcenter.x = (_xsum[ixmax][iymax]+_xsum[ixmax+1][iymax]+
				_xsum[ixmax+1][iymax+1]+
				_xsum[ixmax][iymax+1])/dnmax ;
	
	_mparam.ringcenter.y = (_ysum[ixmax][iymax]+_ysum[ixmax+1][iymax]+
				_ysum[ixmax+1][iymax+1]+
				_ysum[ixmax][iymax+1])/dnmax ;
	
	_mparam.radius = (_rsum[ixmax][iymax]+_rsum[ixmax+1][iymax]+
			  _rsum[ixmax+1][iymax+1]+_rsum[ixmax][iymax+1])/dnmax;
	
	for(i=0; i < _nwords; i++) {
	    totpmtmask[i]=_pmtmask[ixmax][iymax][i] |
		_pmtmask[ixmax+1][iymax][i] |
		_pmtmask[ixmax+1][iymax+1][i] |
		_pmtmask[ixmax][iymax+1][i] ;
	} 
	    
	phsum = 0.0 ;
	nph = 0 ;
	for(i=0; i < npict; i++) {

	    tubenum=cleanpixels[i] ;
	    wordnum=_single_pmt_mask[tubenum].word ;
		
	    // If a tube has contributed to the arc include it's
	    // pulseheight in the sum, and increment the total number
	    // of tubes contributing to arc, nph.

	    if((_single_pmt_mask[tubenum].bit & totpmtmask[wordnum]) > 0) {
		phsum += image[tubenum] ; 
		nph++ ;
	    } 


	} 

	// Calculate fraction of number of centers in the peak bin to
	// the total number of centers.
	if(npict > 3) 
// 	    _mparam.arcstrength = 6.0*dnmax/(double)(npict*(npict-1)*(npict-2)) ;
	    _mparam.arcstrength = dnmax / 4.0 ;
	else {
	    _mparam.arcstrength = 0.0 ;
	} 

	    
	if((_mparam.radius > 0.0) && (nph > 0)) {

	    // Modified 961009 by JB - Acutually while the light per
	    // unit pathlength is proportional to theta^2 or r^2, the
	    // pathlength increases as theta^-1 or r^-1, and the light
	    // per pixel is proportional to the light per azimuth
	    // angle divided by r, so actually the *gain should does
	    // not depend on radius to first order.

	    _mparam.gain = phsum/(double)(nph) ;
		 
	    // Added 961104 by JB - Calculate amount of light in an
	    // anulus about the ring center, approx. correct for
	    // fraction of area covered by pixels and for the linear
	    // dependence of the total amount of light on radius.

	    if (_calibration_is_enabled) {

		// CALIBRATION HERE


		// First calculate arc angle...

		double arcangle = 0.0 ;
		wsum = 0.0 ;
		xcs = 0.0 ;
		ycs = 0.0 ;
		     
		// For all pixels in the cleaned image, see if they
		// also contributed to the muon arc (i.e., were
		// contained in the boxcar averaged combined bitmask
		// for the grid position with the maximum number of
		// ring centers).  For all such pixels that
		// contributed to the muon arc, calculate the centroid
		// of the light distribution relative to the ring
		// center position.

		for(i=0; i<npict; i++) { 
		    tubenum = cleanpixels[i] ;
		    wordnum = _single_pmt_mask[tubenum].word ;

		    if((_single_pmt_mask[tubenum].bit & 
			totpmtmask[wordnum]) > 0) {
			
			dxa = _x[tubenum] - _mparam.ringcenter.x ;
			dya = _y[tubenum] - _mparam.ringcenter.y ;
			wt = image[tubenum] ;
			wsum += wt ;
			xcs += dxa*wt ;
			ycs += dya*wt ;
		    } 

		} 
		     
		if (wsum > 0.0) {
			 
		    // Calculate angle of vector from ring center to
		    // centroid of light distribution in ring
			 
		    xcs /= wsum ;   
		    ycs /= wsum ;   

 		    // wphi = M_PI+atan2(ycs,xcs) ;
		    
		    wphi = atan2(ycs,xcs) ;


// 		    cout << "DEBUG: muon: xcs="<<xcs
// 			 << " ycs="<< ycs << " wphi="<<wphi*180.0/M_PI
// 			 << endl;

		    // Next calculate the length of the vector from
		    // the centroid to the ring center, and divide by
		    // the ring radius.  Define a new quantity muskew
		    // to be the ratio of the two.  Note that muskew
		    // is small for a complete ring.  This may or may
		    // not be useful
			 
		    _mparam.muskew = sqrt(xcs*xcs+ycs*ycs)/_mparam.radius ;
			 
		    // In the spirit of the Hillas moment analysis,
		    // calculate the arclength from the second moment
		    // of the angle phi distribution (the angle of a
		    // vector from the most probably ring center to
		    // each pixel contributing to the ring)
			 
		    wphisq = 0.0 ;
		    wt_sum = 0.0;
		    for(i=0; i<npict; i++) { 
			tubenum = cleanpixels[i] ;
			wordnum = _single_pmt_mask[tubenum].word ;
			if( (_single_pmt_mask[tubenum].bit 
			     & totpmtmask[wordnum]) > 0 ) {
				 
			    dxa = _x[tubenum] - _mparam.ringcenter.x ;
			    dya = _y[tubenum] - _mparam.ringcenter.y ;
			    // phi = M_PI+atan2(dya,dxa) ;
			    phi = atan2(dya,dxa) ;
			    wt = image[tubenum] ;
			    wt_sum += wt;
			    dphi = fabs(phi-wphi) ;
			    dphi = GSL_MIN(dphi,fabs(dphi-2*M_PI)) ;
			    // dphi=gsl_sf_angle_restrict_pos(fabs(phi-wphi));
 			    // dphi = GSL_MIN_DBL(dphi,M_PI-dphi) ; 
			    
			    wphisq += (wt*dphi*dphi) ;
				 
			} 
		    }  // end for all pixels in the arc 
			 
		    arcangle = sqrt( wphisq/wt_sum ) ;
		    _mparam.philo = wphi - arcangle ;
		    _mparam.phihi = wphi + arcangle ;
		    _mparam.phimid = wphi;
		    _mparam.xcs = xcs + _mparam.ringcenter.x;
		    _mparam.ycs = ycs + _mparam.ringcenter.y;
		    
			 
		}
		else {  // if  wsum <= 0.0 
		    wphi = 0.0 ;
		    _mparam.philo = 0.0 ;
		    _mparam.phihi = 0.0 ;
		    _mparam.phimid =0.0 ;
		    arcangle = 1.0e20 ;        /* CHECK THIS */
		} 
		     

		// calculate size
		
		size_clean = 0.0;
		for (i=0; i<npict; i++) {
		    size_clean += image[cleanpixels[i]];
		}


		// =======================================================
		// Now loop through ALL pixels, calculating the total
		// amount of light in an annulus around the ring and
		// attempt to correct for the fraction of light lost
		// at the edge of the camera.  Also, we want to later
		// determine the "smoothness" of the angular
		// distribution of light for muon selection - for that
		// we need to generate a "phi" histogram here. We
		// define two annulli: FINE and COARSE, since we want
		// a more restrictuve annulus for the "smoothness".
		     
		npix = 0 ;
		narc = 0 ;
		stot = 0.0 ;
		sarc = 0.0 ;
		size_raw = 0.0;
		phihist.reset(); 


		double r2_coarse = _mparam.radius + DRAD_COARSE ;
		double r1_coarse = _mparam.radius - DRAD_COARSE ;
		if (r1_coarse < 0.0) 
		    r1_coarse = 0.0 ;
		double r1sq_coarse = r1_coarse*r1_coarse ;
		double r2sq_coarse = r2_coarse*r2_coarse ;

		double r2_fine = _mparam.radius + DRAD_FINE ;
		double r1_fine = _mparam.radius - DRAD_FINE ;
		if (r1_fine < 0.0) 
		    r1_fine = 0.0 ;
		double r1sq_fine = r1_fine*r1_fine ;
		double r2sq_fine = r2_fine*r2_fine ;

		double aanulus_coarse = M_PI*(r2sq_coarse-r1sq_coarse) ;
		double rpsq;

		for(i=0; i<getNumPixels(); i++) {
		    size_raw += image[i];
		    dxa = (_x[i] - _mparam.ringcenter.x) ;
		    dya = (_y[i] - _mparam.ringcenter.y) ;
		    rpsq = dxa*dxa + dya*dya ;

		    // If radius of point falls within the COARSE annulus:

		    if ((rpsq > r1sq_coarse)&&(rpsq < r2sq_coarse)) {

			phi = atan2(dya,dxa) ;
			dphi = fabs(phi-wphi) ;
			dphi = GSL_MIN_DBL(dphi,fabs(dphi-2*M_PI)) ;

			
			// if radius of point falls within the FINE
			// radius, use it for the "smoothness"
			// calculation:

			if ((rpsq > r1sq_fine)&&(rpsq < r2sq_fine)) {
			    
			    // accumulate the phi position,
			    // restricting angles to 0..2*PI: 
			    phihist.accumulate( fmod(phi+M_PI, 2.0*M_PI),
						image[i] );
			}
			     
			// Sum up the total amount of signal falling
			// in the annulus and with an angular distance
			// less than the arclength into one sum called
			// sarc
			     
			if (fabs(dphi) < arcangle) {
			    narc++ ; 
			    sarc += image[i] ;
			} 
			     
			// Sum up all of the signal falling in the
			// annulus into another sum called stot
		  
			npix++ ;
			stot += image[i] ;
		    } // end if rpsq 
		} // end for i
 

		// =======================================================
		// calculate the "smoothness" of the phi histogram for
		// muon selection.  Calculate the average difference
		// between bins (like derivative) and the variance of
		// that quantity.


// 		phihist.save("phihist"); // for DEBUG purposes and plotting

		double delta=0;
		double deltasum=0;
		double deltasum2=0;
		int nphi =0;
		for (i=0; i<phihist.numBins()-1; i++) {

		    delta = fabs(phihist[i+1] - phihist[i]);
		    deltasum += delta;
		    deltasum2 += delta*delta;
		    nphi++;

		}

		double avg;
		if (phihist.numBins()>0) 
		    avg = phihist.sum()/(double)phihist.numBins();
		else
		    avg = 1;
		
		_mparam.smoothness = deltasum/(double)nphi/avg;
		_mparam.smoothness_var = sqrt(deltasum2/(double)nphi)/avg;


		// =======================================================
		// Calculate the radial spread in the image about the
		// ring center

		double rp;
		double sumsigr2=0;
		double totalsignal;

		vector<int>::iterator pix;

		for (pix=cleanpixels.begin();pix!=cleanpixels.end();pix++) {

		    rp += sqrt(pow(_x[*pix] - _mparam.ringcenter.x ,2) +
			pow(_y[*pix] - _mparam.ringcenter.y ,2));
		    
		    if (rp>0) {

			totalsignal += image[*pix];
			sumsigr2 += pow( rp - _mparam.radius,2  )*image[*pix];

		    }
		    
		}

		_mparam.rspread = sqrt(sumsigr2/totalsignal);


		// =======================================================
		// define the muonness: 
		// eventually, this should include "smoothness" factor

		if (_mparam.radius > SATURATION_RADIUS-0.6 
		    && _mparam.radius < SATURATION_RADIUS+0.3) {
		    _mparam.muonness = (stot>0)? stot/size_raw:1e10;
		}
		else {
		    _mparam.muonness = 0.0;
		}


		// =======================================================
		// Calculate the area covered by all pixels with
		// coordinates inside of the annulus, compared with
		// the total area of the annulus.  The ratio
		// (ringfrac1) determines how much of the ring image
		// fell off of the camera.  The arclength (arclen) is
		// also calculated here
		
		if (npix > 0) {

		    if(aanulus_coarse > 1.0e-20) {
			ringfrac1 = (apix * (double)(npix))/aanulus_coarse ;
			if(ringfrac1 > 1.0e-20) {
			    _mparam.mugain = (1.0/36.0)*stot *
				(SATURATION_RADIUS/_mparam.radius)/ringfrac1 ;
			}
			_mparam.ringfrac = ringfrac1;
		    }
			

		    _mparam.arclen = arcangle * SATURATION_RADIUS;   // temporary
								     // until
								     // ring
								     // fit
								     // works
								     // better
// 		    _mparam.arclen = arcangle * _mparam.radius;

		    // NOTE: soal is defined as size over arc length.
		    
		    if(_mparam.arclen > 1.0e-20) {
			
			// ringfrac2 is the fraction of the area in the
			// arc from philo to phihi which falls on the
			// camera.  Note that since dN/dphi scales
			// like the radius, a correction for ring
			// radius is also required.  

			// JB 030617 - TODO: CHECK THIS - I just argued
			// earlier that a correction for radius was
			// not necessary.  Which is correct?

			// KPK: TODO: ringfrac2 is calculated but never used!

			ringfrac2 = (apix * (double)(narc)) / 
			    (GSL_MAX_DBL(1.e-20,aanulus_coarse)*
			     (_mparam.phihi-_mparam.philo)
			     /(2*M_PI));
			_mparam.soal = sarc/_mparam.arclen ;

		    } // end if arclen
		} // end if npix 
		else {
		    _mparam.mugain = 0.0 ;
		} 

	    } // end if fcalibrate 
	} // end if radius 
	else {
	    _mparam.gain = 0.0 ;
	    _mparam.mugain = 0.0 ;
	} 
	     
    } /* end if */

    // As always, set values to some overflow which flag an ivalid
    // event
	
    else {  
	_mparam.invalid = 1;
	_mparam.ringcenter.x = 1000.0 ;
	_mparam.ringcenter.y = 1000.0 ;
	_mparam.radius = 1000.0 ;
	_mparam.gain = 1.0e10 ;
	_mparam.mugain = 1.0e10 ;
	_mparam.arcstrength = 0.0 ;
    } /* end else */

       

}


ostream &operator<<( ostream &stream,const MuonParameterization &p) {

    stream << setw(10) << p.radius <<" "
	   << setw(10) << p.ringcenter.x <<" "
	   << setw(10) << p.ringcenter.y <<" "
	   << setw(10) << p.arcstrength <<" "
	   << setw(10) << p.gain <<" "
	   << setw(10) << p.mugain <<" "
	   << setw(10) << p.muskew <<" "
	   << setw(10) << p.ringfrac <<" "
	   << setw(10) << p.soal <<" "
	   << setw(10) << p.arclen<<" "
	   << setw(10) << p.muonness<<" "
	   << setw(10) << p.smoothness<<" "
	   << setw(10) << p.smoothness_var <<" "
	   << setw(10) << p.rspread <<" "
	;

    return stream;

}
