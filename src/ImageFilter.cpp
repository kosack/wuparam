/**
 *
 * C++ version of Armand Atoyan's HiRes.f program for filtering images
 * on a hexagonal grid.  I converted this code from a fortran program,
 * so it's pretty ugly and I've left the array indexing starting at 1
 * to avoid errors. This means, all arrays must be dimensioned with an
 * extra element, with the 0th element ignored.  The main functions
 * hcut_init() and hcut() (and the wrappers filter_init() and
 * filter_image()) assume that the arrays passed to them are c-indexed
 * (from 0).
 *
 * This code is for academic use only. 
 *
 * FOR REFERENCE, SEE:
 *
 * Atoyan, A., Patera, J., Sahakian, V., and Akhperjanian, A. 2005,
 * Astroparticle Physics, 23, 79
 */

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <valarray>
#include <vector>
#include <fstream>

#include "ImageFilter.h"

using namespace std;


void hcut( double *gpix_c );

void hcut_init( int npix, double *xpix_c,  double *ypix_c, 
		int msub, double afiltr, int nc, double spas);

double dinv(int n, int k, int l );
double f0su3nul(double xin, double yin, int n, double af);
void dofltr0(double s, int n0, double afil);
double f0pass(double xin, double yin, int n, double afil );
void doakl0(int nn0, double tfij[][201], int &k0re, int &k0im);
void derivs(int n,int kp,double width, double dpi);
void dohires(int n0, int msub, double afil, int &maxf);

void filter_test(void);

inline int max0(int a,int b) ;
inline int min0(int a,int b);


//
// Constants
//

const int MAX_PIX = 4000;
const int MAX_HIRES_PIX = 30000;
const int MAX_COEFFS=61;
const int MAX_DERIV=101;


//
// global variables; Common blocks have been converted to global
// anonymous structs with the same name as the original block
//

struct { 
    double re[MAX_COEFFS][MAX_COEFFS];
    double im[MAX_COEFFS][MAX_COEFFS];
    double enumb[MAX_COEFFS][MAX_COEFFS]; 
} Ak0;

struct { 
    double re[MAX_COEFFS][MAX_COEFFS];
    double im[MAX_COEFFS][MAX_COEFFS];
} Akf0;

struct{ 
    double alim;
    double blim;
} Ab;


struct {
    double xpact[2000]; 
    double ypact[2000];
    int ipact[2000]; 
    int jpact[2000];
} Pix;


struct { 
    double fij0[201][201];
    double g2hkp[MAX_PIX] ;
} A00;


struct {
    double x0[MAX_DERIV][MAX_DERIV];
    double xp[MAX_DERIV][MAX_DERIV];
    double xm[MAX_DERIV][MAX_DERIV];   
} D1H;


struct {
    double xx0[MAX_DERIV][MAX_DERIV]; 
    double xxp[MAX_DERIV][MAX_DERIV]; 
    double xxm[MAX_DERIV][MAX_DERIV];
} D2H;

struct { 
    double xhir[MAX_HIRES_PIX]; 
    double yhir[MAX_HIRES_PIX]; 
    int kphir;
} HiR;
		 

struct  {
    double f[MAX_HIRES_PIX];
    double x[MAX_HIRES_PIX]; 
    double y[MAX_HIRES_PIX];
    int k;
} Foto;


double Co[2000], Si[2000];
int Npix=0;
int N0=0, N1=0, N23=0, N13=0;
double Delx=0, Dely=0;
double *Xpix,*Ypix;
double Ypix_shift[MAX_PIX];
double Yshift;
double Dpix;
int Maxf;



namespace FilterParam {
    static int nc;
    static double afiltr;
    static int msub;
    static double spas;
} ;

namespace FilterInit {
    double xcoords[MAX_PIX];
    double ycoords[MAX_PIX];
    double im[MAX_PIX];
};


// int main(int argc, char** argv) {
//     filter_test();
// }


/**
 * Test routine
 */
void 
filter_test() {

    try {

	// ===========================================================
	// load the camera coordinates: 
	// TODO: use the Camera object eventually

	
	// these are numbered from 0, and later via a pointer trick
	// turned into fotran-style arrays:
	valarray<double> coordx(MAX_PIX); 
	valarray<double> coordy(MAX_PIX);
	double x,y;
	int i,npix;


	// ===========================================================
	// Load the camera coordinates:

	i=0;
	ifstream camxfile("coordx.txt");
	while (camxfile && camxfile >> x){
	    coordx[i++] = x;
	}
	camxfile.close();

	npix = i;

	i=0;
	ifstream camyfile("coordy.txt");
	while (camyfile && camyfile >> y ){
	    coordy[i++]= y;
	}
	camyfile.close();

      	if (npix != i) {
	    cout << "Error: size of coordx != size of coordy!" << endl;
	    exit(1);
	}
	cout << "Loaded "<<npix<<" camera coordinates"<<endl;
	
	// ===========================================================
	// Load the image:


	valarray<double> image(npix);

	i=0;
	ifstream imagefile("image.txt");
	while ( imagefile && imagefile >> x){
	    image[i++]=x;
	}
	imagefile.close();


	ofstream outfile("debug_image.txt");
	for (int i=0; i<npix; i++) {
	    outfile << setw(16) << coordx[i]
		    << setw(16) << coordy[i]
		    << setw(16) << image[i]
		    << endl;
	}
	outfile.close();

	cout << "Loaded "<<i<< " image points "<< endl;

	// ===========================================================
	// Filter

	cout << "Filtering... " << endl;

	valarray<double> hxc,hyc;
	valarray<double> fimage(image.size());

	filter_init( npix, coordx, coordy, hxc,hyc, 3,0.25,2,0.0 );
       	filter_image( image,fimage );

	cout << "Done."<<endl;


	// ===========================================================
	// output the smoothed image

	// at this point, the Foto.{x,y,f} arrays have been generated,
	// providing the x/y coordinate of each point in the image and
	// the signal.

	cout << "DEBUG: writing out "<<Foto.k<<" points"<<endl;
	ofstream smfile("debug_smoothed.txt");
	for (int i=1; i<=Foto.k; i++) {
	    smfile  << setw(16) << Foto.x[i]
		    << setw(16) << Foto.y[i]
		    << setw(16) << Foto.f[i]
		    << endl;
	}
	smfile.close();

	cout << "DEBUG: "  <<endl
	     << " size of: hxc="<< hxc.size() <<endl
	     << " size of: hyc="<< hyc.size() <<endl
	     << " size of: fimage="<< fimage.size() <<endl;
	cout.flush();

	cout << "DEBUG: "  <<endl
	     << " max: hxc="<< hxc.max() <<endl
	     << " max: hyc="<< hyc.max() <<endl
	     << " max: fimage="<< fimage.max() <<endl;
	cout <<"DEBUG: xc,yc from 0: "<< endl;
	for (i=0; i<10; i++) {
	    cout <<i << ": " << hxc[i] << " , " << hyc[i]<<endl;
	}
	
    }
    catch (runtime_error &e) {
	cout << "RuntimeError: "<<e.what()<<endl;
    }
    catch(...) {
	cout << "Caught unknown exception." << endl;
    }

}


/**
 *  Prepare the image filter and generate hires camera:
 * 
 *  \param xc array of x camera coords
 *  \param yc array of x camera coords
 *  \param newxc output array of x camera coords, may be larger than original
 *  \param newyc output array of x camera coords
 */
void filter_init( int npix, valarray<double> &xc, valarray<double> &yc,
		  valarray<double> &newxc, valarray<double> &newyc,
		  int msub, double afiltr, int nc, double spas ) {
    
    
    // first copy the valarrays into regular arrays:
    
    for( int i=0; i< npix; i++) {
	FilterInit::xcoords[i] = xc[i];
	FilterInit::ycoords[i] = yc[i];
    }

    // initialize:
    
    hcut_init( npix, FilterInit::xcoords, FilterInit::ycoords,
	       msub, afiltr, nc, spas );


    // now put the hires-camera into the "new" variables:
    
    newxc.resize( HiR.kphir );
    newyc.resize( HiR.kphir );

    for (int i=0; i<HiR.kphir; i++) {
	
	newxc[i] = HiR.xhir[i+1]; // note HiR is fortran-indexed
	newyc[i] = HiR.yhir[i+1] - Yshift;

    }

    cout << "DEBUG: image filter initialized" <<endl;

};

/**
 * Wrapper routine for WUparam integration
 *
 * \param im array of image values
 * \param newim output array of filtered image values
 * \returns maximum value of the filtered image 
 *
 */
double filter_image(  valarray<double> &im, valarray<double> &newim ) {
        
    // put the image into a standard array:

    for (int i=0; i<im.size(); i++) {
	FilterInit::im[i] = im[i];
    }
    
    // filter the image:

    hcut( FilterInit::im );
    
    if (newim.size() != Foto.k) newim.resize(Foto.k);
    
    // Copy the new image into newim 

    // TODO: eventually should just convert everything to valarrays so
    // this doesn't need to happen.  For now, though, the fastest way
    // to put this into a valarray is using a copy constructor.  i
    // measured this to be approximately 200% faster than using a
    // for-loop.

    newim = valarray<double>(&Foto.f[1],Foto.k);

    return Foto.f[Maxf];

}



void hcut_init( int npix, double *xpix_c,  double *ypix_c, 
		int msub, double afiltr, int nc, double spas) {

    cout << "DEBUG: hcut: npix="<<npix<<endl;;

    FilterParam::afiltr = afiltr;
    FilterParam::msub= msub;
    FilterParam::nc = nc;
    FilterParam::spas = spas; //    Spas - the "amplitude" in the filtr

    Npix = npix;
    Xpix = xpix_c-1;	    // pointer hack to allow fortran numbering
    Ypix = ypix_c-1;	    // pointer hack to allow fortran numbering

    const double sq32=sqrt(3.)/2.;
    const double sq3=sqrt(3.);
    
    int nring=11;  
    N0=3*(nring+1);     // number of intervals along each
			    // side of the ffsu3
    N1=N0+1;
    int npix_check=N0*nring+1;

    if (npix != npix_check) {
	cout << "npix="<<npix<<" !=  npix_check="<<npix_check;
	//	throw runtime_error("hcut: wrong number of pixels in camera");
    }

    if (npix > MAX_PIX) {
	cout << "Too many pixels (npix > "<< MAX_PIX<<")"<<endl;
	throw runtime_error("hcut: number of pixels out of bounds");
    }





    // pixel size in degrees (assumes pixels 1 is adjacent to pixel 2)
	
    Dpix = sqrt(pow(Xpix[2]-Xpix[1],2)+pow(Ypix[2]-Ypix[1],2));
    cout << "DEBUG: Dpix="<<Dpix<<endl;


    Ab.blim= (float)N0*Dpix/2.;
    Ab.alim=-Ab.blim;

    double ba=(Ab.blim-Ab.alim);
    Yshift=ba/sq3;

    Delx=Dpix;
    Dely=Delx*sq32;
    N23=2*N0/3;
    N13=N0/3;
    
    int i,j,i1,j1;


    // shift the camera so that 0,0 corresponds with the bottom of the triangle


    for (i=1; i<=Npix; i++) {
	Ypix_shift[i] = Ypix[i] + Yshift;
    }

    int n0hir=N0*FilterParam::msub;
    int n1hir=n0hir+1; //  the number of POINTS along the side of
			  //  the triangle

    double dxhir=Delx/(double)FilterParam::msub;
    double dyhir=Dely/(double)FilterParam::msub;
    int m13=N13*FilterParam::msub;
    int m23=N23*FilterParam::msub;

//    HiR.Kphir - number of pixels in the possibly active region, - hexagonal camera
//    here the coordinates are prepared for further computation

    int jmax;

    for (i1=1; i1<=n1hir; i1++) {
	i=i1-1;
	jmax=n1hir-i;

	for( j1=1; j1<=jmax; j1++) {
	    j=j1-1;
	    if(i+j < m13) continue;
	    if(i > m23+1) continue;
	    if(j > m23+1) continue;
	    HiR.kphir++;
	    HiR.xhir[HiR.kphir] =0.5*dxhir*(float)(i-j);
	    HiR.yhir[HiR.kphir] = dyhir*(float)(i+j);
	}
    }


    int n61=6*N0+1;
    double del=2.0*M_PI/(float)(6*N0);
    double x0=0.0;
    Co[1]=1.0;
    Si[1]=0.0;
    for(i=2; i<=n61; i++) {
	x0=del*(float)(i-1);
	Co[i]=cos(x0);
	Si[i]=sin(x0);
	if(fabs(Co[i])<3.e-7) Co[i]=0.;
	if(fabs(Si[i])<3.e-7) Si[i]=0. ;
    }

}



	
/**
 * Initializes all the variables, and processes the data
 *
 * \todo: make this operate on a single image, not an array of images.
 * \todo: parameters should be arguments
 * \todo: pass in camera coordinates.
 *
 * \param xpix - x camera pixel coordinates
 * \param ypix - y camera pixel coordinates
 *
 * \param msub - the number of subdivisions of the initial interval,
 * the density increases by msub^2
 * \param afiltr - filter parameter
 * \param nc - number of iterations for derivative calculation
 *
 */ 
void hcut(  double *gpix_c ) {
    
    double *gpix = gpix_c-1; // hack to allow fortran numbering

    // done with prep
    //--------------------


    //    int imn=0;
    //    int nsig;
    int kpact=0;
    double xij,yij;
    double xr2,yr2,rr;
    
    int i,i1,j1,ip;
    
    for (i=1;i<=Npix;i++) {
	A00.g2hkp[i]=0.;
    }
        
    //    filling/cleaning the arrays, and redistribution of one-D image
    //    data, GPIX(I) into 2D- data on the grid, ENUMB(I,J)
    
    int j,jmax;

    for (i1=1; i1<=N1; i1++) {
	i=i1-1;
	jmax=N1-i;
        
	for (j1=1;j1<=N1;j1++) {
	    j=j1-1;
// 	    xij=0.5*Delx*(float)(i-j);
// 	    yij=Dely*(float)(i+j);
	    
	    
	    A00.fij0[i1][j1]=0.;
	    Ak0.enumb[i1][j1]=0.;
	    D1H.x0[i1][j1]=0.;
	    D1H.xp[i1][j1]=0.;
	    D1H.xm[i1][j1]=0.;
	    D2H.xx0[i1][j1]=0.;
	    D2H.xxp[i1][j1]=0.;
	    D2H.xxm[i1][j1]=0.;
	    
	    if(i >= N23 || j >= N23) continue;
	    if((i+j) < N13) continue;

	    xij=0.5*Delx*(float)(i-j);
	    yij=Dely*(float)(i+j);

// 	    cout << "DEBUG: deplx="<<Delx
// 		 << " Dely="<<setw(6) <<Dely 
// 		 << " Dpix="<<setw(6) <<Dpix 
// 		 << " xij="<<setw(6) <<xij
// 		 << " yij="<<setw(6) <<yij
// 		 <<" 0.2*Dpix="<<setw(6) <<0.2*Dpix
// 		 << " npix= "<< npix
// 		 << endl;

	    for(ip=1;ip<=Npix;ip++) {
		xr2=pow( (Xpix[ip]-xij),2) ;
		yr2=pow( (Ypix_shift[ip]-yij),2);
		rr=sqrt(xr2+yr2); 

// 		cout << "DEBUG: xpix= "<< setw(6) << xpix[ip]
// 		     << " Ypix="<< setw(6) << Ypix[ip] 
// 		     << " rr="<< setw(6) << rr
// 		     << " gpix="<< setw(6) << gpix[ip]
// 		     << endl;


		if(rr >= 0.2*Dpix) continue; 
		Ak0.enumb[i1][j1]=gpix[ip];
		
		if(gpix[ip] < 0.1) break; 

		kpact++;
		if(kpact >= 2000) 
		    cout << "Large  Kpact="<<kpact<<endl;
				Pix.xpact[kpact]=Xpix[ip];
		Pix.ypact[kpact]=Ypix_shift[ip];
		Pix.ipact[kpact]=i1;
		Pix.jpact[kpact]=j1;
	    } //75

	} // 73
    } //74
    


//    If(Nsig.GT.100) Print*,'I EAS=',ImN,'  KpAct=',KPACT,Nsig


//    assigning initial values to Fij0: Fij0=0 on the boundary;
         
    double cf0=5./36.;
    double s6;
    int ij,k;

    for (k=1; k<=kpact; k++) {
	i1 = Pix.ipact[k];
	j1 = Pix.jpact[k];
	ij=i1+j1-2;
	if (i1 >= 1 || j1==1 || ij == N0) continue;
	s6 = Ak0.enumb[i1+1][j1-1]+Ak0.enumb[i1-1][j1+1]+Ak0.enumb[i1][j1+1];
	s6 = s6+Ak0.enumb[i1][j1-1]+Ak0.enumb[i1+1][j1]+Ak0.enumb[i1-1][j1];
	A00.g2hkp[k]=s6/6. - Ak0.enumb[i1][j1];
	A00.fij0[i1][j1] = Ak0.enumb[i1][j1] - cf0*A00.g2hkp[k];
    }


//    Calculations of the coefficients Akl 

//    double cf2=5./54.;
    double cf3=5./36.;       
    //    double f0max=0;
    int ic;
    int k0re, k0im;
    int ibig, jbig;
    double hwd=0.25;   

    for ( ic= 1; ic<= FilterParam::nc; ic++ ){ // iterate nc times for accuracy

	doakl0(N0,A00.fij0,k0re,k0im);
	derivs(N0,kpact,hwd,Dpix);

	for( k=1; k<=kpact; k++) {
	    ibig=Pix.ipact[k];
	    jbig=Pix.jpact[k];
	    i1=2*ibig-1;
	    j1=2*jbig-1;
	    A00.fij0[ibig][jbig]=Ak0.enumb[ibig][jbig]-cf3*A00.g2hkp[k];
	}

    }

//    calculation of the transform coefficients A_{kl}

    doakl0(N0,A00.fij0,k0re,k0im)   ;

//    filtering-smoothing the coefficients that can be later used in the
//    program F0Pass(,,) for computation of HiR image values, below it
//    is computed with F0SU3NUL , corresponding to sharp cuttof
//    filtering only.

    dofltr0(FilterParam::spas,N0, FilterParam::afiltr);

//    This block sends the initial grid parameters and the filterinf
//    parameter Afiltr, and yealds the high resolution image values 1D
//    array, Ffoto(), and the corresponding coordinates of active pixels
//    only

//    Maxf - is the number of the pixel corresponding to the maximum of
//    the new image; the maximum of the new image corresponds to
//    Ffoto(Maxf); this falue can be fuirther used for %-cuts

    dohires(N0,FilterParam::msub,FilterParam::afiltr, Maxf);


    // Now, we need to shift the y-coordinates back to the original position:

    for (i=1; i<=Foto.k; i++) {
	Foto.y[i] -= Yshift;
    }


//    =============================================


} 




double dinv(int n, int k, int l ){ 

    //    cout << "DEBUG: in dinv() "<< endl;

    if (n == 0 ) return 0;
    if ( (k+l) > n) return 0;
    
    int mk=k/n;
    int ml=l/n;
    int mkl=(k+l)/n;
    int k0=k-n*mk;
    int l0=l-ml*n;
    int kl0=k+l-mkl*n;
    double cl=0.;
    double ck=0.;
    double ckl=0.;
   
    if(l0==0) cl=1.;
    if(k0==0) ck=1.;
    if(kl0==0) ckl=1.;
    return 1./((1.+ckl)*(1.+ck+cl));
}


double f0su3nul(double xin, double yin, int n, double af) {
//      COMMON/AB/ Alim, Blim
//      COMMON/AK0/ AKLRE0(61,61),AKLIM0(61,61),ENUMB(61,61)
    double coy[501], cox[501], six[501];


    double p2=M_PI*2.0;
    int nn=n+1;
    double rn=(double)n;
    int ncut= (int)(rn * (1.- af));
    int nc1=ncut+1;
    double del=Ab.blim-Ab.alim;
    double cen=(Ab.blim+Ab.alim)/2.;
    double x0=(xin-cen)/del;
    double y0=yin/del;
    double sq3=sqrt(3.);
    double  ysq3=y0/sq3;

    if(fabs(x0) >= ysq3) return 0;
    if(y0>=sq3/2.) return 0;

    int ir=(int)(rn*(x0+ysq3));
    int jr=(int)(rn*(ysq3-x0));
    int i1=max0(ir+1,2);
    int j1=max0(jr+1,2);
    double x,y;
    int i,j;
    double s14;
    int m1,m2,n2;
    double xi,yi;

//     int imn2=max0(ir-2,1);
//     int jmn2=max0(jr-2,1);
//     int ipl2=min0(i1+1,nn);
//     int jpl2=min0(j1+1,nn);
//     int ipl3=min0(ir+3,nn);
//     int jpl3=min0(jr+3,nn);
//     int imn1=max0(ir-1,1);
//     int jmn1=max0(jr-1,1);
//     s14 = Ak0.enumb[ipl2][jmn2] + 
// 	Ak0.enumb[imn2][jpl2] +
// 	Ak0.enumb[imn1][jpl3];
//     s14+=Ak0.enumb[ipl3][jmn1];
    
    // enable or disable this?:
    //    s14=-Ak0.enumb[i1-1][j1-1]-Ak0.enumb[i1+2][j1+2];


//     for (m1=1; m1<=4; m1++) {
// 	for( m2=1; m2<=4; m2++) {
//             s14 += Ak0.enumb[i1-2+m1][j1-2+m2];
// 	}
//     }
//     if(s14<1.e-4) return 0;

    x=p2*x0/3.;
    y=p2*y0/sqrt(3.0);
    n2=nn+n;
    xi=0.;
    cox[1]=1.000;
    six[1]=0.0;
    for(i=2; i<=n2; i++) {
	xi=xi+x;
	cox[i]=cos(xi);
	six[i]=sin(xi);
    }
    yi=0.;
    coy[1]=2.0;

    for (j=2; j<=nn; j++) {
	yi=yi+y;
	coy[j]=2.*cos(yi);
    }

    double s0=0.;
    int nhalf=1+n/2;
    double fr;
    int k1,k0;

    for (k1=1;k1<=nhalf; k1++) {
	k0=k1-1;
	if(2*k0>ncut) continue;
	fr=coy[k0+k1]+2.*coy[k1]*cox[3*k0+1];
	s0=s0+fr*Ak0.re[k1][k1];
    }

    double s1=0.;
    int l0,l1,lm;
    int kml, k2l, kpl, l2k;
    double fi;

    for (k1=2; k1<=nc1; k1++) {
	k0=k1-1;
	lm=min0(nn-k0,k0)       ;
	
	for(l1=1;l1<=lm; l1++) {
            l0=l1-1;
            if((k0+l0)>ncut) break;
            kml=k1-l1+1;
            k2l=2*k0+l1;
            l2k=2*l0+k1;
            kpl=k0+l1;
            fr=coy[kpl]*cox[kml]+coy[l1]*cox[k2l]+coy[k1]*cox[l2k];
            fi=coy[kpl]*six[kml]-coy[l1]*six[k2l]+coy[k1]*six[l2k];    
            s1=Ak0.re[k1][l1]*fr - Ak0.im[k1][l1]*fi + s1;
	}
    }

    return s0 + 2.*s1;
}



void dofltr0(double s, int n0, double afil) {
//       COMMON/AKF0/ AFRE0(61,61),AFIM0(61,61)
//       COMMON/AK0/ AKLRE0(61,61),AKLIM0(61,61),AKL.ENUMB0[61,61)

    int n1=n0+1;
    int npas=(int)((float)(n0)*(1.-afil));
    int lpha=10;

    int i,i1,j1;
    int iharm, jmax;
    double r;

    for (i1=1;i1<=n1; i1++) {
	i=i1-1;
	jmax=n1-i;
        for(j1=1;j1<=jmax; j1++) {
            iharm=i+j1-1;
            r=1.-s*pow((float)(iharm)/(float)(npas), lpha);
	    if(r<0.) r=0.;
            Akf0.re[i1][j1]=Ak0.re[i1][j1]*r;
            Akf0.im[i1][j1]=Ak0.im[i1][j1]*r;
	}
    }

}


double f0pass(double xin, double yin, int n, double afil ) {
//       COMMON/AB/ Alim, Blim
//       COMMON/AKF0/ AFRE0(61,61),AFIM0(61,61)
//       COMMON/AK0/ AKLRE0(61,61),AKLIM0(61,61),AKL.ENUMB0[61,61)

//    cout << "DEBUG: in f0pass() "<< endl;

    double coy[501], cox[501], six[501];
    double p2=M_PI*2.0;
    int nn=n+1;
    int nmax=(int)((float)(n)*(1.-afil));
    int nc1=nmax+1;
    double del=Ab.blim-Ab.alim;
    double cen=(Ab.blim+Ab.alim)/2.;
    double     x0=(xin-cen)/del;
    double y0=yin/del;
    double sq3=sqrt(3.);
    double ysq3=y0/sq3;
    if(fabs(x0)>=ysq3)  return 0;
    if(y0>=sq3/2.) return 0;

    double rn=(double)(n);
    int ir=(int)(rn*(x0+ysq3));
    int jr=(int)(rn*(ysq3-x0));
    int i1=max0(ir+1,2);
    int j1=max0(jr+1,2);
    int imn2=max0(ir-2,1);
    int jmn2=max0(jr-2,1);
    int ipl2=min0(i1+1,nn);
    int jpl2=min0(j1+1,nn);
    int ipl3=min0(ir+3,nn);
    int jpl3=min0(jr+3,nn);
    int imn1=max0(ir-1,1);
    int jmn1=max0(jr-1,1);
    double s14=Ak0.enumb[ipl2][jmn2]
	+ Ak0.enumb[imn2][jpl2]
	+ Ak0.enumb[imn1][jpl3];

    s14 += Ak0.enumb[ipl3][jmn1]+s14;
//    S14=-Ak0.enumb[I1-1][J1-1]-Ak0.enumb[I1+2][J1+2]

    int m1,m2;

    for (m1=1;m1<=4;m1++) {
	for (m2=1;m2<=4;m2++) {
	    s14 += Ak0.enumb[i1-2+m1][j1-2+m2];
	}
    }
    if(s14<1.e-4) return 0;

//    IPJ=Int(RN*2.*Ysq3)-IR-JR
//    S3ang=Ak0.enumb0[I1][J1+1]+Ak0.enumb0[I1+1][J1]+Ak0.enumb0[I1+IPJ][J1+IPJ]
//    IF(S3ang<1.e-4) Return

    

    double x=p2*x0/3.;
    double y=p2*y0/1.7320508;
    int n2=nn+n;
    double xi=0.;
    int i;

    cox[1]=1.000;
    six[1]=0.0;
    for(i=2;i<=n2; i++) {
	xi=xi+x;
	cox[i]=cos(xi);
	six[i]=sin(xi);
    }
    
    double yi=0.;
    int j;

    coy[1]=2.0;
    for (j=2; j<nn; j++) {
	yi=yi+y;
	coy[j]=2.*cos(yi);
    }

    double s0=0.;
    int  nhalf=1+n/2;
    int k1,k0;
    double fr;

    for ( k1=1;k1<=nhalf; k1++) {
	k0=k1-1;
	if(2*k0>nmax) continue;
	fr=coy[k0+k1]+2.*coy[k1]*cox[3*k0+1];
	s0=s0+fr*Akf0.re[k1][k1];
    }

   double  s1=0.;
   int l1,lm,l0,kml,k2l,l2k,kpl;
   double fi;
   
   for( k1=2; k1<=nn ;k1++) {
       k0=k1-1;
       lm=min0(nc1-k0,k0) ;
       for(l1=1;l1<=lm;l1++) {
	   l0=l1-1;
	   if(l0+k0>nmax) break;
	   kml=k1-l1+1;
	   k2l=2*k0+l1;
	   l2k=2*l0+k1;
	   kpl=k0+l1;
	   fr=coy[kpl]*cox[kml]+coy[l1]*cox[k2l]+coy[k1]*cox[l2k];
	   fi=coy[kpl]*six[kml]-coy[l1]*six[k2l]+coy[k1]*six[l2k];
	   s1=Akf0.re[k1][l1]*fr - Akf0.im[k1][l1]*fi + s1;
       }
   }

   return s0+2.*s1;

}



void doakl0(int nn0, double tfij[][201], int &k0re, int &k0im) {

    double fdij[121][121];
    double tdkl[121][121] ;
//       COMMON/AK0/ AKLRE0(61,61),AKLIM0(61,61),AKL.ENUMB0[61,61)

    int n0=nn0;
    int n1=n0+1;
    k0re=0;
    k0im=0;

    if(n1>MAX_COEFFS) cout << "ERROR: Insufficient space"<<endl;
    int n6=6*n0 ;
    int nhalf=1+n0/2;

    int i1,j1,lm,i,j;

    for ( i1=1; i1<=n1; i1++) {
	i=i1-1;
	lm=n1-i;
	for (j1=1;j1<=lm;j1++) {
            j=j1-1;
	    tdkl[i1][j1]=dinv(n0,i,j);
            fdij[i1][j1]=2.*tfij[i1][j1]*tdkl[i1][j1];
	}
    }


    double cn=1./(float)(3*n0*n0);
    int k1, k, kipj,kplipj,k2limj,m1,m2,m3,m4,m5,m6;
    double sre, sim,fr;
    int ipj,imj,is;

    for (k1=1;k1<=nhalf;k1++){
	k=k1-1;
	sre=0.;
	sim=0.;
	for(i1=1; i1<=n1; i1++) {
            i=i1-1;
            lm=n1-i;
            for(j1=1;j1<=lm; j1++) {
		if(fabs(fdij[i1][j1])<1.e-7) continue;
		j=j1-1;
		ipj=i+j;
		imj=i-j;
		is=1;
		if(i<j) is=-1;
		kipj=3*k*ipj;
		kplipj=(kipj+kipj);
		k2limj=(2*k+k)*imj*is;
		m1=kplipj-n6*(kplipj/n6) +1;
		m2=kipj-n6*(kipj/n6) +1;
		m3=k2limj-n6*(k2limj/n6) +1;
		fr=Co[m1]+2.*Co[m2]*Co[m3];
		sre=sre+fr*fdij[i1][j1];
	    }
	}
	Ak0.re[k1][k1]=cn*sre*tdkl[k1][k1];
	Ak0.im[k1][k1]=0.0;
	if(fabs(Ak0.re[k1][k1])<1.e-6) k0re=k0re+1;
	k0im=k0im+1;
    }

    int lmax,l1,l,jm,kpl,k2l,l2k,kml,lipj;
    double fi;

    for( k1=2;k1<=n1;k1++ ) {
	k=k1-1;
	lmax=min0(n1-k,k);
	for ( l1=1;l1<=lmax;l1++) {
            l=l1-1;
            sre=0.;
            sim=0.;
	    for( i1=1;i1<=n1; i1++) {
		i=i1-1;
		jm=n1-i;
		for ( j1=1;j1<=jm; j1++) {
		    if(fabs(fdij[i1][j1])<1.e-7) continue;
		    j=j1-1;
		    ipj=i+j;
		    is=1;
		    if(i<j) is=-1;
		    imj=(i-j)*is;
		    kipj=3*k*ipj;
		    lipj=3*l*ipj;
		    kpl=kipj+lipj;
		    k2l=(2*k+l)*imj;
		    l2k=(2*l+k)*imj;
		    kml=k2l-l2k;
                  
		    m1=kpl-n6*(kpl/n6)+1;
		    m2=kipj-n6*(kipj/n6)+1;
		    m3=lipj-n6*(lipj/n6)+1;
		    m5=k2l-n6*(k2l/n6)+1;
		    m6=l2k-n6*(l2k/n6)+1;
		    m4=kml-n6*(kml/n6)+1;

		    fr=Co[m1]*Co[m4]+Co[m3]*Co[m5]+Co[m2]*Co[m6];
		    fi=(Co[m1]*Si[m4]-
			Co[m3]*Si[m5]+Co[m2]*Si[m6])*(float)(is);
		    sre=sre+fr*fdij[i1][j1];
		    sim=sim+fi*fdij[i1][j1];
		}
	    }
            Ak0.re[k1][l1]=cn*sre*tdkl[k1][l1];
            Ak0.im[k1][l1]=-cn*sim*tdkl[k1][l1];
            Ak0.re[l1][k1]=Ak0.re[k1][l1];
            Ak0.im[l1][k1]=-Ak0.im[k1][l1];
            if(fabs(Ak0.re[k1][l1])<1.e-6) k0re=k0re+2;
            if(fabs(Ak0.im[k1][l1])<1.e-6) k0im=k0im+2;
	}
    }

//     cout << "DEBUG: coefficients: ";
//     for ( i=1; i<=n1; i++) {
// 	cout << i << ":" << Ak0.re[i][1] << endl;
//     }

}
      
//    calculation of the coordinates of the triangles in the hexagon
//    devided into 6 K^2 triangles;    
//    


void derivs(int n,int kp,double width, double dpi) {
//       COMMON/PIX/ XpAct(2000),YpAct(2000),IpAct(2000),JpAct(2000)
//       COMMON/D2H/ D2HXX0(101,101),D2HXXP(101,101),D2HXXM(101,101)
//       COMMON/A00/  FIJ0(201,201),A00.G2hkp(MAX_PIX)
//       COMMON/D1H/ D1HX0(101,101),D1HXP(101,101),D1HXM(101,101)

	
    int n0=n;
    //int l0=2*n;
    double h0=dpi/2.;
    double dx=h0*width;
    double amp1=1./width/2.;
    double amp2= pow((1./width),2);
    double dy=sqrt(3.)*dx/2.;
    double c23=2./3.;
    int i,j,k;
    double x0,y0;
    double f1,f2,f3,f4,f5,f6,f02;
    
    for (k=1; k<=kp;k++) {
	i=Pix.ipact[k];
	j=Pix.jpact[k];
	x0=Pix.xpact[k];
	y0=Pix.ypact[k];
	f02 = 2. * f0su3nul(x0,y0,n0,0.);
	
	f1=f0su3nul(x0+dx,y0,n0,0.);
	f2=f0su3nul(x0+dx/2.,y0+dy,n0,0.);
	f3=f0su3nul(x0-dx/2.,y0+dy,n0,0.);
	f4=f0su3nul(x0-dx,y0,n0,0.);
	f5=f0su3nul(x0-dx/2.,y0-dy,n0,0.);
	f6=f0su3nul(x0+dx/2.,y0-dy,n0,0.);
	
	
	D1H.x0[i][j]=(f1-f4)*amp1;
	D1H.xp[i][j]=(f2-f5)*amp1;
	D1H.xm[i][j]=(f6-f3)*amp1;
	D2H.xx0[i][j]=amp2*(f1+f4-f02);
	D2H.xxp[i][j]=amp2*(f2+f5-f02);
	D2H.xxm[i][j]=amp2*(f3+f6-f02);
	A00.g2hkp[k]=c23*(D2H.xx0[i][j]+D2H.xxp[i][j]+D2H.xxm[i][j]);
    }
    
}



void dohires(int n0, int msub, double afil, int &Maxf) {
//       Common/HiR/ XHiR(30000),YHiR(30000),KpHiR
//       Common/Foto/ Ffoto(30000),Xfoto(30000),Yfoto(30000),Kfoto

    double  ffmax=0.;
    //    int n0hir=n0*msub;
    int k;
    double xi,yi,gn;

    Foto.k=0;
    Maxf=0;


    for (k=1; k<=HiR.kphir; k++) {
	xi=HiR.xhir[k];
	yi=HiR.yhir[k];
	gn=f0su3nul(xi,yi,n0,afil);
	// 	gn=f0pass(xi,yi,n0,afil);
//	if(gn <= 1.e-5) continue; // kpk removed this so we get the whole thing
       
	Foto.k++;
	Foto.f[Foto.k]=gn;
        Foto.x[Foto.k]=xi;
        Foto.y[Foto.k]=yi;

	if(gn>ffmax){
	    ffmax=gn;
	    Maxf=Foto.k;
	}
    }
    


}



/**
 * Helper function replacing the fortran intrisic max0
 */
int max0(int a,int b) {
    return (a>=b)? a:b;
}

/**
 * Helper function replacing the fortran intrisic min0
 */
int min0(int a,int b) {
    return (a<=b)?a:b;
}
