//
// wuplot - external plotting utility for plotting 2d data from WUParam, etc.
//
// Karl Kosack


#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <list>
#include <exception>
#include <stdexcept>
#include <unistd.h>
#include <signal.h>

#include "Types.h"
#include "PlotMaker.h"
#include "Image2D.h"
#include "Config.h"
#include "StarCatalog.h"
#include "Exceptions.h"

//#include "wuplot.h" 

using namespace std;

void parsePlot( string filename );
string getStringLiteral( vector<string> &v, int startelement );
void expect( int, vector<string> &v) ;
void readCurve( string filename, vector<double> &xpoints, 
		vector<double> &ypoints );

void prepareQuotedStrings( string &line );

class ParseException : public runtime_error {
public:
    ParseException(string what): runtime_error(what){;}
};


namespace WUPlot {

    bool xdisplay=true;
    string outputfilename;
    string outputtype="ps";
    int marker = 16;
    int smooth = 0;
    
    void showUsage();
    void getCommandLineOptions(int argc , char **argv);

};

const char SPACE_PLACEHOLDER = (char)254;


int
main( int argc , char **argv ) {
    
    cout << "===================================="
	 << "====================================" << endl
	 << "WUPlot - data plotter v"<< VERSION
	 << " <kosack@hbar.wustl.edu> " << endl
	 << "Washington University Physics Department, St. Louis M0"
	 << endl;
    cout << "===================================="
	 << "===================================="
	 << endl << endl;
    
    
    WUPlot::outputfilename = "output.eps";

    WUPlot::getCommandLineOptions( argc, argv );

    if (argc-optind < 1) {
	WUPlot::showUsage();
	exit(1);
    }

    if (WUPlot::outputfilename == argv[optind]) {
	cout << "You can't output to the same file as the input!" << endl;
	exit(1);
    }

    try {
    
	parsePlot( argv[optind] );

    }
    catch( AnalysisException &e) {
	cout << "AnalysisException: "<< e.what() << endl;
	exit(1);
    }
    catch (exception &e) {
	cout << "CRITICAL: "<< e.what() << endl;
	exit(1);
    }
    catch (...) {
	cout << "HUH?: caught an unknown exception!" << endl;
	exit(1);
    }

    if (WUPlot::xdisplay)
	cout <<"Tip: Use the -o <filename> option to "
	     << "output to an EPS file (by default) or other format"<<endl;
    else 
	cout << "Wrote "<<WUPlot::outputtype
	     <<" output to '"<<WUPlot::outputfilename<<"'"<<endl;
    
    if (WUPlot::outputtype=="ps") {
	cout << "Tip: the eps plots generated by WUPlot can be very large. "
	     << endl
	     << "     try converting to PDF using 'ps2pdf "
	     <<WUPlot::outputfilename<<"'"<< endl
	     << "     to drastically reduce the file size." << endl;	    
    }

}


/**
 * parse plot file and process commands 
 * 
 * \todo: change this to use the token list as a queue, shifting off
 * each token as it is processed!  It will be much more readable!
 */
void parsePlot( string filename ) {

    ifstream plotfile( filename.c_str() );
    int linenum=0;
    PlotMaker *pm;

    if (WUPlot::xdisplay == true) 
	pm = new PlotMaker( "X", "x.out" );
    else
	pm = new PlotMaker( WUPlot::outputtype, WUPlot::outputfilename );

    string line;
    int i;

    //================================================================
    // open the plot file
    //================================================================

    if (plotfile.fail()) {
	cout << "** Couldn't open '"<<filename<<"'"<<endl;
	return;
    } 
    
    
    //================================================================
    // read each line, tokenize, and determine what to do...
    //================================================================
    
    while (plotfile && !plotfile.eof()) {
	
	vector<string> tokens;

	linenum++;

	getline( plotfile, line, '\n' );

	// remove white space and comments

	i = line.find("#");
	line = line.substr(0,i);
	removeWhiteSpace( line );

	// ensure that quoted strings take up exactly one token
	// (convert whitespace into special characters). Later, quoted
	// strings must be read back with getStringLiteral to remove
	// the excess formatting and quote characters.
	
	prepareQuotedStrings( line );

	// tokenize

	tokenize( line, tokens, " \t," );
	if (tokens.size() < 1 ) continue;

	// process command: 

	try {
	    // ----------------------------------------------------------
	    // SET
	    // ----------------------------------------------------------
	    if (tokens[0] == "set") {

		expect(2, tokens);

		// -----------------------
		// axes
		// -----------------------
		if (tokens[1] == "axes") {
		    
		    expect( 4, tokens );

		    if (tokens[2] == "geometry") {
			string file = getStringLiteral( tokens, 3 );
			Image2D im2d;
			im2d.load( file );

			pm->setAxisBox( im2d );
			
		    }
		    else if (tokens[2] == "ticcolor") {
			string color = getStringLiteral(tokens, 3);
			pm->setTicMarkColorName( color );
		    }
		    else {
			throw ParseException("Unknown axis type "
						    +tokens[2]);
		    }
		    
		}
		
		// -----------------------
		// contour
		// -----------------------
		
		else if (tokens[1] == "contour") {

		    expect(4, tokens);		    

		    if (tokens[2] == "color") {
			string color = getStringLiteral(tokens, 3);
			pm->setContourColorName( color );
		    }
		    else if (tokens[2] == "style") {
			string style = getStringLiteral(tokens,3);
			pm->setContourLineStyle( style );
			
		    }
		    else if (tokens[2] == "width") {
			double width = atof( tokens[3].c_str() );
			pm->setContourLineWidth( width );
			
		    }
		    else throw ParseException("Unknown attribute "+tokens[2]);

		    
		}
		
		// -----------------------
		// marker
		// -----------------------

		else if (tokens[1] == "marker") {
		    
		    expect( 4, tokens );

		    if (tokens[2] == "type") {
			WUPlot::marker = atoi(tokens[3].c_str());
			pm->setMarkerType(atoi(tokens[3].c_str()));
		    }
		    else if (tokens[2] == "color") {
			string color = getStringLiteral(tokens, 3);
			pm->setMarkerColorName( color );
		    }
		    else if (tokens[2] == "size") {
			pm->setMarkerSize( atof(tokens[3].c_str()));
		    }
		    else throw ParseException("Unknown attribute "+tokens[2]);

		}

		// -----------------------
		// smooth
		// -----------------------

		else if (tokens[1] == "smooth") {

		    expect( 3, tokens );

		    WUPlot::smooth = atoi(tokens[2].c_str());

		}
		
		// -----------------------
		// colormap
		// -----------------------


		else if (tokens[1] == "colormap") {

		    expect( 3, tokens );

		    string colormap = getStringLiteral( tokens, 2);
		    pm->setColorMap( colormap );
		}

		// -----------------------
		// colorscale
		// -----------------------
		
		else if (tokens[1] == "colorscale") {
		    expect( 3,tokens );
		    if (tokens.size() == 3) {
			if (tokens[2]=="auto") {
			    pm->enableAutoScale(true);
			}
			else {
			    throw ParseException("Unknown scale "+tokens[2]);
			}
		    }
		    else {
			expect(4,tokens);
			pm->enableAutoScale(false);
			pm->setColorScale( atof(tokens[2].c_str()),
					   atof(tokens[3].c_str()) );
		    }
		    
		}

		// -----------------------
		// label
		// -----------------------

		else if (tokens[1] == "label") {

		    expect( 4, tokens );

		    if ( tokens[2] == "size" ) {
			pm->setLabelSize(atof( tokens[3].c_str() ));
		    }
		    else if ( tokens[2] == "angle" ) {
			pm->setLabelAngle( atof( tokens[3].c_str() ) );
		    }
		    else if (tokens[2] == "color" ){
			string name = getStringLiteral( tokens, 3 );
			pm->setLabelColorName( name );
		    }
		    else throw ParseException("Unknown attribute "+tokens[2]);

		}

		// -----------------------
		// rotation
		// -----------------------

		else if (tokens[1] == "rotation") {
		    expect( 3, tokens );
		    
		    double theta = atof(tokens[2].c_str());

		    pm->rotate( theta );

		}

		// -----------------------
		// unknown
		// -----------------------		

		else{
		    throw ParseException("Unknown attribute '"+
						tokens[1]+"'");		
		}


		
	    }
	    // ----------------------------------------------------------
	    // PLOT
	    // ----------------------------------------------------------
	    else if (tokens[0] == "plot" ) {

		expect( 2, tokens );
		
		// -----------------------
		// image
		// -----------------------
		if (tokens[1] == "image") {
		    
		    expect( 3, tokens );

		    string file = getStringLiteral( tokens, 2 );
		    Image2D im;
		    im.load(file);
		    if (WUPlot::smooth) {
			for (int i=0; i<WUPlot::smooth; i++) 
			    im.expand();
		    }
		    pm->plot( im );
		

		}

		// -----------------------
		// axes
		// -----------------------
		else if (tokens[1] == "axes" ) {
		    pm->plotAxes();
		}

		// -----------------------
		// contour
		// -----------------------
		else if (tokens[1] == "contour") {

		    double min,max,delta;

		    expect( 6, tokens );

		    min = atof( tokens[2].c_str() );
		    max = atof( tokens[3].c_str() );
		    delta = atof( tokens[4].c_str() );

		    string file = getStringLiteral( tokens, 5 );
		    Image2D im;

		    im.load(file);
		    if (WUPlot::smooth) im.expand();
		    pm->contourPlot( im, min, max, delta );
		
		}
		
		// -----------------------
		// curve
		// -----------------------

		else if (tokens[1] == "curve") {

		    expect( 3, tokens );

		    vector<double> xpoints;
		    vector<double> ypoints;
		    
		    string file = getStringLiteral( tokens, 2 );
		    
		    readCurve( file, xpoints, ypoints );

		    pm->plot( xpoints, ypoints );

		    
		}

		// -----------------------
		// list of marker points
		// -----------------------
		
		else if (tokens[1] == "points") {
		    expect( 3, tokens );

		    vector<double> xpoints;
		    vector<double> ypoints;
		    
		    string file = getStringLiteral( tokens, 2 );
		    
		    readCurve( file, xpoints, ypoints );

		    for (int i=0; i<xpoints.size(); i++) {
			pm->plotMarker( xpoints[i],ypoints[i] );
		    }

		}
		

		// -----------------------
		// title
		// -----------------------
		else if (tokens[1] == "title" ) {

		    expect( 3, tokens );

		    string title = getStringLiteral( tokens, 2 );
		
		    pm->plotTitle( title );

		}

		// -----------------------
		// subtitle
		// -------------9----------
		else if (tokens[1] == "subtitle" ) {

		    expect( 3, tokens );

		    string title = getStringLiteral( tokens, 2 );
		
		    pm->plotSubtitle( title );

		}

		// -----------------------
		// starlist
		// -----------------------

		else if (tokens[1] == "starlist") {

		    expect( 3, tokens );

		    string file = getStringLiteral( tokens, 2 );
		    int n;

		    ifstream starfile( file.c_str() );
		    if (starfile.fail()){
			throw ParseException("Couldn't open "+file);
		    }

		    Star star;
		    starfile >> n;
		    while (starfile) {
			starfile >> star;
			pm->plot( star , WUPlot::marker );
		    }
		    
		    starfile.close();

		}

		// -----------------------
		// label
		// -----------------------

		else if (tokens[1] == "label" ) {
		    
		    expect( 5, tokens );

		    double x = atof(tokens[2].c_str());
		    double y = atof(tokens[3].c_str());
		    string text = getStringLiteral( tokens, 4 );
		   
		    pm->plot( text, x,y  );

		}

		
		// -----------------------
		// ellipse
		// -----------------------

		else if (tokens[1] == "ellipse") {

		    expect( 5, tokens );

		    double x = atof( tokens[2].c_str() );
		    double y = atof( tokens[3].c_str() );
		    double r1 = atof( tokens[4].c_str() );
		    double r2 = atof( tokens[5].c_str() );
		    double theta = atof( tokens[6].c_str() );

		    pm->plotEllipse( x,y,r1,r2,theta );

		}

		// -----------------------
		// Marker
		// -----------------------

		else if (tokens[1] == "marker" ) {
		    expect( 4, tokens );
		    double x = atof( tokens[2].c_str() );
		    double y = atof( tokens[3].c_str() );
		    pm->plotMarker( x,y );
		}

		// -----------------------
		// radecgrid
		// -----------------------

		else if (tokens[1] == "radecgrid") {
		    expect( 4, tokens );
		    Image2D ragrid, decgrid;
		    string ramesh = getStringLiteral( tokens, 2 );
		    string decmesh = getStringLiteral( tokens, 3 );
		    ragrid.load(ramesh);
		    decgrid.load(decmesh);
		    
		    pm->plotRADecGrid( ragrid, decgrid );

		}

		// -----------------------
		// unknown
		// -----------------------
		   
		else {
		    throw ParseException("Unknown plot type: "
						+tokens[1]);
		}
	    

	    }
	    // ----------------------------------------------------------
	    // OTHERS
	    // ----------------------------------------------------------
	    else if (tokens[0] == "push") {
		pm->pushState();
	    }
	    else if (tokens[0] == "pop") {
		pm->popState();
	    }
	    else {
		throw ParseException("Unidentified command '"+tokens[0]+"'");
	    }
       
	}
	catch (ParseException &e) {
	    cout << "WARNING: "<<e.what()<< " at line "
		 << linenum << " (skipped)"<<endl;
	    continue;	    
	}
	catch (MildAnalysisException &e) {
	    cout << "MildAnalysisException: "<<e.what()<< " at line "
		 << linenum << " (skipped)"<<endl;
	    continue;	    
	}

    }
    
    delete pm;
    

}



/**
 * Returns string that was enclosed in quotes and which may span
 * multiple tokens. This is really just a hack so i don't have to
 * write a proper parser.
 */
string 
getStringLiteral( vector<string> &v, int start ) {
    
    string s;
    int  q2;

    for (int i=start; i<v.size(); i++) {
	s.append(v[i]);
	if (i!= v.size()-1) s.append(" ");
    }
    
    if (s[0] != '"') {
	throw ParseException("Expected a quoted string");
	return "";
    }
    
    q2 = s.find_first_of('"',1);
    
    string literal =  s.substr(1,q2-1);

    for (int i=0; i<literal.length(); i++) 
	if (literal[i]==SPACE_PLACEHOLDER) literal[i]=' ';

    return literal;

}



void 
WUPlot::getCommandLineOptions(int argc , char **argv) {

    int ret=0;
    int c;

    while ( ret != -1 ) {

	ret=getopt( argc, argv, "o:T:");
	c = ret;

	switch (c) {
	    
	case 'o':
	    outputfilename = optarg;
	    WUPlot::xdisplay=false;
	    break;

	case 'T':
	    WUPlot::outputtype = optarg;

	default:
	    break;

	}
	
    }
    
}


void
WUPlot::showUsage() {

    cout << "USAGE: wuplot [-o filename] [-T <type>] <inputfile>" << endl;
    cout << "\t-o filename       send output to specified  file"<< endl
	 << "\t                  instead of to the X display" << endl;
    cout << "\t-T type           output file format [ps is default]" << endl;

    cout << "output types: "<<endl
	 <<"X, png, pnm, gif, svg, ai, ps, cgm, fig, pcl, hpgl, regis, or tek"
	 << endl;
	

}


void expect( int n, vector<string> &v) {
    
    if (v.size()<n) {
	 throw ParseException("wrong number of arguments");
    }	        
    

}


void readCurve( string filename, vector<double> &xpoints, 
		vector<double> &ypoints ) {

    string line;
    int linenum=0;
    vector<string> tokens;
    ifstream infile( filename.c_str() );
    if (infile.fail())
	throw ParseException("Couldn't open curve "
				    +filename );

    xpoints.clear();    
    ypoints.clear();
    
    
    while (infile && !infile.eof()) {

	getline( infile, line, '\n' );
	linenum++;
	tokens.clear();
	tokenize( line, tokens, " \t," );

	if (tokens.size() == 0 ) continue;

	if (tokens.size() != 2) {
	    cout << "** skipping Bad data on line "<< linenum 
		 << " of "<<filename << endl;
	    continue;
	}

	xpoints.push_back( atof(tokens[0].c_str()) );
	ypoints.push_back( atof(tokens[1].c_str()) );
	
    }

}


void prepareQuotedStrings( string &line ) {

    vector<int> quotes;

    // find quotes

    for (int i=0; i<line.length(); i++) {
	if (line[i] == '"') quotes.push_back(i);
    }
    
    if (quotes.size() % 2 != 0) {
	throw ParseException("Missing quotation mark");
    }

    for (int i=0; i<quotes.size(); i+=2) {
	for (int j=quotes[i]; j<=quotes[i+1]; j++) {
	    if (line[j] == ' ') {
		line[j]=SPACE_PLACEHOLDER;
	    }
	}
    }
    

}