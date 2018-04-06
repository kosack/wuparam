//
// wu2root  - convert WUParam ntuples to ROOT trees
//
// Karl Kosack, from modified version of Paul Rebillot's script

#include <iostream>
#include <getopt.h>
#include <cstdlib>
#include <fstream>
#include <string>
#include <exception>
#include <stdexcept>
#include <unistd.h>
#include <signal.h>

// ROOT stuff:

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>

// WUPARAM stuff:

#include "Types.h"
#include "Exceptions.h"
#include "Config.h"
#include "ImageAnalyzer.h"
#include "ParamDataReader.h"
#include "Cutter.h"
#include "Log.h"

using namespace std;

namespace WU2Root {
    void getCommandLineOptions(int argc , char **argv);
    void showUsage();
}

int
main( int argc , char **argv ) {

    WU2Root::getCommandLineOptions(argc,argv);
    
    ParamDataReader *reader;
    HillasParameterization param;

    if (argc-optind != 1 ) {
	WU2Root::showUsage();
	exit(1);
    }

    string ntuplefilename, newfilename;
    ntuplefilename = argv[optind];
    newfilename = ntuplefilename+".root";

    TFile *tfile = new TFile(newfilename.c_str(), "RECREATE"); 
    
    
    // ===============================================================
    //define branches:
    // ===============================================================

    TTree *tree = new TTree("params",
			    "WUParam parameter output");

    // simulation parameters:

    tree->Branch("sim_event_num",&param.sim.event_number,"sim_event_num/I");
    tree->Branch("sim_primary",&param.sim.primary_type,"sim_primary/I");
    tree->Branch("sim_energy",&param.sim.primary_energy,"sim_energy/D");
    tree->Branch("sim_impact_x",&param.sim.impact_parameter.x,
		 "sim_impact_x/D");
    tree->Branch("sim_impact_y",&param.sim.impact_parameter.y,
		 "sim_impact_y/D");
    tree->Branch("sim_dircos_x",&param.sim.direction_cos.x,
		 "sim_dircos_x/D");
    tree->Branch("sim_dircos_y",&param.sim.direction_cos.y,
		 "sim_dircos_y/D");
    
    // global information
    
    tree->Branch("tele_id",&param.telescope_id,"tele_id/s");
    tree->Branch("osctime",&param.osctime,"osctime/D");
    tree->Branch("gpstime",&param.gpstime,"gpstime/D");
    tree->Branch("event_num",&param.event_number,"gpstime/I");
    
    // HillasParameterization stuff:

    tree->Branch("cen_x",&param.centroid.x,"cen_x/D");
    tree->Branch("cen_y",&param.centroid.y,"cen_y/D");
    tree->Branch("poo_Ax",&param.point_of_origin_a.x,"poo_Ax/D");
    tree->Branch("poo_Ay",&param.point_of_origin_a.y,"poo_Ay/D");
    tree->Branch("poo_Bx",&param.point_of_origin_b.x,"poo_Bx/D");
    tree->Branch("poo_By",&param.point_of_origin_b.y,"poo_By/D");
    tree->Branch("length",&param.length,"length/D");
    tree->Branch("width",&param.width,"width/D");
    tree->Branch("miss",&param.miss,"miss/D");
    tree->Branch("distance",&param.distance,"distance/D");
    tree->Branch("azwidth",&param.azwidth,"azwidth/D");
    tree->Branch("alpha",&param.alpha,"alpha/D");
    tree->Branch("psi",&param.psi,"psi/D");
    tree->Branch("phi",&param.phi,"phi/D");
    tree->Branch("size",&param.size,"size/D");
    tree->Branch("max1",&param.max[0],"max1/D");
    tree->Branch("max2",&param.max[1],"max2/D");
    tree->Branch("max3",&param.max[2],"max3/D");
    tree->Branch("i_max1",&param.index_of_max[0],"i_max1/I");
    tree->Branch("i_max2",&param.index_of_max[1],"i_max2/I");
    tree->Branch("i_max3",&param.index_of_max[2],"i_max3/I");
    tree->Branch("frac1",&param.frac[0],"frac1/D");
    tree->Branch("frac2",&param.frac[1],"frac2/D");
    tree->Branch("frac3",&param.frac[2],"frac3/D");
    tree->Branch("asym",&param.asymmetry,"asym/D");
    tree->Branch("len_size",&param.length_over_size,"len_size/D");
    tree->Branch("zenith_rad",&param.zenith,"zenith_rad/D");
    tree->Branch("energy_est",&param.energy_estimate,"energy_est/D");
    
    try {
	
    
	try {
	    
	    reader = new ParamDataReader( ntuplefilename );
	    
	    
	    while (1) {
		reader->getNextEventRecord( param );
		tree->Fill();
	    }
	    
	}
	catch (EOFException &e){/* end of file */ }
	tree->Write();

	
    }
    catch (AnalysisException &e){
	cout << "conversion failed: "<< e.what() << endl;
    }
    
    delete reader;
    delete tree;
    delete tfile; 
   
}



void 
WU2Root::getCommandLineOptions(int argc , char **argv) {

    char c=0;

    while ( c != -1 && c!= 255) {

	c=getopt( argc, argv, "");

	switch (c) {

	default:
	    break;
	}
	
    }

}

void WU2Root::showUsage() {

    cerr << "USAGE: wu2root <ntuple>" << endl;
    cerr << "\tOutput will be to *.ntuple.root"<<endl;

}


