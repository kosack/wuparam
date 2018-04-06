
#include "ImageCleaner.h"
#include "Camera.h"
#include "PedestalFinder.h"
#include <vector>
#include <iostream>

using namespace std;

/**
 * The default image cleaner uses picture and boundary theshold values
 * (expressed in standard deviations from the pedestal) to determine
 * the accepted pixels. All tubes above the picture threshold are
 * accepted. Tubes above the boundary threshold with at least one
 * neighboring "picture" tube are also accepted.
 *
 * \param camera reference to the camera object to use for pixel
 * positions and neighborlist.
 * \param peds reference to the Pedestal array for the current run.
 * \param picture the picture threshold in std deviations.
 * \param boundary the boundary threhold in std deviations.
 *
 */
DefaultImageCleaner::
DefaultImageCleaner(Camera &camera, vector<Pedestal> &peds, 
		    double picture,double boundary) 
    : ImageCleaner(camera),_picture_thresh(picture),_boundary_thresh(boundary)
{
    
    _peds = &peds;
    _cleanpixels.reserve(500); // minimize resizing of cleanpixels

}


vector<int>& 
DefaultImageCleaner::
getCleanPixels( const Array_t &image ) {
    
    list<int> *neighs;
    list<int>::iterator iter;
    valarray<bool> picture(image.size());
    vector<Pedestal> &peds = *_peds;

    int i;

    _cleanpixels.clear();

    for (i=0; i<(int)picture.size(); i++) {
	picture[i] = false;
    }

    // Put the "picture" tubes into the list...

    for (i=_cam->getFirstPixel(); i<_cam->getLastPixel(); i++) {
	
	if ( peds[i].type == Pedestal::GOOD &&
	    image[i] > _picture_thresh*peds[i].dispersion  ){ 

	    picture[i] = true;
	    _cleanpixels.push_back(i);
	}

    }

    // The "boundary" tubes must be over the boundary threshold and
    // have a picture tube as a neighbor
    
    for (i=_cam->getFirstPixel(); i<_cam->getLastPixel(); i++) {

	if (picture[i] == false) {
	    if ( (image[i] > _boundary_thresh*peds[i].dispersion)
		 && peds[i].type == Pedestal::GOOD  ) {
		
		neighs = _cam->getNeighborListOfPixel( i );

		for (iter=neighs->begin(); iter != neighs->end(); iter++) {
		    if (picture[*iter]) {
			
			// this tube passes, so add it to the list
			_cleanpixels.push_back(i);
			break;

		    }
		}
	    }
	}
    }

    return _cleanpixels;

}




vector<int>& 
ThresholdImageCleaner::
getCleanPixels( const Array_t &image ) {


    _cleanpixels.clear();
    
    for (int i=_cam->getFirstPixel(); i<_cam->getLastPixel(); i++) {
	if ( image[i] > _threshold ) 
	    _cleanpixels.push_back(i);
    }
    
    return _cleanpixels;
    
}
