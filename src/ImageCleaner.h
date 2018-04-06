#ifndef IMAGECLEANER_H
#define IMAGECLEANER_H

#include <vector>
#include "Types.h"
#include "PedestalFinder.h"

class Camera;

/**
 *  Abstract interface to all cleaning routines. Takes an image as
 *  input and spits out a list of pixels which pass the cleaning
 *  process.  This class can be extended to do more interesting things
 *  (such as island cleaning, etc).
 */
class ImageCleaner {

public: 
    ImageCleaner( Camera &camera ){ _cam = &camera; }
    virtual std::vector<int>& getCleanPixels( const Array_t &image )=0;
					   
protected:
    
    Camera *_cam;
    

};


/**
 * The default cleaner uses a simple 2 threshold process.
 */ 
class DefaultImageCleaner : public ImageCleaner {

public: 
    DefaultImageCleaner(Camera &camera, vector<Pedestal> &peds,
			double picture, double boundary);
    std::vector<int>& getCleanPixels( const Array_t &image );
    
private:
    
    double _picture_thresh;
    double _boundary_thresh;
    std::vector<Pedestal> *_peds;
    std::vector<int> _cleanpixels;

};

class ThresholdImageCleaner : public ImageCleaner {

 public:
    ThresholdImageCleaner( Camera &camera, double threshold ) 
	: ImageCleaner(camera),_threshold(threshold){;}

    std::vector<int>& getCleanPixels( const Array_t &image );
    
    void setThreshold( double t ) {_threshold=t;}


 private:
    double _threshold;
    std::vector<int> _cleanpixels;

};


#endif
