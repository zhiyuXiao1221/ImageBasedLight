#ifndef __FILTERING_H
#define __FILTERING_H

// Assignment 4
// Filtering and convolution

#include "floatimage.h"
#include <iostream>
#include <math.h>
#include "array.h"

using namespace std;

class Filter
{
public:
    std::vector<float> kernel;
    int width;
    int height;

    // function to convolve your filter with an image
    FloatImage Convolve(const FloatImage &im, bool clamp = true) const;

    // Accessors of the filter values
    const float &operator()(int x, int y) const;
    float &operator()(int x, int y);

    //Constructor
    Filter(const vector<float> &fData, const int &fWidth, const int &fHeight);
    Filter(const int &fWidth, const int &fHeight); // kernel with all zero

    // Destructor. Because there is no explicit memory management here, this doesn't do anything
    ~Filter();

    // The following are functions and variables that are not accessible from outside the class
private:
};

// Box Blurring
FloatImage boxBlur(const FloatImage &im, const int &k, bool clamp = true);
FloatImage boxBlur_filterClass(const FloatImage &im, const int &k, bool clamp = true);

// Gradient Filter
FloatImage gradientMagnitude(const FloatImage &im, bool clamp = true);

// Gaussian Blurring
vector<float> gauss1DFilterValues(float sigma, float truncate);
vector<float> gauss2DFilterValues(float sigma, float truncate);
FloatImage gaussianBlur_horizontal(const FloatImage &im, float sigma, float truncate = 3.0, bool clamp = true);
FloatImage gaussianBlur_seperable(const FloatImage &im, float sigma, float truncate = 3.0, bool clamp = true);
FloatImage gaussianBlur_2D(const FloatImage &im, float sigma, float truncate = 3.0, bool clamp = true);
vector<float> gaussWeights(float mu, float sigma);

// Sharpen an FloatImage
FloatImage unsharpMask(const FloatImage &im, float sigma, float truncate = 3.0, float strength = 1.0, bool clamp = true);

// Bilaterial Filtering
FloatImage bilateral(const FloatImage &im, float sigmaRange = 0.1, float sigmaDomain = 1.0, float truncateDomain = 3.0, bool clamp = true);
FloatImage bilaYUV(const FloatImage &im, float sigmaRange = 0.1, float sigmaY = 1.0, float sigmaUV = 4.0, float truncateDomain = 3.0, bool clamp = true);
FloatImage fastBilateral(const FloatImage &im, float sigmaRange = 0.1, float sigmaDomain = 1.0, int truncateDomain = 5.0, float samplingD = 1.0, float samplingR = 0.05);
float trilinear_interpolation(const Array_3D<float> &array, float x, float y, float z);
float clamp(int min_value, int max_value, int x);
// Return impulse image of size kxkx1
FloatImage impulseImg(const int &k);

// Laplacian
vector<FloatImage> gaussPyramid(const FloatImage &im, int levels);
vector<FloatImage> laplacianPyramid(const FloatImage &im,  int levels);
FloatImage localLaplacianFilter(const FloatImage im, int levels, float sigma, float alpha, float beta, int channels);

#endif
