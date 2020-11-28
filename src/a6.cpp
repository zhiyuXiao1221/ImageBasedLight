#include "a6.h"
#include "floatimage.h"
#include <iostream>
#include <math.h>

using namespace std;

// create a new image that is k times bigger than the input by using nearest neighbor interpolation
FloatImage scaleNN(const FloatImage &im, float factor)
{

	int ys, xs;

	// get new FloatImage
	int nWidth = floor(factor * im.width());
	int nHeight = floor(factor * im.height());
	FloatImage im2(nWidth, nHeight, im.channels());

	// get appropriate values (smartAccessor is probably overkill here)
	for (int x = 0; x < nWidth; x++)
		for (int y = 0; y < nHeight; y++)
			for (int z = 0; z < im.channels(); z++)
			{
				ys = round(1 / factor * y);
				xs = round(1 / factor * x);
				im2(x, y, z) = im.smartAccessor(xs, ys, z, true);
			}

	// return new image
	return im2;

}

// helper for upsampleing, scale without interpolation
vector<FloatImage> scaleU(const FloatImage &im, float factor)
{
	// get new FloatImage
	int nWidth = floor(factor * im.width());
	int nHeight = floor(factor * im.height());
	FloatImage im2(nWidth, nHeight, im.channels());
	im2.clear();
	FloatImage im3(im2);
	
	float ys, xs;

	// get appropriate values (smartAccessor is probably overkill here)
	for (int x = 0; x < nWidth; x++)
		for (int y = 0; y < nHeight; y++)
			for (int z = 0; z < im.channels(); z++)
			{
				ys = 1. / factor * y;
				xs = 1. / factor * x;
				if (round(ys) == ys || round(xs) == xs){
					im2(x, y, z) = im.smartAccessor(xs, ys, z, true);
					im3(x, y, z) = 1;
				}


			}

	// return new image
	return vector<FloatImage>{im2, im3};

}

// using bilinear interpolation to assign the value of a location from its neighboring pixel values
float interpolateLin6(const FloatImage &im, float x, float y, int z, bool clamp)
{

	// get the extreme points
	int xf = floor(x);
	int yf = floor(y);
	int xc = xf + 1;
	int yc = yf + 1;

	// compute the distances of the point to the floor-extreme point
	float yalpha = y - yf;
	float xalpha = x - xf;

	// obtain the values at those points
	float tl = im.smartAccessor(xf, yf, z, clamp);
	float tr = im.smartAccessor(xc, yf, z, clamp);
	float bl = im.smartAccessor(xf, yc, z, clamp);
	float br = im.smartAccessor(xc, yc, z, clamp);

	// compute the interpolations on the top and bottom
	float topL = tr * xalpha + tl * (1.0 - xalpha);
	float botL = br * xalpha + bl * (1.0 - xalpha);

	// compute the overall interpolation
	float retv = botL * yalpha + topL * (1.0 - yalpha);

	// return final floar value
	return retv;

}

// create a new image that is k times bigger than the input by using bilinear interpolation
FloatImage scaleLin(const FloatImage &im, float factor)
{

	float ys, xs;

	// get new FloatImage
	
	int nWidth = floor(factor * im.width());
	int nHeight = floor(factor * im.height());
	FloatImage im2(nWidth, nHeight, im.channels());

	// get appropriate values
	for (int x = 0; x < nWidth; x++)
		for (int y = 0; y < nHeight; y++)
			for (int z = 0; z < im.channels(); z++)
			{

				ys = 1 / factor * y;
				xs = 1 / factor * x;

				// im2(x, y, z) = interpolateLin(im, xs, ys, z);
				im2(x, y, z) = interpolateCubic(im, xs, ys, z, true);
			}

	// return new image
	return im2;

}

// using bicubic interpolation to assign the value of a location from its neighboring pixel values
float interpolateCubic(const FloatImage &im, float x, float y, int z, bool clamp)
{
	// Hint: use smartAccessor() to handle coordinates outside the image
	int xi = floor(x);
	int yi = floor(y);

	int xj = xi + 1;
	int yj = yi + 1;

	float dx = x - xi;
	float dy = y - yi;

	float p00, p01, p02, p03
	, p10, p11, p12, p13
	, p20, p21, p22, p23
	, p30, p31, p32, p33;

	// get 16 pixels around point
	p00 = im.smartAccessor(xi-1, yi-1, z, clamp);
    p01 = im.smartAccessor(xi, yi-1, z, clamp);
    p02 = im.smartAccessor(xj, yi-1, z, clamp);
    p03 = im.smartAccessor(xj+1, yi-1, z, clamp);

    p10 = im.smartAccessor(xi-1, yi, z, clamp);
    p11 = im.smartAccessor(xi  , yi, z, clamp);
    p12 = im.smartAccessor(xj, yi, z, clamp);
    p13 = im.smartAccessor(xj+1, yi, z, clamp);

    p20 = im.smartAccessor(xi-1, yj, z, clamp);
    p21 = im.smartAccessor(xi, yj, z, clamp);
    p22 = im.smartAccessor(xj, yj, z, clamp);
    p23 = im.smartAccessor(xj+1, yj, z, clamp);

    p30 = im.smartAccessor(xi-1, yj+1, z, clamp);
    p31 = im.smartAccessor(xi, yj+1, z, clamp);
    p32 = im.smartAccessor(xj, yj+1, z, clamp);
    p33 = im.smartAccessor(xj+1, yj+1, z, clamp);

	const float A = -0.75f;

	float coeffs0x = ((A*(dx + 1) - 5*A)*(dx + 1) + 8*A)*(dx + 1) - 4*A;
	float coeffs1x = ((A + 2)*dx - (A + 3))*dx*dx + 1;
	float coeffs2x = ((A + 2)*(1 - dx) - (A + 3))*(1 - dx)*(1 - dx) + 1;
	float coeffs3x = 1.f - coeffs0x - coeffs1x - coeffs2x;

	float coeffs0y = ((A*(dy + 1) - 5*A)*(dy + 1) + 8*A)*(dy + 1) - 4*A;
	float coeffs1y = ((A + 2)*dy - (A + 3))*dy*dy + 1;
	float coeffs2y = ((A + 2)*(1 - dy) - (A + 3))*(1 - dy)*(1 - dy) + 1;
	float coeffs3y = 1.f - coeffs0y - coeffs1y - coeffs2y;

  	float x0 = coeffs0x*p00 + coeffs1x*p01 + coeffs2x*p02 + coeffs3x*p03;
  	float x1 = coeffs0x*p10 + coeffs1x*p11 + coeffs2x*p12 + coeffs3x*p13;
	float x2 = coeffs0x*p20 + coeffs1x*p21 + coeffs2x*p22 + coeffs3x*p23;
	float x3 = coeffs0x*p30 + coeffs1x*p31 + coeffs2x*p32 + coeffs3x*p33;

	//return final float value
	return coeffs0y*x0 + coeffs1y*x1 + coeffs2y*x2 + coeffs3y*x3;
}

// rotate an image around its center by theta
FloatImage rotate(const FloatImage &im, float theta)
{

	// create rotation matrix
	float yR, xR;

	// center around which to rotate
	float centerX = (im.width() - 1.0) / 2.0;
	float centerY = (im.height() - 1.0) / 2.0;

	// get new image
	FloatImage imR(im.width(), im.height(), im.channels());

	// get appropriate values
	for (int x = 0; x < im.width(); x++)
		for (int y = 0; y < im.height(); y++)
			for (int z = 0; z < im.channels(); z++)
			{
				// compute the x and y values from the original image
				xR = ((float) x - centerX) * cos(theta) + (centerY - (float) y) * sin(theta) + centerX;
				yR = centerY - (-((float) x - centerX) * sin(theta) + (centerY - (float) y) * cos(theta));

				// interpolate the point
				imR(x, y, z) = interpolateLin6(im, xR, yR, z);
			}

	// return new image
	return imR;

}