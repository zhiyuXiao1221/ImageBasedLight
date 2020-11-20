#include "a6.h"
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

// using bilinear interpolation to assign the value of a location from its neighboring pixel values
float interpolateLin(const FloatImage &im, float x, float y, int z, bool clamp)
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

				im2(x, y, z) = interpolateLin(im, xs, ys, z);
			}

	// return new image
	return im2;

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
				imR(x, y, z) = interpolateLin(im, xR, yR, z);
			}

	// return new image
	return imR;

}