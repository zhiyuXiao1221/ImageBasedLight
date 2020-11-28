/*
    CS 89/189 Computational Aspects of Digital Photography C++ basecode.

    Adapted from MIT's 6.815/6.865 basecode, written and designed by:
        Frédo Durand
        Katherine L. Bouman
        Gaurav Chaurasia
        Adrian Vasile Dalca
        Neal Wadhwa

    With additions & modifications by
        Wojciech Jarosz
    for Dartmouth's CS 89/189.
*/

// filtering.cpp
// Assignment 4

#include "filtering.h"
#include "a2.h"
#include <math.h>
#include "array.h"
#include "utils.h"

using namespace std;

/**************************************************************
 //            IMAGE CONVOLUTION AND FILTERING               //
 *************************************************************/

// convolve an image with a box filter of size k by k
FloatImage boxBlur(const FloatImage &im, const int &k, bool clamp)
{
	// create a new empty image
	FloatImage imFilter(im.width(), im.height(), im.channels());
	int sideSize = int((k - 1.0) / 2.0);
	float normalizer = 1.0 / float(k * k);
	float sumVal;

	// for every pixel in the image
	for (int x = 0; x < imFilter.width(); x++)
		for (int y = 0; y < imFilter.height(); y++)
			for (int z = 0; z < imFilter.channels(); z++)
			{
				sumVal = 0.0;

				for (int xBox = -sideSize; xBox <= sideSize; xBox++)
					for (int yBox = -sideSize; yBox <= sideSize; yBox++)
					{
						sumVal += im.smartAccessor(x - xBox, y - yBox, z, clamp);
					}

				// assign the pixel the value from convolution
				imFilter(x, y, z) = sumVal * normalizer;
			}

	return imFilter;
}

// reimeplement the box filter using the filter class.
// check that your results math those in the previous function "boxBlur"
FloatImage boxBlur_filterClass(const FloatImage &im, const int &k, bool clamp)
{
	vector<float> fData(k * k, 1.0 / (k * k));
	Filter boxFilter(fData, k, k);
	FloatImage imFilter = boxFilter.Convolve(im, clamp);
	return imFilter;
}

// uses a Sobel kernel to compute the horizontal and vertical
// components of the gradient of an image and returns the gradient magnitude.
FloatImage gradientMagnitude(const FloatImage &im, bool clamp)
{
	// sobel filtering in x direction
	Filter sobelX(3, 3);
	sobelX(0, 0) = -1.0;
	sobelX(1, 0) = 0.0;
	sobelX(2, 0) = 1.0;
	sobelX(0, 1) = -2.0;
	sobelX(1, 1) = 0.0;
	sobelX(2, 1) = 2.0;
	sobelX(0, 2) = -1.0;
	sobelX(1, 2) = 0.0;
	sobelX(2, 2) = 1.0;

	FloatImage imSobelX = sobelX.Convolve(im, clamp);

	// sobel filtering in y direction
	Filter sobelY(3, 3);
	sobelY(0, 0) = -1.0;
	sobelY(1, 0) = -2.0;
	sobelY(2, 0) = -1.0;
	sobelY(0, 1) = 0.0;
	sobelY(1, 1) = 0.0;
	sobelY(2, 1) = 0.0;
	sobelY(0, 2) = 1.0;
	sobelY(1, 2) = 2.0;
	sobelY(2, 2) = 1.0;

	FloatImage imSobelY = sobelY.Convolve(im, clamp);

	// squared magnitude
	FloatImage magnitude = imSobelX * imSobelX + imSobelY * imSobelY;

	// take the square root
	for (int i = 0; i < magnitude.size(); i++)
	{
		magnitude(i) = sqrt(magnitude(i));
	}

	return magnitude;
}

// create a vector containing the normalized values in a 1D Gaussian filter
vector<float> gauss1DFilterValues(float sigma, float truncate)
{
	// calculate the size of the filter
	int offset = int(ceil(truncate * sigma));
	int filterSize = 2 * offset + 1;

	vector<float> fData(filterSize, 0);

	// compute the un-normalized value of the gaussian
	float normalizer = 0.0;
	for (int i = 0; i < filterSize; i++)
	{
		fData[i] = exp(-pow(i - offset, 2) / (2.0 * pow(sigma, 2)));
		normalizer += fData[i];
	}

	// normalize
	for (int i = 0; i < filterSize; i++)
		fData[i] /= normalizer;

	return fData;
}

// 1D Gaussian for HDR weights
vector<float> gaussWeights(float mu, float sigma)
{
	// calculate the size of the filter
	int filterSize = 2 * mu + 1;

	vector<float> fData(filterSize, 0);

	for (int i = 0; i < filterSize; i++) // for 0-255 range
	{
		fData[i] = 127 * exp(-pow(i - mu, 2) / (2.0 * pow(sigma, 2)));
	}

	return fData;
}

// blur across the rows of an image
FloatImage gaussianBlur_horizontal(const FloatImage &im, float sigma, float truncate, bool clamp)
{
	// filter in the x direction
	vector<float> fData = gauss1DFilterValues(sigma, truncate);
	Filter gaussX(fData, fData.size(), 1);
	FloatImage imFilter = gaussX.Convolve(im, clamp);
	return imFilter;
}

// create a vector containing the normalized values in a 2D Gaussian filter
vector<float> gauss2DFilterValues(float sigma, float truncate)
{
	// compute the filter size
	int offset = int(ceil(truncate * sigma));
	int k = 2 * offset + 1;
	int filterSize = k * k;

	vector<float> fData(filterSize, 0);

	int count = 0;
	float normalizer = 0.0;

	// compute the unnormalized value of the gaussian and put it in a row-major vector
	for (int i = 0; i < k; i++)
	{
		for (int j = 0; j < k; j++)
		{
			fData[count] = exp(-(pow(i - offset, 2) + pow(j - offset, 2)) / (2.0 * pow(sigma, 2)));
			normalizer += fData[count];
			count++;
		}
	}

	// normalize
	for (int i = 0; i < filterSize; i++)
		fData[i] /= normalizer;

	return fData;
}

// Blur an image with a full  full 2D rotationally symmetric Gaussian kernel
FloatImage gaussianBlur_2D(const FloatImage &im, float sigma, float truncate, bool clamp)
{
	// blur using a 2D gaussian filter
	vector<float> fData = gauss2DFilterValues(sigma, truncate);
	int k = sqrt(fData.size());
	Filter gauss(fData, k, k);
	FloatImage imFilter = gauss.Convolve(im, clamp);

	return imFilter;
}

// Use principles of seperabiltity to blur an image using 2 1D Gaussian Filters
FloatImage gaussianBlur_seperable(const FloatImage &im, float sigma, float truncate, bool clamp)
{
	// blur using 2, 1D filters in the x and y directions
	vector<float> fData = gauss1DFilterValues(sigma, truncate);
	Filter gaussX(fData, fData.size(), 1);
	Filter gaussY(fData, 1, fData.size());
	FloatImage imFilter = gaussX.Convolve(im, clamp);
	imFilter = gaussY.Convolve(imFilter, clamp);

	return imFilter;
}
// sharpen an image
FloatImage unsharpMask(const FloatImage &im, float sigma, float truncate, float strength, bool clamp)
{
	// get the low pass image and subtract it from the original image to get the high pass image
	FloatImage lowPass = gaussianBlur_seperable(im, sigma, truncate, clamp);
	FloatImage highPass = im - lowPass;
	FloatImage sharp = im + strength * highPass;

	return sharp;
}

// Denoise an image using bilateral filtering
FloatImage bilateral(const FloatImage &im, float sigmaRange, float sigmaDomain, float truncateDomain, bool clamp)
{
	FloatImage imFilter(im.width(), im.height(), im.channels());

	// calculate the filter size
	int offset = int(ceil(truncateDomain * sigmaDomain));
	int sizeFilt = 2 * offset + 1;
	float sumVal, factorRangeExp, normalizer, factorDomain, factorRange;

	// for every pixel in the image
	for (int x = 0; x < imFilter.width(); x++)
	{
		for (int y = 0; y < imFilter.height(); y++)
		{
			for (int z = 0; z < imFilter.channels(); z++)
			{

				// initilize normalizer and sum value to 0 for every pixel location
				normalizer = 0.0;
				sumVal = 0.0;

				for (int xFilter = 0; xFilter < sizeFilt; xFilter++)
				{
					for (int yFilter = 0; yFilter < sizeFilt; yFilter++)
					{

						// calculate the distance between the 2 pixels (in range)
						factorRangeExp = 0.0;
						for (int z1 = 0; z1 < imFilter.channels(); z1++)
						{
							factorRangeExp +=
								pow(im.smartAccessor(x + xFilter - offset, y + yFilter - offset, z1, clamp) - im.smartAccessor(x, y, z1, clamp), 2);
						}

						// calculate the exonentiated weighting factor from the domain and range
						factorDomain =
							exp(-(pow(xFilter - offset, 2) + pow(yFilter - offset, 2)) / (2.0 * pow(sigmaDomain, 2)));
						factorRange = exp(-factorRangeExp / (2.0 * pow(sigmaRange, 2)));

						// increase the normalizer by the weighting amount
						normalizer += factorDomain * factorRange;
						sumVal += factorDomain * factorRange * im.smartAccessor(x + xFilter - offset, y + yFilter - offset, z, clamp);
					}
				}

				// set pixel in filtered image to weighted sum of values in the filter region
				imFilter(x, y, z) = sumVal / normalizer;
			}
		}
	}

	return imFilter;
}
// Fast bilateral based on a truncated kernel
// In implementing this method, we used the following information:
// http://people.csail.mit.edu/sparis/publi/2006/tr/Paris_06_Fast_Bilateral_Filter_MIT_TR_low-res.pdf
//http://people.csail.mit.edu/sparis/bf/#code
FloatImage fastBilateral(const FloatImage &im, float sigmaRange, float sigmaDomain, int truncateDomain, float samplingD, float samplingR)
{
	FloatImage output(im.width(), im.height(), im.channels());
	//Step1:Downsample
	int width = im.width();
	int height = im.height();
	int kernelSize = int(ceil(truncateDomain * sigmaDomain));
	float sigma_r = sigmaRange / samplingR;
	float sigma_s = sigmaDomain / samplingD;
	int padding_xy = (kernelSize - 1) / 2;
	int padding_z = padding_xy;
	for (int channel = 0; channel < im.channels(); channel++)
	{
		float im_Min = im(0, 0, channel);
		float im_Max = im(0, 0, channel);
		for (int i = 0; i < width; i++)
		{
			for (int j = 0; j < height; j++)
			{
				if (im(i, j, channel) < im_Min)
					im_Min = im(i, j, channel);
				if (im(i, j, channel) > im_Max)
					im_Max = im(i, j, channel);
			}
		}
		float im_delta = im_Max - im_Min;
		int downSample_width = ((width - 1) / samplingD) + 1 + 2 * padding_xy;
		int downSample_height = ((height - 1) / samplingD) + 1 + 2 * padding_xy;
		int downSample_depth = (im_delta / samplingR) + 1 + 2 * padding_z;
		Array_3D<float> w(downSample_width, downSample_height, downSample_depth, 0.0);
		Array_3D<float> iw(downSample_width, downSample_height, downSample_depth, 0.0);
		for (int i = 0; i < width; i++)
		{
			for (int j = 0; j < height; j++)
			{
				int downSample_x = (1.0 * i / samplingD + 0.5) + padding_xy;
				int downSample_y = (1.0 * j / samplingD + 0.5) + padding_xy;
				int downSample_z = ((im(i, j, channel) - im_Min) / samplingR + 0.5) + padding_z;

				w(downSample_x, downSample_y, downSample_z) += 1.0;
				iw(downSample_x, downSample_y, downSample_z) += im(i, j, channel);
			}
		}
		Array_3D<float> kernel(kernelSize, kernelSize, kernelSize, 0.0);
		//get the midpoint coordinate of kernel
		int mid = (kernelSize - 1) / 2;
		float sumVal = 0.0;
		for (int i = 0; i < kernelSize; i++)
		{
			for (int j = 0; j < kernelSize; j++)
			{
				for (int k = 0; k < kernelSize; k++)
				{
					float rr = pow(i - mid, 2) / pow(sigma_s, 2) + pow(j - mid, 2) / pow(sigma_s, 2) + pow(k - mid, 2) / pow(sigma_r, 2);
					kernel(i, j, k) = exp(-rr * 0.5);
					sumVal += exp(-rr * 0.5);
				}
			}
		}
		//normalize the weights
		for (int i = 0; i < kernelSize; i++)
		{
			for (int j = 0; j < kernelSize; j++)
			{
				for (int k = 0; k < kernelSize; k++)
				{
					kernel(i, j, k) /= sumVal;
				}
			}
		}
		//Step2: Three-dimensional convolution
		//Convolve both iw and w
		for (int i = padding_xy; i < downSample_width - padding_xy; i++)
		{
			for (int j = padding_xy; j < downSample_height - padding_xy; j++)
			{
				for (int k = padding_z; k < downSample_depth - padding_z; k++)
				{
					int minX = i - padding_xy;
					int minY = j - padding_xy;
					int minZ = k - padding_z;
					float val1 = 0.0, val2 = 0.0;
					for (int m = 0; m < kernelSize; m++)
					{
						for (int n = 0; n < kernelSize; n++)
						{
							for (int p = 0; p < kernelSize; p++)
							{
								val1 += kernel(m, n, p) * iw(minX + m, minY + n, minZ + p);
								val2 += kernel(m, n, p) * w(minX + m, minY + n, minZ + p);
							}
						}
					}
					iw(i, j, k) = val1;
					w(i, j, k) = val2;
				}
			}
		}

		//Step 3:triliner interpolation
		for (int i = 0; i < width; i++)
		{
			for (int j = 0; j < height; j++)
			{
				float Z = im(i, j, channel) - im_Min;
				float x = i / samplingD + padding_xy;
				float y = j / samplingD + padding_xy;
				float z = Z / samplingR + padding_z;
				float IW = trilinear_interpolation(iw, x, y, z);
				float W = trilinear_interpolation(w, x, y, z);
				output(i, j, channel) = IW / (W + exp(-10));
			}
		}
	}

	return output;
}

float trilinear_interpolation(const Array_3D<float> &array, float x, float y, float z)
{
	int x_size = array.x_size();
	int y_size = array.y_size();
	int z_size = array.z_size();

	int x_index = clamp(0, x_size - 1, x);
	int xx_index = clamp(0, x_size - 1, x_index + 1);

	int y_index = clamp(0, y_size - 1, y);
	int yy_index = clamp(0, y_size - 1, y_index + 1);

	int z_index = clamp(0, z_size - 1, z);
	int zz_index = clamp(0, z_size - 1, z_index + 1);

	float x_alpha = x - x_index;
	float y_alpha = y - y_index;
	float z_alpha = z - z_index;

	return (1.0f - x_alpha) * (1.0f - y_alpha) * (1.0f - z_alpha) * array(x_index, y_index, z_index) +
		   x_alpha * (1.0f - y_alpha) * (1.0f - z_alpha) * array(xx_index, y_index, z_index) +
		   (1.0f - x_alpha) * y_alpha * (1.0f - z_alpha) * array(x_index, yy_index, z_index) +
		   x_alpha * y_alpha * (1.0f - z_alpha) * array(xx_index, yy_index, z_index) +
		   (1.0f - x_alpha) * (1.0f - y_alpha) * z_alpha * array(x_index, y_index, zz_index) +
		   x_alpha * (1.0f - y_alpha) * z_alpha * array(xx_index, y_index, zz_index) +
		   (1.0f - x_alpha) * y_alpha * z_alpha * array(x_index, yy_index, zz_index) +
		   x_alpha * y_alpha * z_alpha * array(xx_index, yy_index, zz_index);
}
float clamp(int min_value, int max_value, int x)
{
	return std::max(std::min(x, (max_value)), (min_value));
}

// Bilaterial Filter an image seperately for the Y and UV components of an image
FloatImage bilaYUV(const FloatImage &im, float sigmaRange, float sigmaY, float sigmaUV, float truncateDomain,
				   bool clamp)
{
	//convert from RGB to YUV
	FloatImage imYUV = rgb2yuv(im);

	FloatImage bilY = bilateral(imYUV, sigmaRange, sigmaY, truncateDomain, clamp);
	FloatImage bilUV = bilateral(imYUV, sigmaRange, sigmaUV, truncateDomain, clamp);

	// put the Y and UV parts of the image back into one image
	for (int i = 0; i < im.width(); i++)
	{
		for (int j = 0; j < im.height(); j++)
		{
			imYUV(i, j, 0) = bilY(i, j, 0);
			imYUV(i, j, 1) = bilUV(i, j, 1);
			imYUV(i, j, 2) = bilUV(i, j, 2);
		}
	}

	// convert from YUV back to RGB
	FloatImage bilRGB = yuv2rgb(imYUV);

	return bilRGB;
}

/**************************************************************
 //                 FILTER CLASS FUNCTIONS                  //
 *************************************************************/

// write a convolution function for the filter class
FloatImage Filter::Convolve(const FloatImage &im, bool clamp) const
{
	FloatImage imFilter(im.width(), im.height(), im.channels());

	int sideW = int((width - 1.0) / 2.0);
	int sideH = int((height - 1.0) / 2.0);
	float sumVal;

	// for every pixel in the image
	for (int x = 0; x < imFilter.width(); x++)
		for (int y = 0; y < imFilter.height(); y++)
			for (int z = 0; z < imFilter.channels(); z++)
			{
				sumVal = 0.0;

				for (int xFilter = 0; xFilter < width; xFilter++)
					for (int yFilter = 0; yFilter < height; yFilter++)
					{

						// sum the image pixel values weighted by the filter
						// TODO: Give a hint in the document to use operator()(xFilter, yFilter)
						sumVal += operator()(xFilter, yFilter) * im.smartAccessor(x - xFilter + sideW, y - yFilter + sideH, z, clamp);
					}

				// assign the pixel the value from convolution
				imFilter(x, y, z) = sumVal;
			}

	return imFilter;
}

/**************************************************************
 //               DON'T EDIT BELOW THIS LINE                //
 *************************************************************/

// Create an image of 0's with a value of 1 in the middle. This function
// can be used to test that you have properly set the kernel values in your
// Filter object. Make sure to set k to be larger than the size of your kernel
FloatImage impulseImg(const int &k)
{
	// initlize a kxkx1 image of all 0's
	FloatImage impulse(k, k, 1);

	// set the center pixel to have intensity 1
	int center = floor(k / 2);
	impulse(center, center, 0) = 1;

	return impulse;
}

Filter::Filter(const vector<float> &fData, const int &fWidth, const int &fHeight)
{
	// TODO: check that width*height = length of filterVals and that width and height are odd

	kernel = fData;
	width = fWidth;
	height = fHeight;
}

Filter::Filter(const int &fWidth, const int &fHeight)
{
	width = fWidth;
	height = fHeight;
	kernel = std::vector<float>(width * height, 0);
}

const float &Filter::operator()(int x, int y) const
{
	if (x < 0 || x >= width)
		throw OutOfBoundsException();
	if (y < 0 || y >= height)
		throw OutOfBoundsException();

	return kernel[x + y * width];
}

float &Filter::operator()(int x, int y)
{
	if (x < 0 || x >= width)
		throw OutOfBoundsException();
	if (y < 0 || y >= height)
		throw OutOfBoundsException();

	return kernel[x + y * width];
}
Filter::~Filter() {}
