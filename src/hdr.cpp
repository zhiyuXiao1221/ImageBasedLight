// hdr.cpp
// Assignment 5

#include "hdr.h"
#include "filtering.h"
#include "a2.h"
#include "utils.h"
#include <math.h>
#include <algorithm>

using namespace std;

/**************************************************************
 //                       HDR MERGING                        //
 *************************************************************/

// Generate a weight image that indicates which pixels are good to use in hdr
FloatImage computeWeight(const FloatImage &im, float epsilonMini, float epsilonMaxi)
{
	// should return an image with pixel values = 1 if im(...) falls
	// in the range [epsilonMini, epsilonMaxi]
	FloatImage output(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.width(); i++)
	{
		for (int j = 0; j < im.height(); j++)
		{
			for (int c = 0; c < im.channels(); c++)
			{
				if (epsilonMaxi < im(i, j, c) || im(i, j, c) < epsilonMini)
					output(i, j, c) = 0.0;
				else
					output(i, j, c) = 1.0;
			}
		}
	}
	return output;
}

// Compute the multiplication factor between a pair of images
float computeFactor(const FloatImage &im1, const FloatImage &w1, const FloatImage &im2, const FloatImage &w2)
{
	vector<float> ratio;
	int median = 0.0;
	// populate list with ratios of im1 and im2 for all valid pixels
	for (int i = 0; i < im1.width(); i++)
	{
		for (int j = 0; j < im1.height(); j++)
		{
			for (int c = 0; c < im1.channels(); c++)
			{
				//use valid pixels
				if (w1(i, j, c) == 1.0 && w2(i, j, c) == 1.0)
				{
					float current_ratio = im2(i, j, c) / (im1(i, j, c) + 1e-10);
					ratio.push_back(current_ratio);
				}
			}
		}
	}
	// return midpoint of ratio list
	if (ratio.size() % 2 == 0)
	{
		median = (ratio.size() / 2 + ratio.size() / 2 - 1) / 2;
	}
	else
	{
		median = ratio.size() / 2;
	}
	std::nth_element(ratio.begin(), ratio.begin() + median, ratio.end());

	return ratio[median];
}

// Merge images to make a single hdr image
FloatImage makeHDR(vector<FloatImage> &imSeq, float epsilonMini, float epsilonMaxi)
{
	vector<float> factors;
	float currentFactor;
	float lastFactor;
	vector<FloatImage> weights;
	FloatImage output(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
	for (int i = 0; i < ((int)imSeq.size()); i++)
	{
		// invert gamma correction
		imSeq[i] = changeGamma(imSeq[i], 1.0 / 2.2, 1.0f);
		//two special cases
		if (i == 0)
			weights.push_back(computeWeight(imSeq[i], epsilonMini, 1.0));
		else if (i == (int)imSeq.size() - 1)
			weights.push_back(computeWeight(imSeq[i], 0.0, epsilonMaxi));
		//other cases
		else
			weights.push_back(computeWeight(imSeq[i], epsilonMini, epsilonMaxi));
	}
	//get the first and second factor
	factors.push_back(1.0);
	float secondFactor = computeFactor(imSeq[0], weights[0], imSeq[1], weights[1]);
	factors.push_back(secondFactor);
	lastFactor = secondFactor;
	//get other factors
	for (int i = 1; i < (int)imSeq.size() - 1; i++)
	{
		currentFactor = computeFactor(imSeq[i], weights[i], imSeq[i + 1], weights[i + 1]);
		factors.push_back(currentFactor * lastFactor);
		lastFactor *= currentFactor;
	}

	for (int i = 0; i < output.width(); i++)
	{
		for (int j = 0; j < output.height(); j++)
		{
			for (int c = 0; c < output.channels(); c++)
			{
				float sum = 0.0;
				float weightSum = 0.0;
				for (int k = 0; k < (int)imSeq.size(); k++)
				{
					weightSum += weights[k](i, j, c);
					sum += (weights[k](i, j, c) * (imSeq[k](i, j, c)) / factors[k]);
				}
				if (weightSum == 0)
					output(i, j, c) = 0;
				else
					output(i, j, c) = sum / weightSum;
			}
		}
	}
	return output; // change this
}

/**************************************************************
 //                      TONE MAPPING                        //
 *************************************************************/

// Tone map an hdr image
FloatImage toneMap(const FloatImage &im, float targetBase, float detailAmp, bool useBila, float sigmaRange)
{
	// UNCOMMENT THIS LINE!
	// add gamma correction back into the image right before returning
	FloatImage output(im.width(), im.height(), im.channels());
	FloatImage lumi = lumiChromi(im)[0];
	FloatImage color = lumiChromi(im)[1];
	FloatImage lumi_log(im.width(), im.height(), im.channels());
	FloatImage lumi_base(im.width(), im.height(), im.channels());
	lumi_log = log10FloatImage(lumi);

	float sigmaDomain = im.width() > im.height() ? im.width() / 50 : im.height() / 50;
	if (useBila)
	{
		lumi_base = bilateral(lumi_log, sigmaRange, sigmaDomain);
	}
	else
	{
		lumi_base = gaussianBlur_2D(lumi_log, sigmaDomain);
	}

	FloatImage lumi_detail = lumi_log - lumi_base;
	//Shrink the dynamic range in lumi_base
	float lumi_range = lumi_base.max() - lumi_base.min();
	float k = log10(targetBase) / lumi_range;
	//get the output_log
	lumi = detailAmp * lumi_detail + k * (lumi_base - lumi_base.max());
	//get output
	lumi = exp10FloatImage(lumi);
	output = lumiChromi2rgb(lumi, color);
	//decode
	output = changeGamma(output, 1.0f, 1 / 2.2);
	return output;
}

// Tone Mapping Helpful Functions

// image --> log10FloatImage
FloatImage log10FloatImage(const FloatImage &im)
{
	// Taking a linear image im, transform to log10 scale.
	// To avoid infinity issues, make any 0-valued pixel be equal the the minimum
	// non-zero value. See image_minnonzero(im).
	FloatImage output(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.width(); i++)
	{
		for (int j = 0; j < im.height(); j++)
		{
			for (int c = 0; c < im.channels(); c++)
			{
				if (im(i, j, c) == 0.0)
				{
					output(i, j, c) = image_minnonzero(im);
				}
				output(i, j, c) = log10(im(i, j, c));
			}
		}
	}
	return output; // change this
}

// FloatImage --> 10^FloatImage
FloatImage exp10FloatImage(const FloatImage &im)
{
	// take an image in log10 domain and transform it back to linear domain.
	// see pow(a, b)
	FloatImage output(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.width(); i++)
	{
		for (int j = 0; j < im.height(); j++)
		{
			for (int c = 0; c < im.channels(); c++)
			{
				output(i, j, c) = pow(10, im(i, j, c));
			}
		}
	}
	return output;
}

// min non-zero pixel value of image
float image_minnonzero(const FloatImage &im)
{
	float min = 100.0;
	for (int i = 0; i < im.width(); i++)
	{
		for (int j = 0; j < im.height(); j++)
		{
			for (int c = 0; c < im.channels(); c++)
			{
				if (im(i, j, c) != 0.0 && min > im(i, j, c))
				{
					min = im(i, j, c);
				}
			}
		}
	}
	return min;
}
