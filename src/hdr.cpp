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
#ifdef ANY_ASSIGNMENT
	return FloatImage(); // change this
#else
	FloatImage weight(im.width(), im.height(), im.channels());
	weight.clear(0.f);

	for (int x = 0; x < im.width(); x++)
		for (int y = 0; y < im.height(); y++)
			for (int z = 0; z < im.channels(); z++)
				if ((im(x, y, z) >= epsilonMini) && (im(x, y, z) <= epsilonMaxi))
					weight(x, y, z) = 1;
	return weight;
#endif
}

// Compute the multiplication factor between a pair of images
float computeFactor(const FloatImage &im1, const FloatImage &w1, const FloatImage &im2, const FloatImage &w2)
{
	vector<float> ratio;
	float ratioVal;

	FloatImage im1tmp = im1 + 1e-10;

	for (int x = 0; x < im1.width(); x++)
		for (int y = 0; y < im1.height(); y++)
			for (int z = 0; z < im1.channels(); z++)
				if (w1(x, y, z) > 0.1 && w2(x, y, z) > 0.1)
				{
					ratioVal = im2(x, y, z) / im1tmp(x, y, z);
					ratio.push_back(ratioVal);
				}

	sort(ratio.begin(), ratio.end());
	int midpoint = floor(float(ratio.size()) / 2.0);

	float factor = ratio[midpoint];

	return factor;
}

// Merge images to make a single hdr image
FloatImage makeHDR(vector<FloatImage> &imSeq, float epsilonMini, float epsilonMaxi)
{
	// invert gamma correction
	for (int i = 0; i < ((int)imSeq.size()); i++)
		imSeq[i] = changeGamma(imSeq[i], 1.0 / 2.2, 1.0f);

	vector<float> k(1, 1);
	FloatImage out(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
	FloatImage wTotal(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
	FloatImage w1 = wTotal;
	FloatImage w2 = wTotal;

	for (int i = 0; i < ((int)imSeq.size()) - 1; i++)
	{
		bool first = (i == 0);
		bool last = (i == ((int)imSeq.size()) - 2);
		w1 = computeWeight(imSeq[i], epsilonMini, first ? 1.0 : epsilonMaxi);
		w2 = computeWeight(imSeq[i + 1], last ? 0.0 : epsilonMini, epsilonMaxi);

		float factor = computeFactor(imSeq[i], w1, imSeq[i + 1], w2);
		k.push_back(factor * k[i]);

		out = out + (1.0 / k[i]) * w1 * imSeq[i];
		wTotal = wTotal + w1;

		if (last)
		{
			out = out + (1.0 / k[i + 1]) * w2 * imSeq[i + 1];
			wTotal = wTotal + w2;
		}
	}

	for (int x = 0; x < wTotal.width(); x++)
		for (int y = 0; y < wTotal.height(); y++)
			for (int z = 0; z < wTotal.channels(); z++)
				if (wTotal(x, y, z) < 1.0)
				{
					wTotal(x, y, z) = 1;
					out(x, y, z) = (1.0 / k[0]) * imSeq[0](x, y, z);
				}

	out = out / wTotal;
	return out;
}

/**************************************************************
 //                      TONE MAPPING                        //
 *************************************************************/

// Tone map an hdr image
FloatImage toneMap(const FloatImage &im, float targetBase, float detailAmp, bool useBila, bool useFast, float sigmaRange)
{
	// initiate baselumi. This is Unnecessary but otherwise have issues with scope?
	FloatImage baselumi(im.width(), im.height(), im.channels());

	// load in luminance and chrominance
	vector<FloatImage> lc = lumiChromi(im);

	// get the log of the luminance avoiding zero-valued pixels
	FloatImage loglumi = log10FloatImage(lc[0]);

	// compute the base
	int maxdim = max(loglumi.width(), loglumi.height());
	float sigma = ((float)maxdim) / 50;

	if (useBila)
	{ // bilateral filter
		if (useFast)
		{
			baselumi = fastBilateral(loglumi, 0.4, sigma, 3, 10, 0.05);
		}
		else
		{
			cout << "!" << endl;
			baselumi = bilateral(loglumi, sigmaRange, sigma);
		}
	}
	else
	{ // normal gaussian filter
		baselumi = gaussianBlur_seperable(loglumi, sigma);
	}

	// compute difference (detail)
	FloatImage detail = loglumi - baselumi;

	// range
	float largeRange = baselumi.max() - baselumi.min();
	float k = log10(targetBase) / largeRange;

	// final log luminance
	FloatImage outLog = detail * detailAmp + baselumi * k - baselumi.max() * k;

	// compute image
	FloatImage finallumi = exp10FloatImage(outLog);

	// finally, go back finaloglumi+chromi --> rgb
	FloatImage final = lumiChromi2rgb(finallumi, lc[1]);

	// gamma encode 1.0f --> 1/2.2
	final = changeGamma(final, 1.0f, 1 / 2.2);

	// return
	return final;
}
// Tone Mapping Helpful Functions
// image --> log10FloatImage
FloatImage log10FloatImage(const FloatImage &im)
{
	// Taking a linear image im, transform to log10 scale.
	// To avoid infinity issues, make any 0-valued pixel be equal the the minimum
	// non-zero value. See image_minnonzero(im).
	float minf = image_minnonzero(im);
	FloatImage loglumi(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.size(); i++)
		loglumi(i) = log10(max(minf, im(i)));
	return loglumi;
}

// FloatImage --> 10^FloatImage
FloatImage exp10FloatImage(const FloatImage &im)
{
	// take an image in log10 domain and transform it back to linear domain.
	// see pow(a, b)
	FloatImage exp10(im.width(), im.height(), im.channels());
	for (int i = 0; i < im.size(); i++)
		exp10(i) = pow(10.0, im(i));
	return exp10;
}

// min non-zero pixel value of image
float image_minnonzero(const FloatImage &im)
{
	float minf = 2;
	long imsize = im.size();
	for (int i = 0; i < imsize; i++)
		if (im(i) > 0)
			minf = min(minf, im(i));

	return minf;
}
