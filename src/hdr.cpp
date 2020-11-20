// hdr.cpp
// Assignment 5


#include "hdr.h"
#include "filtering.h"
#include "a2.h"
#include "a6.h"
#include "utils.h"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <math.h>
#include <algorithm>
#include<fstream>

using namespace std;
using namespace Eigen;

/**************************************************************
 //                       HDR MERGING                        //
 *************************************************************/

// get intensity weight
float w(float z){
	return z > 127 ? 255 - z : z - 0;
}

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
void writeToCSVfile(string name, MatrixXf matrix)
{
    ofstream file(name.c_str());
    file << matrix.format(CSVFormat);
}

// calculate CRF
vector<VectorXf> calibrateCRF(vector<FloatImage> &imSeq, vector<float> &exposures, float smooth)
 {
	 vector<float> times;
	 vector<FloatImage> scaledSeq;
	 vector<VectorXf> crfs;
	 float fact = 0.01;
	 for (int i = 0; i < (int)imSeq.size(); i++){
		// downsample images
		scaledSeq.push_back(scaleNN(imSeq[i], fact));
		// scaledSeq.push_back(imSeq[i]);
		// log exposures
		times.push_back(log2(exposures[i]));
	}
	
	int pixelSize = scaledSeq[0].height() * scaledSeq[0].width();
	cout << "pixel size = " << pixelSize << endl;

	MatrixXf A(scaledSeq.size() * pixelSize + 256, 256 + pixelSize);
	VectorXf b(scaledSeq.size() * pixelSize + 256), x;
	cout << "size of A = " << A.rows() << ", " << A.cols() << endl;

	for (int z = 0; z < scaledSeq[0].channels(); z++){ // loop through all channels

		A.setZero();
		b.setZero();	
		int row = 0;

		for (int i = 0; i < (int)scaledSeq.size(); i++){
			for (int iy = 0; iy < scaledSeq[0].height(); iy++){
				for (int ix = 0; ix < scaledSeq[0].width(); ix++){
					
					// get pixel intensity
					float Z = floor(scaledSeq[i](ix, iy, z) * 255);

					A(row, Z) = w(Z);
					A(row, 256+iy*scaledSeq[0].width()+ix) = -w(Z);
					b[row] = w(Z) * times[i];

					row += 1;
					
				}		
			}
		}

		A(row, 127) = 1; // set f^-1(Z_mid) = 1
		row += 1;

		for (int p = 0; p < 255; p++){
			A(row, p) = smooth * w(p+1);
			A(row, p+1) = -2 * smooth * w(p+1);
			A(row, p+2) = smooth * w(p+1);
			// if (row >= 446 && row < 450) cout << row << " " << p+1  << " " << smooth * w(p+1) << endl;
			
			row += 1;
			
		}

		// writeToCSVfile(DATA_DIR "/output/A" + to_string(z) + ".csv", A);
		// writeToCSVfile(DATA_DIR "/output/b" + to_string(z) + ".csv", b);
		
		x = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b); // solve equation
		double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
   		cout << "The relative error is:\n" << relative_error << endl;
		
		x = x.head(256);
		crfs.push_back(x);
		// IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";");
		// cout << x.format(CommaInitFmt) << endl;
	}
	return crfs;
 }

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
	for (int i = 0; i < ((int) imSeq.size()); i++)
		imSeq[i] = changeGamma(imSeq[i], 1.0 / 2.2, 1.0f);

	vector<float> k(1, 1);
	FloatImage out(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
	FloatImage wTotal(imSeq[0].width(), imSeq[0].height(), imSeq[0].channels());
	FloatImage w1 = wTotal;
	FloatImage w2 = wTotal;

	for (int i = 0; i < ((int) imSeq.size()) - 1; i++)
	{
		bool first = (i == 0);
		bool last = (i == ((int) imSeq.size()) - 2);
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
FloatImage toneMap(const FloatImage &im, float targetBase, float detailAmp, bool useBila, float sigmaRange)
{
	// initiate baselumi. This is Unnecessary but otherwise have issues with scope?
	FloatImage baselumi(im.width(), im.height(), im.channels());

	// load in luminance and chrominance
	vector<FloatImage> lc = lumiChromi(im);

	// get the log of the luminance avoiding zero-valued pixels
	FloatImage loglumi = log10FloatImage(lc[0]);

	// compute the base
	int maxdim = max(loglumi.width(), loglumi.height());
	float sigma = ((float) maxdim) / 50;

	if (useBila)
	{   // bilateral filter
		baselumi = bilateral(loglumi, sigmaRange, sigma);
	}
	else
	{   // normal gaussian filter
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
