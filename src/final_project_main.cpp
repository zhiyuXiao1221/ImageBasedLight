#include "hdr.h"
#include "panoramicTrans.h"
#include "a2.h"
#include "a6.h"
#include "utils.h"
#include "filtering.h"
#include <sstream>

using namespace std;
using namespace Eigen;

//testPanoramicTrans
void testPanoramicTrans()
{
	const FloatImage im(DATA_DIR "/input/final_project/hdr/sphere_hdr_robertson.hdr");
	FloatImage output1 = sphere2Latlong(im);
	output1.write(DATA_DIR "/output/HDR_Latlong/robertson.hdr");
	FloatImage output2 = sphere2Latlong(im, true);
	output2.write(DATA_DIR "/output/HDR_Latlong/robertson_trans.hdr");
}

void testCompositing()
{
	const FloatImage scene(DATA_DIR "/input/final_project/differentialRendering/1/bg.png");
	const FloatImage mask(DATA_DIR "/input/final_project/differentialRendering/1/mask.png");
	const FloatImage withObj(DATA_DIR "/input/final_project/differentialRendering/1/withObjects.png");
	const FloatImage withoutObj(DATA_DIR "/input/final_project/differentialRendering/1/withoughtObjts.png");

	// create composite using default c value
	// c value used for changing shadows and reflection values
	FloatImage output = composite(scene, mask, withObj, withoutObj, 1.0);
	output.write(DATA_DIR "/output/Composite-default.png");

	// create composite images with varying c values
	// images are generated in 0.5 c value step increments
	int count = 0;
	for (float i = 0.0; i < 5.0; i += 0.5)
	{
		ostringstream ss;
		ss << DATA_DIR "/output/compositeResults/Composite-" << count << ".png";
		string filename = ss.str();
		output = composite(scene, mask, withObj, withoutObj, i);
		output.write(filename);
		count++;
	}
}
void testFastBilateral()
{
	FloatImage im(DATA_DIR "/input/final_project/lens-3-med.png");
	Timer timer;
	timer.reset();
	//Different downsample size
	//The bigger the sample size, the faster the time, but the courser the image

	// When sampleD =1, took 2.80500 seconds
	// When sampleD =2, took 0.98400 seconds
	// When sampleD =4,  took 0.57900 seconds
	// When sampleD =8,  took 0.10000 seconds
	// When sampleD =16, took 0.06900 seconds

	FloatImage fastBilateral1 = fastBilateral(im, 0.1, 1.0, 3, 1.0, 0.05);
	printf("When sampleD =1, took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral2 = fastBilateral(im, 0.1, 1.0, 3, 2.0, 0.05);
	printf("When sampleD =2, took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral3 = fastBilateral(im, 0.1, 1.0, 3, 4.0, 0.05);
	printf("When sampleD =4,  took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral4 = fastBilateral(im, 0.1, 1.0, 3, 8.0, 0.05);
	printf("When sampleD =8,  took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral5 = fastBilateral(im, 0.1, 1.0, 3, 16.0, 0.05);
	printf("When sampleD =16, took %3.5f seconds\n", timer.elapsed() / 1000.f);

	fastBilateral1.write(DATA_DIR "/output/fastBilateralResults/fastBilateral1.png");
	fastBilateral2.write(DATA_DIR "/output/fastBilateralResults/fastBilateral2.png");
	fastBilateral3.write(DATA_DIR "/output/fastBilateralResults/fastBilateral3.png");
	fastBilateral4.write(DATA_DIR "/output/fastBilateralResults/fastBilateral4.png");
	fastBilateral5.write(DATA_DIR "/output/fastBilateralResults/fastBilateral5.png");
}
void CompareTwoBilateral()
{
	// 	Fast bilateral took 0.84500 seconds
	// Bilateral took 2.76100 seconds
	Timer timer;
	timer.reset();
	FloatImage im(DATA_DIR "/input/final_project/lens-3-med.png");
	FloatImage fastBilateral1 = fastBilateral(im, 0.1, 1.0, 3, 2.0, 0.05);
	printf("Fast bilateral took %3.5f seconds\n", timer.elapsed() / 1000.f);
	fastBilateral1.write(DATA_DIR "/output/fastBilateralResults/Compare_fastBilateral_noDownsample.png");
	// Perform bilaterial filtering on an RGB image
	timer.reset();
	FloatImage rgbBilatIm = bilateral(im);
	rgbBilatIm.write(DATA_DIR "/output/fastBilateralResults/Compare_Bilateral.png");
	printf("Bilateral took %3.5f seconds\n", timer.elapsed() / 1000.f);
}
void testMakeNaiveHdr_Room()
{
	vector<string> filenames;
	int nImages = 8;

	// load an image sequence
	vector<FloatImage> imSeq;
	for (int i = 1; i <= nImages; i++)
	{
		ostringstream ss;
		ss << DATA_DIR "/input/final_project/indoor/imgs/room-" << i << ".png";
		string filename = ss.str();
		imSeq.push_back(FloatImage(filename));
	}

	// generate an hdr image
	FloatImage hdr = makeHDR(imSeq, 0.1, 0.99);
	// save out HDR image
	hdr.write(DATA_DIR "/output/room3-out.hdr");
	// tone map with fast bilaterial
	Timer timer;
	timer.reset();
	//toneMap(hdr, 20, 1.5, true, true, 0.1).write(DATA_DIR "/output/room-fastbilateral.png");
	toneMap(hdr, 100, 1.5, true, true, 0.05).write(DATA_DIR "/output/room-fastbilateral.png");
	printf("Bilateral took %3.5f seconds\n", timer.elapsed() / 1000.f);
}

// test CRF calculation
void testCRF(){
	int nImages = 8;
	float smooth = 100;
	vector<FloatImage> imSeq;
	for (int i = 1; i <= nImages; i++)
	{
		ostringstream ss;
		ss << DATA_DIR "/input/final_project/indoor/scene1/Angle1/sphere_" << i << ".jpg";
		// ss << DATA_DIR "/input/final_project/indoor/scene1/Angle2/sphere_" << i << ".jpg";
		string filename = ss.str();
		imSeq.push_back(FloatImage(filename));
	}
	
	vector<float> exposures{1/80., 1/40., 1/20., 1/10., 1/5., 1/2., 0.77, 1.6};

	vector<FloatImage> sampled = sampleFromHist(imSeq);
	// vector<FloatImage> sampled = sampleDown(imSeq, 0.01);
	
	// for (int i = 0; i < nImages; i++){
	// 	sampled[i].write(DATA_DIR "/output/hist-" + to_string(i) + ".png");
	
	vector<VectorXf> crfs = calibrateCRF(sampled, exposures, smooth);
	FloatImage map = getRadiance(imSeq, exposures, crfs, true);
	map.write(DATA_DIR "/output/test-hdr.hdr");

	FloatImage scaledNN = scaleNN(map, 0.2);
	toneMap(scaledNN, 300, 1.5).write(DATA_DIR "/output/test-tonemap.png");
}

void testLaplacian(){
	FloatImage im(DATA_DIR "/input/final_project/Cambridge2.png");
	float levels = 3;
	float sigma = log(2.5);
	float alpha = 0.1;
	float beta = 0;
	int channels = 3;

	// vector<FloatImage> gPyramid = gaussPyramid(im, levels);
	// // vector<FloatImage> imSeq = laplacianPyramid(im, levels);
	// for (int i = 0; i < levels; i++){
	// 	gPyramid[i].write(DATA_DIR "/output/gauss-" + to_string(i) + ".png");
	// }

	// localLaplacianFilter(im, levels, sigma, alpha, beta, channels);
	
}

int main()
{
	//testPanoramicTrans();
	//testCompositing();
	//testFastBilateral();
	//CompareTwoBilateral();
<<<<<<< HEAD
	testMakeNaiveHdr_Room();
=======
	//testMakeNaiveHdr_Room();
	// testCRF();
	// testLaplacian();
>>>>>>> 2d5cd6f551b52d9d8ac838a737652ba99d557abd
}
