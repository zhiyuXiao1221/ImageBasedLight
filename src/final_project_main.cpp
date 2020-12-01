#include "hdr.h"
#include "panoramicTrans.h"
#include "a2.h"
#include "utils.h"
#include "filtering.h"
#include "corner_detect.h"
#include <sstream>

using namespace std;

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
	const FloatImage scene(DATA_DIR "/input/final_project/differentialRendering/3/scene.png");
	const FloatImage mask(DATA_DIR "/input/final_project/differentialRendering/3/mask.png");
	const FloatImage withObj(DATA_DIR "/input/final_project/differentialRendering/3/withObjects.png");
	const FloatImage withoutObj(DATA_DIR "/input/final_project/differentialRendering/3/withoutObjects.png");

	// create composite using default c value
	// c value used for changing shadows and reflection values
	FloatImage output = composite(scene, mask, withObj, withoutObj, 1.0);
	output.write(DATA_DIR "/output/Composite-default-2.png");

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
		ss << DATA_DIR "/input/final_project/indoor/scene1/Angle2/sphere_" << i << ".jpg";
		string filename = ss.str();
		imSeq.push_back(FloatImage(filename));
	}

	// generate an hdr image
	FloatImage hdr = makeHDR(imSeq, 0.1, 0.9);
	// save out HDR image
	hdr.write(DATA_DIR "/output/room3-out.hdr");
}

void testCornerDetection() {
	// image of real life scene 
	const FloatImage scene(DATA_DIR "/input/final_project/corner_detection/scene_with_checkerboard.jpg");

	// settings for tuning corner detectors
	int windowSize = 2;
	float threshold = 100;
	vector<vector<int>> corners = detectCorners(scene, threshold, windowSize, false, true);

	// mark in red where corners are detected
	FloatImage output = highlightCorners(scene, corners);
	output.write(DATA_DIR "/output/corner_detection/corner_detect.png");
	
	// average nearby detected corners within a given range
	int range = 10;
	vector<vector<int>> avgCorners = avgLocalPoints(corners, range);
	FloatImage outputAvgLocal = highlightCorners(scene, avgCorners);
	outputAvgLocal.write(DATA_DIR "/output/corner_detection/corner_detect_avg_local_pts.png");
}

int main()
{
	//testPanoramicTrans();
	// testCompositing();
	// testFastBilateral();
	//CompareTwoBilateral();
	//testMakeNaiveHdr_Room();
	testCornerDetection();
}
