#include "hdr.h"
#include "panoramicTrans.h"
#include "a2.h"
#include "a6.h"
#include "utils.h"
#include "filtering.h"
#include <sstream>
#include <string>

using namespace std;

//testPanoramicTrans
void testPanoramicTrans()
{
	const FloatImage im(DATA_DIR "/input/final_project/indoor/indoor_sphere/4.png");
	FloatImage output = sphere2Latlong(im);
	output.write(DATA_DIR "/output/test2Result.png");
}

void testCompositing()
{
	const FloatImage scene(DATA_DIR "/input/final_project/differentialRendering/scene.png");
	const FloatImage mask(DATA_DIR "/input/final_project/differentialRendering/mask.png");
	const FloatImage withObj(DATA_DIR "/input/final_project/differentialRendering/withObjects.png");
	const FloatImage withoutObj(DATA_DIR "/input/final_project/differentialRendering/withoutObjects.png");

	// create composite using default c value
	// c value used for changing shadows and reflection values
	FloatImage output = composite(scene, mask, withObj, withoutObj, 1.0);
	output.write(DATA_DIR "/output/testComposite-default.png");

	// create composite images with varying c values
	// images are generated in 0.5 c value step increments
	int count = 0;
	for (float i = 0.0; i < 5.0; i += 0.5)
	{
		ostringstream ss;
		ss << DATA_DIR "/output/compositeResults/testComposite-" << count << ".png";
		string filename = ss.str();
		output = composite(scene, mask, withObj, withoutObj, i);
		output.write(filename);
		count++;
	}
}
void testFastBalateral()
{
	FloatImage im(DATA_DIR "/input/final_project/Cambridge2.png");
	Timer timer;
	timer.reset();
	//Different downsample size
	//The bigger the sample size, the faster the time, but the courser the image

	// When sampleD =1, took 16.88900 seconds
	// When sampleD =2, took 4.30800 seconds
	// When sampleD =4,  took 1.20200 seconds
	// When sampleD =8,  took 0.48600 seconds

	FloatImage fastBilateral1 = fastBilateral(im, 5, 0.1, 16.0, 1.0, 0.05);
	printf("When sampleD =1, took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral2 = fastBilateral(im, 5, 0.1, 16.0, 2.0, 0.05);
	printf("When sampleD =2, took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral3 = fastBilateral(im, 5, 0.1, 16.0, 4.0, 0.05);
	printf("When sampleD =4,  took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral4 = fastBilateral(im, 5, 0.1, 16.0, 8.0, 0.05);
	printf("When sampleD =8,  took %3.5f seconds\n", timer.elapsed() / 1000.f);

	timer.reset();
	FloatImage fastBilateral5 = fastBilateral(im, 5, 0.1, 16.0, 16.0, 0.05);
	printf("When sampleD =16, took %3.5f seconds\n", timer.elapsed() / 1000.f);

	fastBilateral1.write(DATA_DIR "/output/fastBilateralResults/fastBilateral1.png");
	fastBilateral2.write(DATA_DIR "/output/fastBilateralResults/fastBilateral2.png");
	fastBilateral3.write(DATA_DIR "/output/fastBilateralResults/fastBilateral3.png");
	fastBilateral4.write(DATA_DIR "/output/fastBilateralResults/fastBilateral4.png");
	fastBilateral5.write(DATA_DIR "/output/fastBilateralResults/fastBilateral5.png");
}

// a function to test scaling
void testScaling(){
	// load in the image and print out original size
	float fact = 0.10;
	const FloatImage bostonim(DATA_DIR "/input/final_project/Cambridge2.png");
	printf("Boston image is %dx%dx%d\n", bostonim.width(), bostonim.height(), bostonim.channels());

	// scale using NN interpolation and print the size of the new image
	FloatImage scaledNN = scaleNN(bostonim, fact);
	scaledNN.write(DATA_DIR "/output/Cambridge2-scaled-NN.png");
	printf("Scaled-NN image is %dx%dx%d\n", scaledNN.width(), scaledNN.height(), scaledNN.channels());
}

// test CRF calculation
void testCRF(){
	int nImages = 8;
	float smooth = 50;
	vector<FloatImage> imSeq;
	for (int i = 0; i < nImages; i++)
	{
		ostringstream ss;
		// ss << DATA_DIR "/input/final_project/indoor/" << i << ".jpg";
		ss << "/home/sb/Desktop/ph5/scene-" << i << ".png";
		string filename = ss.str();
		imSeq.push_back(FloatImage(filename));
	}

	vector<float> exposures{1/80., 1/50., 1/15., 1/6., 0.3, 0.8, 2, 4};
	// calibrateCRF(imSeq, exposures, smooth);

	vector<FloatImage> sampled = sampleFromHist(imSeq);
	// for (int i = 0; i < nImages; i++){
	// 	sampled[i].write(DATA_DIR "/output/hist-" + to_string(i) + ".png");
	// }
	calibrateCRF(sampled, exposures, smooth);

}

int main()
{
	//testPanoramicTrans();
	//testCompositing();
	// testFastBalateral();
	// try { testScaling(); }        catch(...) {cout << "testScaling Failed!" << endl;}
	testCRF();

}
