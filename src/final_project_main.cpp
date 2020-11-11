#include "hdr.h"
#include "panoramicTrans.h"
#include "a2.h"
#include <sstream>

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
	for (float i = 0.0; i < 5.0; i += 0.5){
		ostringstream ss;
		ss << DATA_DIR "/output/compositeResults/testComposite-" << count << ".png";
		string filename = ss.str();
		output = composite(scene, mask, withObj, withoutObj, i);
		output.write(filename);
		count++;
	}


}

int main()
{
	testPanoramicTrans();
	testCompositing();
}
