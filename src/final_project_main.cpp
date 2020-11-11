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
int main()
{
	testPanoramicTrans();
}
