#include "hdr.h"
#include "panoramicTrans.h"
#include "a2.h"
#include <sstream>

using namespace std;

//testPanoramicTrans
void testPanoramicTrans()
{
	const FloatImage im(DATA_DIR "/input/final_project/test1.png");
	FloatImage output = sphere2Latlong(im);
	output.write(DATA_DIR "/output/test1.png");
}
int main()
{
	testPanoramicTrans();
}
