#include "corner_detect.h"
#include "a2.h"
#include "filtering.h"
#include <math.h>

using namespace std;

// corners are detected in scene and returned as a list of {x, y} coordinate pairs
vector<vector<int>> detectCorners(const FloatImage &im, float threshold, int windowSize, bool useGaussianBlur, bool clamp) {
	vector<vector<int>> corners;
	vector<int> corner;

	// image is grayscaled
	FloatImage scene = color2gray(im);

	if(useGaussianBlur) {
		float sigma = 3.0;
		scene = gaussianBlur_seperable(scene, sigma);
	}

	// sobel filter of scene for x and y components calculated 
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

	FloatImage imSobelX = sobelX.Convolve(scene, clamp);
	FloatImage imSobelY = sobelY.Convolve(scene, clamp);

	// Harris corner detection constant referenced from source below
	// https://docs.opencv.org/master/dc/d0d/tutorial_py_features_harris.html
	const float harris_corner_const = 0.04;

	// calculate rolling window of x and y values
	// calculated values that pass threshold are recognized as a corner
	float i_xx, i_xy, i_yy, i_xx_yy_sum, shift_intensity;
	for (int x = 0; x < scene.width(); x++)
	{
		for (int y = 0; y < scene.height(); y++)
		{
			i_xx = 0;
			i_xy = 0;
			i_yy = 0;

			for(int shiftX = -windowSize/2; shiftX < windowSize/2; shiftX++){
				for(int shiftY = -windowSize/2; shiftY < windowSize/2; shiftY++){
					i_xx += pow(imSobelX.smartAccessor(x + shiftX, y + shiftY, 0, clamp), 2);
					i_xy += imSobelX.smartAccessor(x + shiftX, y + shiftY, 0, clamp) * imSobelY.smartAccessor(x + shiftX, y + shiftY, 0, clamp);
					i_yy += pow(imSobelY.smartAccessor(x + shiftX, y + shiftY, 0, clamp), 2);
				}
			}
			float determinant = i_xx * i_yy - pow(i_xy, 2);
			i_xx_yy_sum = i_xx + i_yy;
			shift_intensity = determinant - harris_corner_const * pow(i_xx_yy_sum, 2);

			if(shift_intensity >= threshold){
				corner.push_back(x);
				corner.push_back(y);
				corners.push_back(corner);
				corner.clear();
			}
		}
	}

	return corners;
}

// scene and list of corner coordinated are passed 
// returns an image with every detected corner highlighted in red
FloatImage highlightCorners(const FloatImage &im, vector<vector<int>> points){
	int markSize = 5;

	FloatImage output(im);

	for(auto& point : points){
		for (int x = -markSize/2; x < markSize/2; x++){
			for (int y = -markSize/2; y < markSize/2; y++){
				if(point[0] + x >= 0 && point[0] + x < im.width() && point[1] + y >= 0 && point[1] + y < im.height()){
					output(point[0] + x, point[1] + y, 0) = 255;
					output(point[0] + x, point[1] + y, 1) = 0;
					output(point[0] + x, point[1] + y, 2) = 0;
				}
			}
		}
	}
	return output;
}

// coordinates that are within a range are averaged into a single point
// returns a list of vector coordinates with consolidated clusters of points
vector<vector<int>> avgLocalPoints(const vector<vector<int>> points, int range){
	vector<vector<int>> avgPoints;

	bool localPointFound = false;

	for(auto& point : points){
		if(avgPoints.size() > 0){
			for(int i = 0; i < avgPoints.size(); i++){
				if(abs(avgPoints[i][0] - point[0]) <= range && abs(avgPoints[i][1] - point[1]) <= range){
					avgPoints[i][0] += point[0];
					avgPoints[i][0] /= 2;
					avgPoints[i][1] += point[1];
					avgPoints[i][1] /= 2;
					localPointFound = true;
				}
			}
			if(!localPointFound){
				avgPoints.push_back(point);
			}
			localPointFound = false;
		} else {
			avgPoints.push_back(point);
		}
	}

	return avgPoints;
}