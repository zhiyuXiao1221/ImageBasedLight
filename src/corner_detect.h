#include "floatimage.h"

std::vector<std::vector<int>> detectCorners(const FloatImage &im, float threshold, int windowSize, bool useGaussianBlur, bool clamp);
FloatImage highlightCorners(const FloatImage &im, std::vector<std::vector<int>> points);
std::vector<std::vector<int>> avgLocalPoints(const std::vector<std::vector<int>> points, int range);