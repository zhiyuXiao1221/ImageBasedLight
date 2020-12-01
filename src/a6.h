#pragma once

#include "floatimage.h"

// A6 Functions: Scaling and Rotating an Image
FloatImage scaleNN(const FloatImage &im, float factor);
FloatImage scaleLin(const FloatImage &im, float factor);
float interpolateLin6(const FloatImage &im, float x, float y, int z, bool clamp=false);
float interpolateCubic(const FloatImage &im, float x, float y, int z, bool clamp=false);
FloatImage rotate(const FloatImage &im, float theta);

std::vector<FloatImage> scaleU(const FloatImage &im, float factor);