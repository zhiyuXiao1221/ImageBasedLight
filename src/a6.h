#pragma once

#include "floatimage.h"

// A6 Functions: Scaling and Rotating an Image
FloatImage scaleNN(const FloatImage &im, float factor);
FloatImage scaleLin(const FloatImage &im, float factor);
float interpolateLin(const FloatImage &im, float x, float y, int z, bool clamp=false);
FloatImage rotate(const FloatImage &im, float theta);