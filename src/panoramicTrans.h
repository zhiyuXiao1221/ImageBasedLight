#pragma once

#include "floatimage.h"
#include <iostream>
#include <math.h>

FloatImage computeXY(const FloatImage &im);
FloatImage sphere2Latlong(const FloatImage &im);
float interpolateLin(const FloatImage &im, float x, float y, int z, bool clamp);
