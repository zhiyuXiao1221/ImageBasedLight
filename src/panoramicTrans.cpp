#include "panoramicTrans.h"
#include "a2.h"
#include "utils.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
using namespace std;
FloatImage computeXY(const FloatImage &im)
{
    //get the square image
    int width = im.width();
    float radius = (float)width / 2.0;
    //Translate the position
    // int x0 = i - im.width() / 2;
    // int y0 = -j + im.height() / 2;
};
FloatImage sphere2Latlong(const FloatImage &im)
{
    float radius = (float)im.width() / 2.0;
    float factor = (float)im.width() / M_PI;
    FloatImage output(2 * im.width(), im.height(), im.channels());
    for (int i = 0; i < output.width(); i++)
    {
        for (int j = 0; j < output.height(); j++)
        {
            for (int c = 0; c < output.channels(); c++)
            {
                double phis = (float)j / factor;
                double thetas = (float)j / factor;
                //get the reflection vector
                float R_x = sin(thetas) * cos(phis);
                float R_y = sin(thetas) * sin(phis);
                float R_z = cos(thetas);
                //the view direction is(0,0,1)
                float V_x = 0.0;
                float V_y = 0.0;
                float V_z = 1.0;
                //get the normal vector
                float N_x = V_x + R_x;
                float N_y = V_y + R_y;
                float N_z = V_z + R_z;
                //standalize the normal vector
                float N_length = sqrt(N_x * N_x + N_y * N_y + N_z * N_z);
                N_x /= N_length;
                N_y /= N_length;
                N_z /= N_length;
                //now get the original x and original y
                output(i, j, c) = interpolateLin(im, N_x, N_y, c, true);
            }
        }
    }
    return output;
};
float interpolateLin(const FloatImage &im, float x, float y, int z, bool clamp)
{
    // Hint: use smartAccessor() to handle coordinates outside the image
    //get the 4 neighboring pixels
    float p_00 = im.smartAccessor(floor(x), floor(y), z, clamp);
    float p_01 = im.smartAccessor(ceil(x), floor(y), z, clamp);
    float p_10 = im.smartAccessor(floor(x), ceil(y), z, clamp);
    float p_11 = im.smartAccessor(ceil(x), ceil(y), z, clamp);
    float x_dis = x - floor(x);
    float y_dis = y - floor(y);
    //lerp
    float lerp_1 = p_00 + (p_01 - p_00) * x_dis;
    float lerp_2 = p_10 + (p_11 - p_10) * x_dis;
    float lerp_3 = lerp_1 + (lerp_2 - lerp_1) * y_dis;
    //return final float value
    return lerp_3;
}