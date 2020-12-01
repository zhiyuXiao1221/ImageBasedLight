## Final Project - Image-based Lighting
> Group 5   Zhiyu Xiao**|**Haowen Liu**|** Eammon Littler

### Motivation and Introduction 

Investigate the application of HDR in image-based lighting (IBL), which is a rendering technique in computer graphics that allows us to insert synthetic 3D models into photographs in a realistic way.

In our project, we took a series of photos with different exposures to synthesize the HDR image, and used the HDR image for image-based lighting. When synthesizing HDR images, we consider camera calibration, which makes the synthesized HDR images more accurate. Then, we used fast bilateral filter and Laplacian filter to tonemap the background images, and finally we rendered our model and scene in Maya to generate image illumination based images. We also tried rendering custom models in OpengL.

### How to Use

Before running this code, be sure to download XXX

### Core Functions in this Project

```c++
//Convert Mirrored sphere domain to equirectangular domain
FloatImage sphere2Latlong(const FloatImage &im, bool translation)
  
//Differential render compositing funtion
FloatImage composite(const FloatImage &scene, const FloatImage &mask, const FloatImage &withObj, const FloatImage &withoutObj, const float c)
  
//Fast bilateral function
FloatImage fastBilateral(const FloatImage &im, float sigmaRange, float sigmaDomain, int truncateDomain, float samplingD, float samplingR)
  
//Make image smoother when upsampleing in fast bilateral function
float trilinear_interpolation(const Array_3D<float> &array, float x, float y, float z)
```



### Problems and Challenges

1. When taking photos with different exposures, the mirror ball occupies too little space, so the mirror ball looks blurry. But the good news is that the rendered results look good.
2. In the fast bilateral filtering method, this method may fail to run when the image is very large and does not downsample.





### Task Division

1. **Recovering HDR Image**

   Shoot and crop scenes with different exposures[Zhiyu Xiao]

   Response function estimation[Haowen Liu]

   LDR image Merging [Haowen Liu]

2. **Panoramic Transformation**

   Converted from mirror ball to latitude-longitude[ Zhiyu Xiao]

3. **Tone Mapping**

   Fast bilateral filter[Zhiyu Xiao]

   Local laplacian filters [Haowen Liu]

4. **Modeling and Rendering**

   Model and render the scene in Maya [Zhiyu Xiao]

   Use rendered images to perform "differential render" compositing [Eammon Littler]

   OpenGL based local rendering [Eammon Littler]

### Reference

1. [Programming Project #4: Image-Based Lighting CS498](https://courses.engr.illinois.edu/cs498dh3/fa2014/projects/ibl/ComputationalPhotography_ProjectIBL.html)
2. [Advanced High Dynamic Range Imaging Book](http://advancedhdrbook.com/)
3. [A Fast Approximation of the Bilateral Filter using a Signal Processing Approach](http://people.csail.mit.edu/sparis/publi/2006/tr/Paris_06_Fast_Bilateral_Filter_MIT_TR.pdf)
4. [Rendering Synthetic Objects into Real Scenes, by Paul Debevec, Siggraph 1998](http://www.pauldebevec.com/Research/IBL/)
5. [Recovering High Dynamic Range Radiance Maps from Photographs, by Paul Debevec, Siggraph 1997](http://www.pauldebevec.com/Research/HDR/)
6. [Local Laplacian Filters: Edge-aware Image Processing with a Laplacian Pyramid](https://people.csail.mit.edu/sparis/publi/2015/cacm/Paris_15_Local_Laplacian_Filters.pdf)
7. [Harris Corner Detection](https://www.pauldebevec.com/Research/IBL/debevec-siggraph98.pdff)





