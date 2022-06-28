# Final Project - Image-based Lighting
The purpose of this project is to investigate the application of HDR imaging in image-based lighting (IBL) using environment maps, specifically in relighting virtual objects.
[The report pfd link is :](https://github.com/zhiyuXiao1221/ImageBasedLight/blob/main/CS273%20Final%20Project.pdf)
## Main Tasks

1. Captured LDR images of chrome sphere for IBL
2. Recover HDR image from LDR images with camera response curve calibration
3. Tone-mapping with fast bilateral filter and local Laplacian filter
4. Rendered IBL scene in Maya and custome OpenGL implementation
5. Marker based corner detection for rendering reference frame

## Prerequisite

Install the `Eigen` library.

## Functions Implemented

Some of the main functions we implemented on our own are listed here. There are some test cases available in `final_project_main.cpp`.

```c++
//panoramicTrans.cpp
//Convert Mirrored sphere domain to equirectangular domain
FloatImage sphere2Latlong(const FloatImage &im, bool translation)
  
//Differential render compositing funtion
FloatImage composite(const FloatImage &scene, const FloatImage &mask, const FloatImage &withObj, const FloatImage &withoutObj, const float c)

//filtering.cpp
//Fast bilateral function
FloatImage fastBilateral(const FloatImage &im, float sigmaRange, float sigmaDomain, int truncateDomain, float samplingD, float samplingR)
  
//Make image smoother when upsampleing in fast bilateral function
float trilinear_interpolation(const Array_3D<float> &array, float x, float y, float z)
    
//Local Laplacian Filter
FloatImage localLaplacianFilter(const FloatImage im, int levels, float sigma, float alpha, float beta, int channels)

//hdr.ccp
//Calculate CRF
vector<VectorXf> calibrateCRF(vector<FloatImage> &scaledSeq, vector<float> &exposures, float smooth)
    
//HDR merging
FloatImage getRadiance(vector<FloatImage> &imSeq, vector<float> exposures, vector<VectorXf> crfs, bool robertson)

// compositing
FloatImage composite(const FloatImage &scene, const FloatImage &mask, const FloatImage &withObj, const FloatImage &withoutObj, float c);

// corner detection
vector<vector<int>> detectCorners(const FloatImage &im, float threshold, int windowSize, bool useGaussianBlur, bool clamp);
```

## Problems and Challenges

1. Reflections on chrome ball are not very clear when cropped out as we took the LDR images from a far distance. The good news is that the IBL results seem to be unaffected.
2. Fast bilateral filtering may fail to run when the image is very large due to memory issues.
3. Local Laplacian filter may take a very long time and memory to run for larger images.
4. Edges of diffuse virtual objects composited into scene has color artifacts.
5. OpenGL libraries and code proved to be difficult to merge into main codebase. Additional repository with standalone rendering code can be found [here](https://github.com/ealitt/Scene-Rendering/tree/master).
6. Corner detection has not been tested with real life markers, only on ideal virtual markers. Further tuning for regional searches in an images would be necessary.

## Additional Materials 
Detailed description of some implementations as well as the validation of results can be found in the Google Doc link submitted. Task Division and a more comprehensive list of references can be found in the document as well.  

## References

1. [Programming Project #4: Image-Based Lighting CS498](https://courses.engr.illinois.edu/cs498dh3/fa2014/projects/ibl/ComputationalPhotography_ProjectIBL.html)
2. [Advanced High Dynamic Range Imaging Book](http://advancedhdrbook.com/)
3. [A Fast Approximation of the Bilateral Filter using a Signal Processing Approach](http://people.csail.mit.edu/sparis/publi/2006/tr/Paris_06_Fast_Bilateral_Filter_MIT_TR.pdf)
4. [Rendering Synthetic Objects into Real Scenes, by Paul Debevec, Siggraph 1998](http://www.pauldebevec.com/Research/IBL/)
5. [Recovering High Dynamic Range Radiance Maps from Photographs, by Paul Debevec, Siggraph 1997](http://www.pauldebevec.com/Research/HDR/)
6. [Local Laplacian Filters: Edge-aware Image Processing with a Laplacian Pyramid](https://people.csail.mit.edu/sparis/publi/2015/cacm/Paris_15_Local_Laplacian_Filters.pdf)
7. [Harris Corner Detection](https://docs.opencv.org/master/dc/d0d/tutorial_py_features_harris.html)
8. [Diffuse irradiance](https://learnopengl.com/PBR/IBL/Diffuse-irradiance)
9. [OpenGL Rendering](https://learnopengl.com/PBR/IBL/Specular-IBL)


