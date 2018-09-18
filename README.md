# fastGMatch
  * The programs are licenced by GPL v3. 
  * Fast image matching method using Gaussian smoothing by sliding discrete Fourier 
  * Sample images (Graffiti and Boat) are contained for check the program.
  
## Programs
  * C programs for the matching program
  * M-files to check the former's results.

## Compile of matching program
  `cd src`
  `make`

## Usage of the matching program
 1. `cd src`  
 1. `./fastGMatch command imageName templateName`
  * command = 0 : Search a part that is similar to template data (templateName) from input image (imageName)
  * command = 1 : Make and output template data (templateName) from input image (imageName)
  * Folders for input images and output images are specified by macro variables (WORKDIR and WORKOUTDIR, respectively) in fastGMatch.h

## Example
  * Make template  
  `./fastGMatch 1 img1 img1.template`
  * Matching  
  `./fastGMatch 0 img1 img1.template`

## Files for for the matching program
  * `fastGMatch.c` : Main program and function for matching.
  * `fastGMatch.h` : Parameters to use this program and header file for fastGMatch.c.
  * `calFeature.c` : Programs to calculate features.
  * `calFeature.h` : Header file for calFeature.c.
  * `utilities.c` : Utility functions (File input/output). 
  * `utilities.h` : Header file for utilities.c.

## Usage of the checking result program
 1.  `cd matlab`
 1.  `matlab &`
 1. In the matlab command shell type  
  `mainMatchCheck`
 1.  Then, it show the template image with center, scale, and main direction and several best matched pairs of the (transformed) template image and input image with centers, scales, and main directions. The images are also stored in a folder.

## Files for for the checking result program
  * `mainMatchCheck.m` : Main program. It include parameters to specify image and the number of output
  * `readFeatures.m` : Read features from a file describing features in an image.
  * `drawMatch.m` : To draw a circle and line in an image to show scale and main direction, respectively

## Parameters of mainMatchCheck.m
  * workDir : images for input
  * workOutDir : images and data for input and output
  * inpImgName : Name of input (ex. `img1`)
  * tmplImgName : Name of template (ex. `img1` ('template' is added.)
  * stOrd : The order to start output the matching result.
  * enOrd : The order to end output the matching result. (If we set `stOrd = 1` and `edOrd = 4` for example, then from the best matching result to the 4-th best matching result are shown.)