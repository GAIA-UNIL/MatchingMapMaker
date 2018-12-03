# MatchingMapMaker

The MatchingMapMaker (M<sup>3</sup>) is a MATLAB function used to detect shifts between two images using convolutions.
The function will operate a weighted sum of absolute errors for each lag vector given as input.

The main code is written in C/C++, but the code is interfaced to MATLAB using mex.

The function **movsae2** use 4 parameters:
1. an odd size kernel, which can be full of one to compute the Sum of absolute error, or a normalized version to have the mean
2. reference image
3. compared image
4. a matrix, where each line represents a potential lag vector.

and return for each pixel:
1. the lag with the lowest error
2. the lowest error

### Usage in MATLAB:

`[lagIndex,quality]=movsae2(kernel,ref,shiftedImage,lagVector);`

### Installation

0. test with `mex -setup c++` that the mingw-w64 compiler is installed, in other cases install it from [here](https://www.mathworks.com/matlabcentral/fileexchange/52848-matlab-support-for-mingw-w64-c-c-compiler)
1. run `compile`
2. run `test` (to test if it works! )