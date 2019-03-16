# aruco-pattern-write

Comments/Bugs/Problems: amy.tabb@ars.usda.gov

Basic code for generating a gird of aruco patterns.

~March 2019.  



# Underlying ideas; how and when to cite this work

This README file is to accompany code produced by Amy Tabb as a companion to a paper:
**TODO**



If you use this code in project that results in a publication, please cite at a minimum the paper above.  Otherwise, there are no restrictions in your use of this code.  However, no guarantees are expressed or implied.

# Building

This README covers instructions for building the code.

## Dependencies

This code uses the OpenCV 4.0 and OpenCV 4.0 extra modules and is written in C++.

To get the OpenCV4.0 extra modules to build, our experience is that you need to build *both* OpenCV and the extra modules together from source.  Instructions are here:

[OpenCV contributed modules on Github](https://github.com/opencv/opencv_contrib)

This code has been tested on Ubuntu 16.04 and Ubuntu 18.04.  You are welcome to convert it to Windows, but I have not.  While OpenCV is available from distribution repositories, my long experience with it is has always been to build from the source to get the best results.

**Is getting all of this to work on your system too much of a pain and you are interested in a Docker release?  Let me know!  The squeaky wheel gets the grease.  Email above.**

## Building 


 1. Clone the git repository.  `cd` to the desired directory, then from the command line  -- ` git clone ` the repository name.
 
 2. If you use the [Eclipse CDT IDE](https://www.eclipse.org/cdt/), you can import the project directly after cloning.  You may have to alter the include directory for opencv4.  

2. **Note the OpenCV 4.x requires a >C++11 compiler!  I enable the GNU version with a flag : `-std=gnu++111.  There other ways of doing this, as well as other choices of compiler.**

3. Required libraries are: opencv_core, opencv_imgcodecs, opencv_aruco.


## Running

The program takes one argument, which is the directory where the aruco pattern and its specification file is written.

## Output

The output of the program is:

1. one image file, in .png format, of aruco patterns laid out in a grid pattern.  The name is determined based on the dimensions of the grid.  

2. one text file called specification_file.txt, that records the values used to create the grid.

## Altering

Alter the Create() function to output images with different dimensions as well as grids with more, or fewer aruco patterns.


	