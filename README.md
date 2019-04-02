# aruco-pattern-write

Comments/Bugs/Problems: amy.tabb@ars.usda.gov

Basic code for generating a gird of aruco patterns.

~March 2019.  


# Underlying ideas; how and when to cite this work

This README file is to accompany code produced by Amy Tabb as a companion to a paper:
To be announced. Email for a draft.



Code release citation:

```
@misc{tabb2019code_using,
	title = {Code from: {Using} cameras for precise measurement of two-dimensional plant features},
	url = {https://data.nal.usda.gov/dataset/code-using-cameras-precise-measurement-two-dimensional-plant-features},
	urldate = {2019-04-01},
	publisher = {Ag Data Commons},
	author = {Tabb, Amy},
	year = {2019},
}
```


If you use this code in project that results in a publication, please cite at a minimum the paper above.  Otherwise, there are no restrictions in your use of this code.  However, no guarantees are expressed or implied.

# Related repositories

[amy-tabb/camera-as-scanner-data](https://github.com/amy-tabb/camera-as-scanner-data)

[amy-tabb/camera-as-scanner](https://github.com/amy-tabb/camera-as-scanner)

## Dependencies

This code uses the OpenCV 4.0 and OpenCV 4.0 extra modules and is written in C++.

### Tested operating system

This code has been tested on Ubuntu 16.04 and Ubuntu 18.04.  You are welcome to convert it to Windows, but I have not.  While OpenCV is available from distribution repositories, my long experience with it is has always been to build from the source to get the best results.

### OpenCV 4.0

To get the OpenCV4.0 extra modules to build, our experience is that you need to build *both* OpenCV and the extra modules together from source.  Instructions are here:

[OpenCV contributed modules on Github](https://github.com/opencv/opencv_contrib)

**Is getting all of this to work on your system too much of a pain and you are interested in a Docker release?  Let me know!  The squeaky wheel gets the grease.  Email above.**

## Building 

 1. Clone the git repository.  `cd` to the desired directory, then from the command line  -- ` git clone ` the repository name.
 
 2. If you use the [Eclipse CDT IDE](https://www.eclipse.org/cdt/), you can import the project directly after cloning.  You may have to alter the include directory for opencv4.  

2. **Note the OpenCV 4.x requires a >C++11 compiler!  I enable the GNU version with a flag : `-std=gnu++11.`  There other ways of doing this, as well as other choices of compiler.**

3. Required libraries are: opencv_core, opencv_imgcodecs, opencv_aruco.

4. Build the project.


## Running

The program takes one argument, which is the directory where the aruco pattern and its specification file is written.

Example:

```
./aruco-pattern-write-project /home/your-name/read_directory
```

An example is given in this repository via the directory `SampleInput`.

## Input   

Within the read directory used as an argument above, you will need to create a text file called `specification_file.txt`.  The example from directory `SampleInput` is given below. 

All of these values can be changed to reflect the wishes of the user, but highly suggest not changing `arc_code`, which is the code used to generate the dictionary of aruco patterns within OpenCV.  Also, `squaresX`*`squaresY`< `1000` given that the dictionary we have used is `cv::aruco::DICT_6X6_1000`.  

	```
	squaresX 12
	squaresY 15
	squareLength 200
	markerLength 100
	margins 100
	arc_code 11
	```

## Output

The output of the program is:

1. one image file, in .png format, of aruco patterns laid out in a grid pattern.  The name is determined based on the dimensions of the grid.  

2. one text file called `specification_file.txt`, that records the values used to create the grid.

## Altering

Alter the Create() function to output images with different dimensions as well as grids with more, or fewer aruco patterns.

	
