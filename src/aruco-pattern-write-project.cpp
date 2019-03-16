//============================================================================
// Name        : aruco-pattern-read-write-project.cpp
// Author      : Amy Tabb
// Version     :
// Copyright   : MIT
// Description :
//============================================================================


#include <sys/stat.h>
#include <sys/time.h>
#include <chrono>
#include <opencv2/core/core.hpp>
#include <opencv2/core/mat.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/highgui/highgui_c.h>
#include <opencv2/calib3d/calib3d.hpp>
#include <opencv2/aruco.hpp>
#include <opencv2/aruco/charuco.hpp>
#include <sys/stat.h>
#include <sys/time.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <chrono>


//#include "Includes.hpp"


using std::cout;
using std::endl;
using std::string;

using namespace cv;


template<class T>
T FromString(const std::string& s)
{
	std::istringstream stream (s);
	T t;
	stream >> t;
	return t;
}

template<class T>
string ToString(T arg)
{
	std::ostringstream s;

	s << arg;

	return s.str();

}


void EnsureDirHasTrailingBackslash(string& write_directory){
	int n_letters = write_directory.size();
	bool eval =  (write_directory[n_letters - 1] == '/');
	cout << "Last character compare " << write_directory << " " <<  eval << endl;
	if (eval == false){
		write_directory = write_directory + "/";
	}

}


int Create(string write_directory){
	/// mainCreate creates the backstop, some in as well as the individual markers.

	std::ofstream out;
	string filename = write_directory + "specification_file.txt";
	out.open(filename.c_str());

	int arc_code =  cv::aruco::DICT_6X6_1000;
	Ptr<aruco::Dictionary> dictionary =
			aruco::getPredefinedDictionary(aruco::PREDEFINED_DICTIONARY_NAME(arc_code));

	Size imageSize;
	int squaresX = 12;
	int squaresY = int(float(squaresX)*1.25);
	int squareLength = 200;
	int markerLength = 100;
	int margins = squareLength - markerLength;

	out << "squaresX " << squaresX << endl;
	out << "squaresY " << squaresY << endl;
	out << "squareLength " << squareLength << endl;
	out << "markerLength " << markerLength << endl;
	out << "margins " << margins << endl;
	out << "arc_code " << arc_code << endl;
	out.close();
	imageSize.width = squaresX * squareLength;
	imageSize.height = squaresY * squareLength;


	Mat markerImg;
	Mat boardImage = Mat::zeros(imageSize.height, imageSize.width, CV_8UC1);
	boardImage.setTo(255);
	int x0, y0;

	for (int x = 0, count = 0; x < squaresX; x++){
		for (int y = 0; y < squaresY; y++, count++){
			cout << count << endl;
			aruco::drawMarker(dictionary, count, markerLength, markerImg, 1);

			/// where to place?
			x0 = x*squareLength + margins/2;
			y0 = y*squareLength + margins/2;

			Rect R = Rect(x0, y0, markerLength, markerLength);
			cout << "x0, y0 " << x0 << ", " << y0 << endl;

			markerImg.copyTo(boardImage(R));
		}
	}


	filename = write_directory + "aruco_image" + ToString<int>(squaresX) + "by" + ToString<int>(squaresY) + ".png";
	cv::imwrite(filename, boardImage);


	return 0;
}


int main(int argc, char **argv){
	// create this in  a directory.

	if (argc > 1){
		string write_directory = argv[1];

		EnsureDirHasTrailingBackslash(write_directory);

		struct stat info;
		if( stat( write_directory.c_str(), &info ) != 0 ){
			cout << "Path to write directory is wrong and/or cannot access " << write_directory << endl;
			exit(1);
		}


		Create(write_directory);

	}	else {
		cout << "Please enter a write directory.  Call this program with programname write-directory." << endl;
		exit(1);
	}
}






