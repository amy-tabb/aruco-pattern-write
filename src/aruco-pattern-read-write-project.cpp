//============================================================================
// Name        : TwoD-aruco-create-read-grape-project.cpp
// Author      : Amy Tabb
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
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
#include <iostream>
#include <iomanip>
#include <string>
#include <utility>
#include <vector>
#include <chrono>

//#include "VertexClass.hpp"
//#include "EdgeClass.hpp"
#include "Includes.hpp"
#include "DirectoryFunctions.hpp"



using std::list;
using std::cout;
using std::endl;
using std::string;



template<class T>
std::string FormatWithCommas(T value)
{
	int_type_t uvalue = value;
	bool negative = false;
	if (value < 0){
		negative = true;
		uvalue = -value;
	}


	string s;
	int cnt = 0;
	do
	{
		s.insert(0, 1, char('0' + uvalue % 10));
		uvalue /= 10;
		if (++cnt == 3 && uvalue)
		{
			s.insert(0, 1, ',');
			cnt = 0;
		}
	} while (uvalue);

	if (negative){
		s = "-" + s;
	}
	return s;
}


void WritePlyFile(string outfile, vector<double>& points, vector<int>& faces, vector<int>& color){

	// each vertex needs a color ....

	cout << "Writing to " << outfile << endl;
	std::ofstream out;
	out.open(outfile.c_str());

	out << "ply" << endl;
	out << "format ascii 1.0" << endl;
	out << "element vertex " << points.size()/3 << endl;
	out << "property float x" << endl;
	out << "property float y" << endl;
	out << "property float z" << endl;
	out << "property uchar red" << endl;
	out << "property uchar green" << endl;
	out << "property uchar blue" << endl;
	out << "property uchar alpha" << endl;
	out << "element face " << faces.size()/2 << endl;
	out << "property list uchar int vertex_indices"<< endl;
	out << "end_header" << endl;


	for (uint i = 0; i < (points.size()/3); i++){
		out << points[3*i] << " " <<  points[3*i + 1]  << " " <<  points[3*i + 2]  << " ";
		out << color[0] << " " <<  color[1] << " " << color[2] << " 255" << endl;
	}

	for (int_type_t i = 0; i < (faces.size()/4); i++){
		out << "3 " << faces[4*i]  << " " <<  faces[4*i + 1]   << " " <<  faces[4*i + 2]  << endl; //<< " " << faces[4*i + 3] << endl;;
		out << "3 " << faces[4*i]  << " " <<  faces[4*i + 2]   << " " <<  faces[4*i + 3]  << endl;
	}

	out << endl;

	out.close();
}


void AxisAngleToDCM(vector<vector<double> >& R, vector<double>& axis, float theta){
	// assume axis normalized


	double xs, ys, zs, xC, yC, zC, xyC, yzC, zxC;
	double x, y, z;
	double c, s, C;
	x = axis[0];  y = axis[1];  z = axis[2];
	c = cos(theta);  s = sin(theta); C = 1-c;
	xs = x*s;  ys = y*s; zs = z*s;
	xC = x*C; yC = y*C; zC = z*C;
	xyC = x*yC; yzC = y*zC; zxC = z*xC;
	R[0][0] = x*xC+c;
	R[0][1] = xyC-zs;
	R[0][2] = zxC+ys;

	R[1][0] = xyC+zs;
	R[1][1] = y*yC+c;
	R[1][2] = yzC-xs;

	R[2][0]=  zxC-ys;
	R[2][1] = yzC+xs;
	R[2][2] = z*zC+c;
}


void AxisAngleToDCM(Matrix3d& R, Vector3d& axis, float theta){
	// assume axis normalized


	double xs, ys, zs, xC, yC, zC, xyC, yzC, zxC;
	double x, y, z;
	double c, s, C;
	x = axis(0);  y = axis(1);  z = axis(2);
	c = cos(theta);  s = sin(theta); C = 1-c;
	xs = x*s;  ys = y*s; zs = z*s;
	xC = x*C; yC = y*C; zC = z*C;
	xyC = x*yC; yzC = y*zC; zxC = z*xC;
	R(0, 0) = x*xC+c;
	R(0,1) = xyC-zs;
	R(0,2) = zxC+ys;

	R(1, 0) = xyC+zs;
	R(1, 1) = y*yC+c;
	R(1, 2) = yzC-xs;

	R(2, 0)=  zxC-ys;
	R(2, 1) = yzC+xs;
	R(2, 2) = z*zC+c;
}


//}
using namespace cv;

Vector3d CrossProduct(Point2f& pa, Point2f& pb){
	Vector3d c;

	vector<double> a(3);
	vector<double> b(3);

	a[0] = pa.x;
	a[1] = pa.y;
	a[2] = 1;


	b[0] = pb.x;
	b[1] = pb.y;
	b[2] = 1;

	c(0) = a[1]*b[2] - b[1]*a[2];
	c(1) = -a[0]*b[2] + b[0]*a[2];
	c(2) = a[0]*b[1] - b[0]*a[1];

	return c;

}


Vector3d CrossProduct(Vector3d& a, Vector3d& b){

	Vector3d c;

	c(0) = a(1)*b(2) - b(1)*a(2);
	c(1) = -a(0)*b(2) + b(0)*a(2);
	c(2) = a(0)*b(1) - b(0)*a(1);

	return c;

}

void FitRectangleToSamples(vector<vector<Point2d> >& samples_from_mid, VectorXd& rect_representation){

	/// define Xd before this function.
	//normx normy disp_line0 disp_line1 disp_line2 disp_line3
	/// 6 variables.


	int total_number_samples = 0;

	for (int m = 0; m < 4; m++){
		total_number_samples+= samples_from_mid[m].size();
	}


	//cout << "Total number of samples " << total_number_samples << endl;
	MatrixXd M(total_number_samples, 6);

	int matrix_counter = 0;

	for (int m = 0; m < 4; m++){

		for (int n = 0; n < int(samples_from_mid[m].size()); n++){
			switch (m)
			{
			case 0:{
				M(matrix_counter, 0) = samples_from_mid[m][n].x;
				M(matrix_counter, 1) = samples_from_mid[m][n].y;
				M(matrix_counter, 2) = 1;
				M(matrix_counter, 3) = 0;
				M(matrix_counter, 4) = 0;
				M(matrix_counter, 5) = 0;
			} break;
			case 1:{
				M(matrix_counter, 0) = -samples_from_mid[m][n].y;
				M(matrix_counter, 1) = samples_from_mid[m][n].x;
				M(matrix_counter, 2) = 0;
				M(matrix_counter, 3) = 1;
				M(matrix_counter, 4) = 0;
				M(matrix_counter, 5) = 0;
			} break;
			case 2:{
				M(matrix_counter, 0) = samples_from_mid[m][n].x;
				M(matrix_counter, 1) = samples_from_mid[m][n].y;
				M(matrix_counter, 2) = 0;
				M(matrix_counter, 3) = 0;
				M(matrix_counter, 4) = 1;
				M(matrix_counter, 5) = 0;
			} break;
			case 3:{
				M(matrix_counter, 0) = -samples_from_mid[m][n].y;
				M(matrix_counter, 1) = samples_from_mid[m][n].x;
				M(matrix_counter, 2) = 0;
				M(matrix_counter, 3) = 0;
				M(matrix_counter, 4) = 0;
				M(matrix_counter, 5) = 1;
			} break;
			}

			matrix_counter++;
		}

	}

	cout << "M" << endl << M << endl;

	JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
	cout << "Its singular values are:" << endl << svd.singularValues() << endl;
	//	cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
	cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;

	rect_representation = svd.matrixV().col(5);


}


void FitRectangleToSamples2(vector<vector<Point2d> >& samples_from_mid, VectorXd& rect_representation){

	/// define Xd before this function.
	//normx normy disp_line0 disp_line1 disp_line2 disp_line3
	/// 6 variables.


	int total_number_samples = 0;

	for (int m = 0; m < 4; m++){
		total_number_samples+= samples_from_mid[m].size();
	}


	cout << "Total number of samples " << total_number_samples << endl;
	MatrixXd M0(samples_from_mid[0].size() + samples_from_mid[2].size(), 4);

	int matrix_counter = 0;

	for (int m = 0; m < 4; m = m + 2){

		for (int n = 0; n < int(samples_from_mid[m].size()); n++){
			switch (m)
			{
			case 0:{
				M0(matrix_counter, 0) = samples_from_mid[m][n].x;
				M0(matrix_counter, 1) = samples_from_mid[m][n].y;
				M0(matrix_counter, 2) = 1;
				M0(matrix_counter, 3) = 0;
				//M(matrix_counter, 4) = 0;
				//M(matrix_counter, 5) = 0;
			} break;

			case 2:{
				M0(matrix_counter, 0) = samples_from_mid[m][n].x;
				M0(matrix_counter, 1) = samples_from_mid[m][n].y;
				M0(matrix_counter, 2) = 0;
				M0(matrix_counter, 3) = 1;
				//M(matrix_counter, 4) = 0;
				//M(matrix_counter, 5) = 0;
			} break;
			}

			matrix_counter++;
		}

	}


	//cout << "M" << endl << M0 << endl;

	JacobiSVD<MatrixXd> svd(M0, ComputeThinU | ComputeThinV);
	//cout << "Its singular values are:" << endl << svd.singularValues() << endl;
	//	cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
	//cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;

	VectorXd rect_representation0 = svd.matrixV().col(3);




	//cout << "Total number of samples " << total_number_samples << endl;
	MatrixXd M1(samples_from_mid[1].size() + samples_from_mid[3].size(), 4);

	matrix_counter = 0;

	for (int m = 1; m < 4; m = m + 2){

		for (int n = 0; n < int(samples_from_mid[m].size()); n++){
			switch (m)
			{
			case 1:{
				M1(matrix_counter, 0) = samples_from_mid[m][n].x;
				M1(matrix_counter, 1) = samples_from_mid[m][n].y;
				M1(matrix_counter, 2) = 1;
				M1(matrix_counter, 3) = 0;
				//M(matrix_counter, 4) = 0;
				//M(matrix_counter, 5) = 0;
			} break;

			case 3:{
				M1(matrix_counter, 0) = samples_from_mid[m][n].x;
				M1(matrix_counter, 1) = samples_from_mid[m][n].y;
				M1(matrix_counter, 2) = 0;
				M1(matrix_counter, 3) = 1;
				//M(matrix_counter, 4) = 0;
				//M(matrix_counter, 5) = 0;
			} break;
			}

			matrix_counter++;
		}

	}


	//cout << "M" << endl << M0 << endl;

	JacobiSVD<MatrixXd> svd1(M1, ComputeThinU | ComputeThinV);
	//cout << "Its singular values are:" << endl << svd1.singularValues() << endl;
	//	cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
	//cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd1.matrixV() << endl;

	VectorXd rect_representation1 = svd1.matrixV().col(3);

	rect_representation = VectorXd(8, 1);
	for (int i = 0; i < 4; i++){
		rect_representation(i) = rect_representation0(i);
	}

	for (int i = 0; i < 4; i++){
		rect_representation(i + 4) = rect_representation1(i);
	}



	/// then use this to do the other ones.

	//	for (int m = 0; m < 4; m++){
	//
	//		for (int n = 0; n < int(samples_from_mid[m].size()); n++){
	//			switch (m)
	//			{
	//			case 0:{
	//				M(matrix_counter, 0) = samples_from_mid[m][n].x;
	//				M(matrix_counter, 1) = samples_from_mid[m][n].y;
	//				M(matrix_counter, 2) = 1;
	//				M(matrix_counter, 3) = 0;
	//				M(matrix_counter, 4) = 0;
	//				M(matrix_counter, 5) = 0;
	//			} break;
	//			case 1:{
	//				M(matrix_counter, 0) = -samples_from_mid[m][n].y;
	//				M(matrix_counter, 1) = samples_from_mid[m][n].x;
	//				M(matrix_counter, 2) = 0;
	//				M(matrix_counter, 3) = 1;
	//				M(matrix_counter, 4) = 0;
	//				M(matrix_counter, 5) = 0;
	//			} break;
	//			case 2:{
	//				M(matrix_counter, 0) = samples_from_mid[m][n].x;
	//				M(matrix_counter, 1) = samples_from_mid[m][n].y;
	//				M(matrix_counter, 2) = 0;
	//				M(matrix_counter, 3) = 0;
	//				M(matrix_counter, 4) = 1;
	//				M(matrix_counter, 5) = 0;
	//			} break;
	//			case 3:{
	//				M(matrix_counter, 0) = -samples_from_mid[m][n].y;
	//				M(matrix_counter, 1) = samples_from_mid[m][n].x;
	//				M(matrix_counter, 2) = 0;
	//				M(matrix_counter, 3) = 0;
	//				M(matrix_counter, 4) = 0;
	//				M(matrix_counter, 5) = 1;
	//			} break;
	//			}
	//
	//			matrix_counter++;
	//		}
	//
	//	}
	//
	//	cout << "M" << endl << M << endl;
	//
	//	JacobiSVD<MatrixXd> svd(M, ComputeThinU | ComputeThinV);
	//	cout << "Its singular values are:" << endl << svd.singularValues() << endl;
	//	//	cout << "Its left singular vectors are the columns of the thin U matrix:" << endl << svd.matrixU() << endl;
	//	cout << "Its right singular vectors are the columns of the thin V matrix:" << endl << svd.matrixV() << endl;
	//
	//	rect_representation = svd.matrixV().col(5);


}


//Function for average
double avg ( vector<double>& v )
{
	double return_value = 0.0;
	int n = v.size();

	for ( int i=0; i < n; i++)
	{
		return_value += v[i];
	}

	return ( return_value / double(n));
}
//****************End of average funtion****************


//Function for variance
double variance ( vector<double>& v , double mean )
{
	double sum = 0.0;
	double temp =0.0;
	double var =0.0;
	int n = v.size();

	for ( int j =0; j < n; j++)
	{
		temp = pow(v[j] - mean,2);
		sum += temp;
	}

	return var = sum/(v.size() -1);
}

bool reverse_sort( pair<double, int> p0, pair<double, int> p1){

	if (p0.first > p1.first){
		return true;
	}	else {
		return false;
	}
}

void local_undistort_points(vector<Point2f>& distorted_points, vector<Point2f>& undistorted_points, Mat& cameraMatrix, Mat& distCoeffs){

	Point2f dist_norm;
	Point2f undist_norm;
	Point2f undist_px;

	double r_squared;
	double coeff;
	cout << "camera matrix " << endl << cameraMatrix << endl << "dist Coeffs " << endl << distCoeffs << endl;
	for (int m = 0; m < int(distorted_points.size() ); m++){

		dist_norm.x = (distorted_points[m].x - cameraMatrix.at<double>(0,2))/cameraMatrix.at<double>(0,0);
		dist_norm.y = (distorted_points[m].y - cameraMatrix.at<double>(1,2))/cameraMatrix.at<double>(1,1);

		r_squared = pow(dist_norm.x, 2) + pow(dist_norm.y, 2);
		coeff = 1 + distCoeffs.at<double>(0, 0)*r_squared + distCoeffs.at<double>(0, 1)*r_squared*r_squared;
		cout << "coeff " << coeff << endl;

		//coeff = 1; /// test
		undist_norm.x = dist_norm.x*coeff;
		undist_norm.y = dist_norm.y*coeff;

		undist_px.x = cameraMatrix.at<double>(0,0)*undist_norm.x + cameraMatrix.at<double>(0,2);
		undist_px.y = cameraMatrix.at<double>(1,1)*undist_norm.y + cameraMatrix.at<double>(1,2);



		undistorted_points.push_back(undist_px);
		cout << "dist original " << distorted_points[m] << endl;
		cout << "dist norm " << dist_norm << endl;
		cout << "undist norm " << undist_norm << endl;
		cout << "undist px " << undist_px << endl;

	}
}

bool sort_by_x( pair<Point2f, Point2f> p0, pair<Point2f, Point2f> p1){

	if (p0.first.x < p1.first.x){
		return true;
	}	else {
		return false;
	}
}

//https://docs.opencv.org/3.4.1/d9/d6d/tutorial_table_of_content_aruco.html
int mainCreate(int argc, char **argv){
	/// mainCreate creates the backstop, some in as well as the individual markers.

	/// Print out a file the registers the markers on the board
	/// everything is written to the current directory ... probably not the best design.
	std::ofstream out;
	string filename = "specification_file.txt";
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


	int id_start_number  = squaresX*squaresY + 100;
	out << "squaresX " << squaresX << endl;
	out << "squaresY " << squaresY << endl;
	out << "squareLength " << squareLength << endl;
	out << "markerLength " << markerLength << endl;
	out << "margins " << margins << endl;
	out << "id_start_number " << id_start_number << endl;
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


	filename = "backstop" + ToString<int>(squaresX) + "by" + ToString<int>(squaresY) + ".png";
	cv::imwrite(filename, boardImage);

	Mat markerImageLarge = Mat::zeros(squareLength, squareLength, CV_8UC1);
	markerImageLarge.setTo(255);
	for (int id = id_start_number, count  = 0; count < 8; count++, id++){

		aruco::drawMarker(dictionary, id, markerLength, markerImg, 1);
		Rect R = Rect(margins/2, margins/2, markerLength, markerLength);
		markerImg.copyTo(markerImageLarge(R));
		filename = "markerID" + ToString<int>(id) + ".png";
		cv::imwrite(filename, markerImageLarge);
	}

	return 0;
}

//
//
//
//
//	Ptr<aruco::CharucoBoard> board = aruco::CharucoBoard::create(squaresX, squaresY, (float)squareLength,
//			(float)markerLength, dictionary);
//
//	// show created board
//	Mat boardImage;
//	board->draw(imageSize, boardImage, margins, 2);

//	if(showImage) {
//		imshow("board", boardImage);
//		waitKey(0);
//	}

//imwrite(out, boardImage);


/**
 */
static bool readDetectorParameters(string filename, Ptr<aruco::DetectorParameters> &params) {
    FileStorage fs(filename, FileStorage::READ);
    if(!fs.isOpened())
        return false;
    fs["adaptiveThreshWinSizeMin"] >> params->adaptiveThreshWinSizeMin;
    fs["adaptiveThreshWinSizeMax"] >> params->adaptiveThreshWinSizeMax;
    fs["adaptiveThreshWinSizeStep"] >> params->adaptiveThreshWinSizeStep;
    fs["adaptiveThreshConstant"] >> params->adaptiveThreshConstant;
    fs["minMarkerPerimeterRate"] >> params->minMarkerPerimeterRate;
    fs["maxMarkerPerimeterRate"] >> params->maxMarkerPerimeterRate;
    fs["polygonalApproxAccuracyRate"] >> params->polygonalApproxAccuracyRate;
    fs["minCornerDistanceRate"] >> params->minCornerDistanceRate;
    fs["minDistanceToBorder"] >> params->minDistanceToBorder;
    fs["minMarkerDistanceRate"] >> params->minMarkerDistanceRate;
    //fs["doCornerRefinement"] >> params->doCornerRefinement;
    fs["cornerRefinementWinSize"] >> params->cornerRefinementWinSize;
    fs["cornerRefinementMaxIterations"] >> params->cornerRefinementMaxIterations;
    fs["cornerRefinementMinAccuracy"] >> params->cornerRefinementMinAccuracy;
    fs["markerBorderBits"] >> params->markerBorderBits;
    fs["perspectiveRemovePixelPerCell"] >> params->perspectiveRemovePixelPerCell;
    fs["perspectiveRemoveIgnoredMarginPerCell"] >> params->perspectiveRemoveIgnoredMarginPerCell;
    fs["maxErroneousBitsInBorderRate"] >> params->maxErroneousBitsInBorderRate;
    fs["minOtsuStdDev"] >> params->minOtsuStdDev;
    fs["errorCorrectionRate"] >> params->errorCorrectionRate;
    return true;
}



int mainRead(int argc, char **argv){
	Eigen::initParallel();




	string read_dir = "";
	string write_dir = "";
	string dir_decorator = "";

	if (argc >= 3){
		read_dir = argv[1];
		write_dir = argv[2];
	}	else {
		cout << "Not enough arguments.  prog_name input_filename write_filename number_dilations." << endl;
		exit(1);
	}

	vector<string> image_names;

	ReadDirectory(read_dir, image_names);


	Ptr<aruco::DetectorParameters> detectorParams = aruco::DetectorParameters::create();

	// detectorParams->ccornerRefinementMethod = aruco::CORNER_REFINE_SUBPIX; // do corner refinement in markers
	bool readOk = readDetectorParameters(string("../src/detector_params.yml"), detectorParams);
	if(!readOk) {
		cout << "Invalid detector parameters file" << endl;
		return 0;
	}

	Ptr<aruco::Dictionary> dictionary =
			aruco::getPredefinedDictionary(aruco::PREDEFINED_DICTIONARY_NAME(cv::aruco::DICT_6X6_1000));



	for (int i = 0, in = image_names.size(); i < in; i++){
		Mat image, imageCopy;
		string filename = read_dir + "/" + image_names[i];
		cout << "filename " << filename << endl;
		image = imread(filename.c_str());

		vector< int > ids;
		vector< vector< Point2f > > corners, rejected;
		vector< Vec3d > rvecs, tvecs;

		// detect markers and estimate pose
		aruco::detectMarkers(image, dictionary, corners, ids, detectorParams, rejected);

		// draw results
		image.copyTo(imageCopy);
		if(ids.size() > 0) {
			aruco::drawDetectedMarkers(imageCopy, corners, ids);
			filename = write_dir + "/" + image_names[i];
			cv::imwrite(filename, imageCopy);
		}
	}

	return 0;

}

int main(int argc, char **argv){

	if (argc == 3){
		mainRead(argc, argv);
	}	else {
		mainCreate(argc, argv);
	}
}






