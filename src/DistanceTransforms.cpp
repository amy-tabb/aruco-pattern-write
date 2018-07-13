//#include "SkelGraph.hpp"
#include "DistanceTransforms.hpp"
//#include "SubSteps.hpp"

double f(int_type_t x, int_type_t i, int_type_t y, vector<vector<double> >& g){
	return pow((double(x) - double(i)), 2) + pow(g[i][y], 2);
}

double fwith_square(int_type_t x, int_type_t i, int_type_t y, vector<vector<double> >& g){

	return pow((double(x) - double(i)), 2) + pow(g[i][y], 2);
}

// for the non-grid version, scan
double f(int_type_t x, int_type_t i, double* g, int_type_t* line){

	return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

//// for the z-dominant edge distance version, scan
//double f(int_type_t x, int_type_t i, vector<double*>& g, int_type_t* line, ReconstructionStructure& RS){
//
//	int_type_t xi, yi, zi;
//	RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
//
//	return pow((double(x) - double(i)), 2) + pow(g[zi][xi*RS.number_voxels_per_dim[1] + yi], 2);
//	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
//}

// for the z-dominant edge distance version, scan
double f(int_type_t x, int_type_t i, vector<double*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){

	return pow((double(x) - double(i)), 2) + pow(g[z][x_pos_in_line[i]*ysize + y], 2);
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}


// for the z-dominant edge distance version, scan
float f(int_type_t x, int_type_t i, vector<float*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){

	return pow((float(x) - float(i)), 2) + pow(g[z][x_pos_in_line[i]*ysize + y], 2);
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

// TODO double check this one..
// for the z-dominant edge distance version, scan
//uint_dist_type f(int_type_t x, int_type_t i, vector< uint_dist_type*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){
//
//	uint_dist_type first_term = 0;
//	x > i ? first_term = (x - i)*(x - i) : first_term = (i - x)*(i - x);
//	//return pow((float(x) - float(i)), 2) + pow(g[z][x_pos_in_line[i]*ysize + y], 2);
//
//	return first_term + pow(g[z][x_pos_in_line[i]*ysize + y], 2);
//	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
//}

uint_dist_type f(int_type_t x, int_type_t i, vector<uint_dist_type*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){
	int64_t val = pow((int64_t(x) - int64_t(i)), 2) + pow(int64_t(g[z][x_pos_in_line[i]*ysize + y]), 2);
	// downgrade
	return uint_dist_type(val);
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}



// for the z-dominant edge distance version, scan
double f(int_type_t x, int_type_t i, double* g){

	return pow((double(x) - double(i)), 2) + pow(g[i], 2);
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

uint_dist_type f(int_type_t x, int_type_t i, vector<uint_dist_type*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize, uint_dist_type big_number){
	int64_t val = pow((int64_t(x) - int64_t(i)), 2) + pow(int64_t(g[z][x_pos_in_line[i]*ysize + y]), 2);
	// downgrade

	if (val > big_number){
		return big_number;
	}	else {
		return uint_dist_type(val);
	}
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

//
//// for the z-dominant edge distance version, scan
//double f(int_type_t x, int_type_t i, int_type_t y, vector<double*>& g, int_type_t z, ReconstructionStructure& RS){
//
//	return pow((double(x) - double(i)), 2) + pow(g[z][i*RS.number_voxels_per_dim[1] + y], 2);
//	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
//}
//
//// for the z-dominant edge distance version, scan
//double f_with_square(int_type_t x, int_type_t i, int_type_t x0, int_type_t y0, vector<double*>& g, ReconstructionStructure& RS){
//
//	return pow((double(x) - double(i)), 2) + (g[i][x0*RS.number_voxels_per_dim[1] + y0]); /// g already squared
//	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
//}

// for the z-dominant edge distance version, scan
double f_with_square(int_type_t x, int_type_t i, int_type_t x0, int_type_t y0, vector<double*>& g, int_type_t* z_pos_in_line, int_type_t ysize){

	return pow((double(x) - double(i)), 2) + (g[z_pos_in_line[i]][x0*ysize + y0]); /// g already squared
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

// for the z-dominant edge distance version, scan
double f_with_square(int_type_t x, int_type_t i, double* g){

	return pow((double(x) - double(i)), 2) + (g[i]); /// g already squared
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}

float f_with_square(int_type_t x, int_type_t i, int_type_t x0, int_type_t y0, vector<float*>& g, int_type_t* z_pos_in_line, int_type_t ysize){

	return pow((float(x) - float(i)), 2) + (g[z_pos_in_line[i]][x0*ysize + y0]); /// g already squared
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}


uint_dist_type f_with_square(int_type_t x, int_type_t i, int_type_t x0, int_type_t y0, vector<uint_dist_type*>& g, int_type_t* z_pos_in_line, int_type_t ysize){
	uint_dist_type result = 0;
	x > i ? result =  (x - i)*(x - i) + (g[z_pos_in_line[i]][x0*ysize + y0]) :  result = (i - x)*(i - x) + (g[z_pos_in_line[i]][x0*ysize + y0]);

	return result;
	//return pow((float(x) - float(i)), 2) + (g[z_pos_in_line[i]][x0*ysize + y0]); /// g already squared
	//return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
}
//
//// for the z-dominant edge distance version, scan
//double Sep_with_square(int_type_t i, int_type_t u, int_type_t x0, int_type_t y0, vector<double*>& g,
//		ReconstructionStructure& RS, double big_number){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	double s = int(double(u*u) - double(i*i)
//			+ g[u][x0*RS.number_voxels_per_dim[1] + y0] - g[i][x0*RS.number_voxels_per_dim[1] + y0])/int(2*(double(u - i)));
//
//	//	if (s < 0){
//	//		cout << "Error!  Sep less than zero " << s << endl;
//	//		exit(1);
//	//	}
//
//	if (s < 0){
//
//		if (fabs(s) < 0.00001){
//			s = 0;
//		}	else {
//			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
//			if (abs(g[u][x0*RS.number_voxels_per_dim[1] + y0] - big_number) < 0.01){
//				//cout << "Triggering big number case .... " << endl;
//				return u;
//			}	else {
//				cout << "u, i " << u << ", " << i << endl;
//				cout << "gu, gi " <<  g[u][x0*RS.number_voxels_per_dim[1] + y0] << ", " << g[i][x0*RS.number_voxels_per_dim[1] + y0] << endl;
//				double numer = int(double(u*u) - double(i*i)
//						+ g[u][x0*RS.number_voxels_per_dim[1] + y0] - g[i][x0*RS.number_voxels_per_dim[1] + y0]);
//				double denom = int(2*(double(u - i)));
//				cout << "number, denom " << numer << ", " << denom << endl;
//
//
//				cout << "Error!  Sep less than zero " << s << endl;
//				exit(1);
//			}
//		}
//	}
//	return s;
//}
//





// for the z-dominant edge distance version, scan
double Sep_with_square(int_type_t i, int_type_t u, int_type_t x0, int_type_t y0, vector<double*>& g,
		int_type_t* z_pos_in_line, int_type_t ysize, double big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double g_value_u = g[z_pos_in_line[u]][x0*ysize + y0];
	double g_value_i = g[z_pos_in_line[i]][x0*ysize + y0];



	//	double s = int(double(u*u) - double(i*i)
	//			+ g[u][x0*RS.number_voxels_per_dim[1] + y0] - g[i][x0*RS.number_voxels_per_dim[1] + y0])/int(2*(double(u - i)));

	double s = int(double(u*u) - double(i*i)
			+ g_value_u - g_value_i)/int(2*(double(u - i)));

	//	if (s < 0){
	//		cout << "Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}

	if (s < 0){

		if (fabs(s) < 0.00001){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if (fabs(g_value_u - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				return u;
			}	else {
				cout << "u, i " << u << ", " << i << endl;
				cout << "gu, gi " <<  g_value_u << ", " << g_value_i << endl;
				double numer = int(double(u*u) - double(i*i)
						+ g_value_u - g_value_i);
				double denom = int(2*(double(u - i)));
				cout << "number, denom " << numer << ", " << denom << endl;


				cout << "Error!  Sep_sq less than zero " << s << endl;
				exit(1);
			}
		}
	}
	return s;
}

// for the z-dominant edge distance version, scan
double Sep_with_square(int_type_t i, int_type_t u, double* g, double big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double g_value_u = g[u];
	double g_value_i = g[i];

	double s = int(double(u*u) - double(i*i)
			+ g_value_u - g_value_i)/int(2*(double(u - i)));

	if (s < 0){

		if (fabs(s) < 0.00001){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if (fabs(g_value_u - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				return u;
			}	else {
				cout << "u, i " << u << ", " << i << endl;
				cout << "gu, gi " <<  g_value_u << ", " << g_value_i << endl;
				double numer = int(double(u*u) - double(i*i)
						+ g_value_u - g_value_i);
				double denom = int(2*(double(u - i)));
				cout << "number, denom " << numer << ", " << denom << endl;


				cout << "Error!  Sep_sq less than zero " << s << endl;
				exit(1);
			}
		}
	}
	return s;
}

// for the z-dominant edge distance version, scan
float Sep_with_square(int_type_t i, int_type_t u, int_type_t x0, int_type_t y0, vector<float*>& g,
		int_type_t* z_pos_in_line, int_type_t ysize, float big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	float g_value_u = g[z_pos_in_line[u]][x0*ysize + y0];
	float g_value_i = g[z_pos_in_line[i]][x0*ysize + y0];



	//	double s = int(double(u*u) - double(i*i)
	//			+ g[u][x0*RS.number_voxels_per_dim[1] + y0] - g[i][x0*RS.number_voxels_per_dim[1] + y0])/int(2*(double(u - i)));

	float s = int(float(u*u) - float(i*i)
			+ g_value_u - g_value_i)/int(2*(float(u - i)));

	//	if (s < 0){
	//		cout << "Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}

	if (s < 0){

		if (fabs(s) < 0.01){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if (fabs(g_value_u - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				return u;
			}	else {
				cout << "u, i " << u << ", " << i << endl;
				cout << "gu, gi " <<  g_value_u << ", " << g_value_i << endl;
				float numer = int(float(u*u) - float(i*i)
						+ g_value_u - g_value_i);
				float denom = int(2*(float(u - i)));
				cout << "number, denom " << numer << ", " << denom << endl;


				cout << "Error!  Sep_sq less than zero " << s << endl;
				exit(1);
			}
		}
	}
	return s;
}

// for the z-dominant edge distance version, scan
uint_dist_type Sep_with_square(int_type_t i, int_type_t u, int_type_t x0, int_type_t y0, vector<uint_dist_type*>& g,
		int_type_t* z_pos_in_line, int_type_t ysize, uint_dist_type big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	int64_t g_value_u = g[z_pos_in_line[u]][x0*ysize + y0];
	int64_t g_value_i = g[z_pos_in_line[i]][x0*ysize + y0];



	//	double s = int(double(u*u) - double(i*i)
	//			+ g[u][x0*RS.number_voxels_per_dim[1] + y0] - g[i][x0*RS.number_voxels_per_dim[1] + y0])/int(2*(double(u - i)));

	int64_t s = (int64_t(u*u) - int64_t(i*i)
			+ g_value_u - g_value_i)/(2*((int64_t(u) - int64_t(i))));

	//	if (s < 0){
	//		cout << "Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}

	if (s < 0){
		if (g_value_u == big_number){
			return u;
		}

		if (g_value_i >= big_number){
			return i;
		}

		if (g_value_u >= big_number){
			return u;
		}

		cout << "u, i " << u << ", " << i << endl;
		cout << "gu, gi " <<  g_value_u << ", " << g_value_i << endl;
		int64_t numer = (int64_t(u*u) - int64_t(i*i)
				+ g_value_u - g_value_i);
		int64_t denom = (2*((int64_t(u) - int64_t(i))));
		cout << "number, denom " << numer << ", " << denom << endl;
		cout << "big number " << big_number << endl;

		cout << "Error!  Sep_sq less than zero " << s << endl;
		exit(1);
	}
	return uint_dist_type(s);
}




//double f(int_type_t x, int_type_t i, int_type_t y, vector<vector<double> >& g){
//	return pow((double(x) - double(i)), 2) + pow(g[i][y], 2);
//}

// for the non-grid version, scan
double f_with_square(int_type_t x, int_type_t i, double* g, int_type_t* line){

	return pow((double(x) - double(i)), 2) + g[line[i]];
}

// for the non-grid version, scan
double f(int_type_t x, int_type_t i, double* g, int_type_t* line, int_type_t big_number){

	if (line[i] == big_number){
		// this value is going to be zero.
		return pow((double(x) - double(i)), 2);
	}	else {
		return pow((double(x) - double(i)), 2) + pow(g[line[i]], 2);
	}
}



// NOte: here g is the squared distance
double f(int_type_t z, int_type_t i, int_type_t x, int_type_t y, vector<vector<vector<double> > >& g){
	return pow((double(z) - double(i)), 2) + g[x][y][i];
}

double Sep(int_type_t i, int_type_t u, int_type_t y, vector<vector<double> >& g){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));


	//	if (s < 0){
	//		cout << "Error!  Sep less than zero " << s << endl;
	//
	//	}

	if (s < 0){

		if (fabs(s) < 0.00001){
			s = 0;
		}	else {
			cout << "i, u, y, gu, gy " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			cout << "83 Error!  Sep less than zero " << s << endl;
			exit(1);
		}
	}
	return s;
}

double Sep(int_type_t i, int_type_t u, int_type_t y, vector<vector<double> >& g, double big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));


	//	if (s < 0){
	//		cout << "Error!  Sep less than zero " << s << endl;
	//
	//	}

	if (s < 0){

		if (fabs(s) < 0.00001){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if (abs(g[u][y] - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				return u;
			}	else {
				cout << "Big number " << big_number << ", " << abs(g[u][y] - big_number) << endl;
				cout << "115 Error!  Sep less than zero " << s << endl;
				exit(1);
			}
		}
	}
	return s;
}

// scan 4, non-grid version
double Sep(int_type_t i, int_type_t u, double* g, int_type_t* line){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + pow(g[line[u]], 2) - pow(g[line[i]], 2))/int(2*(double(u - i)));

	if (s < 0){
		cout << "100 Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}
//
//// scan 4, semi=-grid version with z-dominant indexing.
//double Sep(int_type_t i, int_type_t u, vector<double*>& g, int_type_t* line, ReconstructionStructure& RS){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	int_type_t xi, yi, zi;
//	RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
//	double g_value_u = g[zi][xi*RS.number_voxels_per_dim[1] + yi];
//
//	RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
//	double g_value_i = g[zi][xi*RS.number_voxels_per_dim[1] + yi];
//
//	double s = int(double(u*u) - double(i*i) + pow(g_value_u, 2) - pow(g_value_i, 2))/int(2*(double(u - i)));
//
//	if (s < 0){
//		cout << "123 Error!  Sep less than zero " << s << endl;
//		exit(1);
//	}
//	return s;
//}

//pow(g[z][x_pos_in_line[i]*ysize + y], 2);
// scan 4, semi=-grid version with z-dominant indexing.
double Sep(int_type_t i, int_type_t u, vector<double*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize, double big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	//int_type_t xi, yi, zi;
	//RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
	double g_value_u = g[z][x_pos_in_line[u]*ysize + y];

	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
	double g_value_i = g[z][x_pos_in_line[i]*ysize + y];

	double s = int(double(u*u) - double(i*i) + pow(g_value_u, 2) - pow(g_value_i, 2))/int(2*(double(u - i)));

	//	if (s < 0){
	//		cout << "i, u " << i << ", " << u << endl;
	//		cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
	//		cout << "123 Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}
	//	return s;


	if (s < 0){
		if (fabs(s) < 0.00001){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if ((g_value_u - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				// shouldn't this be i?
				return u;
			}	else {

				cout << "i, u " << i << ", " << u << endl;
				cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
				cout << "269 Error!  Sep less than zero " << s << endl;
				exit(1);
			}
		}
	}


	return s;
}


double Sep(int_type_t i, int_type_t u,  double* g, double big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	//int_type_t xi, yi, zi;
	//RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
	double g_value_u = g[u];

	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
	double g_value_i = g[i];

	double s = int(double(u*u) - double(i*i) + pow(g_value_u, 2) - pow(g_value_i, 2))/int(2*(double(u - i)));

	//	if (s < 0){
	//		cout << "i, u " << i << ", " << u << endl;
	//		cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
	//		cout << "123 Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}
	//	return s;


	if (s < 0){
		if (fabs(s) < 0.00001){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if ((g_value_u - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				// shouldn't this be i?
				return u;
			}	else {

				cout << "i, u " << i << ", " << u << endl;
				cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
				cout << "269 Error!  Sep less than zero " << s << endl;
				exit(1);
			}
		}
	}


	return s;
}

float Sep(int_type_t i, int_type_t u, vector<float*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize, float big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	//int_type_t xi, yi, zi;
	//RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
	float g_value_u = g[z][x_pos_in_line[u]*ysize + y];

	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
	float g_value_i = g[z][x_pos_in_line[i]*ysize + y];

	float s = int(float(u*u) - float(i*i) + pow(g_value_u, 2) - pow(g_value_i, 2))/int(2*(float(u - i)));

	//	if (s < 0){
	//		cout << "i, u " << i << ", " << u << endl;
	//		cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
	//		cout << "123 Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}
	//	return s;


	if (s < 0){
		if (fabs(s) < 0.001){
			s = 0;
		}	else {
			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
			if ((g_value_u - big_number) < 0.01){
				//cout << "Triggering big number case .... " << endl;
				// shouldn't this be i?
				return u;
			}	else {

				cout << "i, u " << i << ", " << u << endl;
				cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
				cout << "269 Error!  Sep less than zero " << s << endl;
				exit(1);
			}
		}
	}


	return s;
}

///// TODO go through this one extensively.
//uint_dist_type Sep(int_type_t i, int_type_t u, vector<uint_dist_type*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize,
//		uint_dist_type big_number){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	uint_dist_type g_value_u = g[z][x_pos_in_line[u]*ysize + y];
//
//	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
//	uint_dist_type g_value_i = g[z][x_pos_in_line[i]*ysize + y];
//	// go big -- don't know what might happen here.
//	int64_t positive_terms = u*u + pow(g_value_u, 2);
//	int64_t negative_terms = i*i + pow(g_value_i, 2);
//	// use signed to deal with possible negative.
//	int64_t s = (positive_terms - negative_terms)/(2*((int(u) - int(i))));
//	/// possibility -- the denominator would flip this ... ???
//	if ( s < 0 ){
//		if (g_value_u == big_number){
//			return u;
//		}	else {
//			cout << "i, u " << i << ", " << u << endl;
//			cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
//			cout << "Big number " << big_number << endl;
//			cout << "269 Error!  Sep less than zero " << s << endl;
//			exit(1);
//		}
//	}	else {
//		// downward cast back to unsigned int.  Should be fine.
//		return uint_dist_type(s);
//	}
//
////	if (s < 0){
////		if (fabs(s) < 0.001){
////			s = 0;
////		}	else {
////			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
////			if ((g_value_u - big_number) < 0.01){
////				//cout << "Triggering big number case .... " << endl;
////				// shouldn't this be i?
////				return u;
////			}	else {
////
////				cout << "i, u " << i << ", " << u << endl;
////				cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
////				cout << "269 Error!  Sep less than zero " << s << endl;
////				exit(1);
////			}
////		}
////	}
////
////
////	return s;
//}


//uint_dist_type Sep(int_type_t i, int_type_t u, vector<uint_dist_type*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize, uint_dist_type big_number){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	//int_type_t xi, yi, zi;
//	//RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
//	int64_t g_value_u = g[z][x_pos_in_line[u]*ysize + y];
//
//	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
//	int64_t g_value_i = g[z][x_pos_in_line[i]*ysize + y];
//
//	int64_t s = (int64_t(u*u) - int64_t(i*i) + pow(int64_t(g_value_u), 2) - pow(int64_t(g_value_i), 2))/int64_t(2*((int64_t(u) - int64_t(i))));
//
//	//	if (s < 0){
//	//		cout << "i, u " << i << ", " << u << endl;
//	//		cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
//	//		cout << "123 Error!  Sep less than zero " << s << endl;
//	//		exit(1);
//	//	}
//	//	return s;
//
//
//	if (s < 0){
//		if (fabs(s) < 0.001){
//			s = 0;
//		}	else {
//			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
//			if ((g_value_u - big_number) < 0.01){
//				//cout << "Triggering big number case .... " << endl;
//				// shouldn't this be i?
//				return u;
//			}	else {
//
//				cout << "i, u " << i << ", " << u << endl;
//				cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
//				cout << "693 Error!  Sep less than zero " << s << endl;
//				cout << "Big number " << big_number << endl;
//				exit(1);
//			}
//		}
//	}
//
//	/// downward cast -- should be okay.
//	return uint_dist_type(s);
//}

uint_dist_type Sep(int_type_t i, int_type_t u, vector<uint_dist_type*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize, uint_dist_type big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	//int_type_t xi, yi, zi;
	//RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
	int64_t g_value_u = g[z][x_pos_in_line[u]*ysize + y];

	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
	int64_t g_value_i = g[z][x_pos_in_line[i]*ysize + y];

	int64_t s = (int64_t(u*u) - int64_t(i*i) + pow(int64_t(g_value_u), 2) - pow(int64_t(g_value_i), 2))/int64_t(2*((int64_t(u) - int64_t(i))));

	//	if (s < 0){
	//		cout << "i, u " << i << ", " << u << endl;
	//		cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
	//		cout << "123 Error!  Sep less than zero " << s << endl;
	//		exit(1);
	//	}
	//	return s;
	if (g_value_u >= big_number){
		return u;
	}

	if (s < 0){
		if (fabs(s) < 0.001){
			s = 0;
		}	else {

			if (g_value_u >= big_number){
				return u;
			}

			//			if (g_value_i >= big_number){
			//				return u;
			//			}

			cout << "i, u " << i << ", " << u << endl;
			cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
			cout << "693 Error!  Sep less than zero " << s << endl;
			cout << "Big number " << big_number << endl;
			exit(1);
		}
	}

	/// downward cast -- should be okay.
	return uint_dist_type(s);
}


float Sep(int_type_t i, int_type_t u, vector<float*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){
	/// big number flag is -1
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	//int_type_t xi, yi, zi;
	//RS.ReturnXYZIndicesFromIndex(line[u], xi, yi, zi);
	float g_value_u = g[z][x_pos_in_line[u]*ysize + y];

	//RS.ReturnXYZIndicesFromIndex(line[i], xi, yi, zi);
	float g_value_i = g[z][x_pos_in_line[i]*ysize + y];

	if (g_value_u < 0){
		return u;
	}	else {
		float s = int(float(u*u) - float(i*i) + pow(g_value_u, 2) - pow(g_value_i, 2))/int(2*(float(u - i)));

		//	if (s < 0){
		//		cout << "i, u " << i << ", " << u << endl;
		//		cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
		//		cout << "123 Error!  Sep less than zero " << s << endl;
		//		exit(1);
		//	}
		//	return s;


		if (s < 0){
			if (fabs(s) < 0.00001){
				s = 0;
			}	else {


				cout << "i, u " << i << ", " << u << endl;
				cout << "gi, gu " << g_value_i << ", " << g_value_u << endl;
				cout << "269 Error!  Sep less than zero " << s << endl;
				exit(1);
			}
		}


		return s;
	}


}

//
//double Sep(int_type_t i, int_type_t u, int_type_t y, vector<double*>& g, int_type_t z,  ReconstructionStructure& RS){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	//double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));
//
//	double s = int(double(u*u) - double(i*i) + pow(g[z][u*RS.number_voxels_per_dim[1] + y], 2) - pow(g[z][i*RS.number_voxels_per_dim[1] + y], 2))/int(2*(double(u - i)));
//
//	//	if (s < 0){
//	//		cout << "Error!  Sep less than zero " << s << endl;
//	//		exit(1);
//	//	}
//
//	if (s < 0){
//
//		if (fabs(s) < 0.00001){
//			s = 0;
//		}	else {
//			cout << "149 Error!  Sep less than zero " << s << endl;
//			exit(1);
//		}
//	}
//	return s;
//}


//
//double Sep(int_type_t i, int_type_t u, int_type_t y, vector<double*>& g, int_type_t z,  ReconstructionStructure& RS, double big_number){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	//double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));
//
//	double s = int(double(u*u) - double(i*i) + pow(g[z][u*RS.number_voxels_per_dim[1] + y], 2) - pow(g[z][i*RS.number_voxels_per_dim[1] + y], 2))/int(2*(double(u - i)));
//
//	//	if (s < 0){
//	//		cout << "Error!  Sep less than zero " << s << endl;
//	//		exit(1);
//	//	}
//
//	//	if (s < 0){
//	//
//	//		if (fabs(s) < 0.00001){
//	//			s = 0;
//	//		}	else {
//	//			cout << "149 Error!  Sep less than zero " << s << endl;
//	//			exit(1);
//	//		}
//	//	}
//
//	if (s < 0){
//
//		if (fabs(s) < 0.00001){
//			s = 0;
//		}	else {
//			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
//			if (abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number) < 0.01){
//				//cout << "Triggering big number case .... " << endl;
//				// shouldn't this be i?
//				return u;
//			}	else {
//				//cout << "Big number " << big_number << ", " << abs(g[u][y] - big_number) << endl;
//				cout << "i " << i << " u " << u << endl;
//				cout << "Big number -- u side " << big_number << ", " << abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number) << endl;
//				cout << " i side              " << g[z][i*RS.number_voxels_per_dim[1] + y] << endl;
//				cout << "229 Error!  Sep less than zero " << s << endl;
//				exit(1);
//			}
//		}
//	}
//
//
//	return s;
//}
//
//double SepTemp(int_type_t i, int_type_t u, int_type_t y, vector<double*>& g, int_type_t z,  ReconstructionStructure& RS, double big_number, double big_number_temp){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	//double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));
//
//	double s = int(double(u*u) - double(i*i) + pow(g[z][u*RS.number_voxels_per_dim[1] + y], 2) - pow(g[z][i*RS.number_voxels_per_dim[1] + y], 2))/int(2*(double(u - i)));
//
//	//	if (s < 0){
//	//		cout << "Error!  Sep less than zero " << s << endl;
//	//		exit(1);
//	//	}
//
//	//	if (s < 0){
//	//
//	//		if (fabs(s) < 0.00001){
//	//			s = 0;
//	//		}	else {
//	//			cout << "149 Error!  Sep less than zero " << s << endl;
//	//			exit(1);
//	//		}
//	//	}
//
//	if (s < 0){
//
//		if (fabs(s) < 0.00001){
//			s = 0;
//		}	else {
//			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
//			if (abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number) < 0.01 || abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number_temp) < 0.01 ){
//				//cout << "Triggering big number case .... " << endl;
//				// shouldn't this be i?
//				//return u;
//				return i;
//			}	else {
//				//cout << "Big number " << big_number << ", " << abs(g[u][y] - big_number) << endl;
//				cout << "i " << i << " u " << u << endl;
//				cout << "Big number -- u side " << big_number << ", " << abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number) << endl;
//				cout << " i side              " << g[z][i*RS.number_voxels_per_dim[1] + y] << endl;
//				cout << "229 Error!  Sep less than zero " << s << endl;
//				exit(1);
//			}
//		}
//	}
//
//
//	return s;
//}

//double Sep(int_type_t i, int_type_t u, int_type_t y, vector<double*>& g, int_type_t z,  ReconstructionStructure& RS, double big_number){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	//double s = int(double(u*u) - double(i*i) + pow(g[u][y], 2) - pow(g[i][y], 2))/int(2*(double(u - i)));
//
//	double s = int(double(u*u) - double(i*i) + pow(g[z][u*RS.number_voxels_per_dim[1] + y], 2) - pow(g[z][i*RS.number_voxels_per_dim[1] + y], 2))/int(2*(double(u - i)));
//
//	//	if (s < 0){
//	//		cout << "Error!  Sep less than zero " << s << endl;
//	//		exit(1);
//	//	}
//
//	//	if (s < 0){
//	//
//	//		if (fabs(s) < 0.00001){
//	//			s = 0;
//	//		}	else {
//	//			cout << "149 Error!  Sep less than zero " << s << endl;
//	//			exit(1);
//	//		}
//	//	}
//
//	if (s < 0){
//
//		if (fabs(s) < 0.00001){
//			s = 0;
//		}	else {
//			//cout << "i, u, y, gu, gi " << i << ", " << u << ", " << y << ", " << g[u][y] << ", " << g[i][y] << endl;
//			if (abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number) < 0.01){
//				//cout << "Triggering big number case .... " << endl;
//				// shouldn't this be i?
//				return u;
//			}	else {
//				//cout << "Big number " << big_number << ", " << abs(g[u][y] - big_number) << endl;
//				cout << "i " << i << " u " << u << endl;
//				cout << "Big number -- u side " << big_number << ", " << abs(g[z][u*RS.number_voxels_per_dim[1] + y] - big_number) << endl;
//				cout << " i side              " << g[z][i*RS.number_voxels_per_dim[1] + y] << endl;
//				cout << "229 Error!  Sep less than zero " << s << endl;
//				exit(1);
//			}
//		}
//	}
//
//
//	return s;
//}

//// NOte: here g is the squared distance
//double Sep(int_type_t i, int_type_t u, int_type_t x, int_type_t y, vector<vector<vector<double> > >& g){
//	if (u == i){
//		cout << "Error!  Divide by zero " << i << ", " << u << endl;
//		exit(1);
//	}
//
//	double s = int(double(u*u) - double(i*i) + g[x][y][u] - g[x][y][i])/int(2*(double(u - i)));
//
//	if (s < 0){
//		cout << "Error!  Sep less than zero " << s << endl;
//		exit(1);
//	}
//	return s;
//}


// scan 4, non-grid version
double Sep_with_square(int_type_t i, int_type_t u, double* g, int_type_t* line){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + g[line[u]] - g[line[i]])/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}

// scan 4, non-grid version
double Sep(int_type_t i, int_type_t u, double* g, int_type_t* line, int_type_t big_number){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double a= 0;
	double b = 0;
	if (line[u] != big_number){
		a = pow(g[line[u]], 2);
	}

	if (line[i] != big_number){
		b = pow(g[line[i]], 2);
	}

	double s = int(double(u*u) - double(i*i) + a - b)/int(2*(double(u - i)));

	if (s < 0){
		cout << "Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}

// NOte: here g is the squared distance
double Sep(int_type_t i, int_type_t u, int_type_t x, int_type_t y, vector<vector<vector<double> > >& g){
	if (u == i){
		cout << "Error!  Divide by zero " << i << ", " << u << endl;
		exit(1);
	}

	double s = int(double(u*u) - double(i*i) + g[x][y][u] - g[x][y][i])/int(2*(double(u - i)));

	if (s < 0){
		cout << "227 Error!  Sep less than zero " << s << endl;
		exit(1);
	}
	return s;
}




// z dominant version here.
void ReturnXYZIndicesFromIndex(int_type_t index, int_type_t& x, int_type_t& y, int_type_t& z, int_type_t xsize, int_type_t ysize){
	z = index/(xsize*ysize);
	int_type_t interim = index % (xsize*ysize);

	x = interim/ysize;
	y = interim % ysize;
}



void DistanceTransformSqdParallelTubeXIV(vector<vector< int_type_t*> >& d_gamma_indices_per_thread,
		vector<int_type_t>& count_d_gamma_per_thread, vector<uint_dist_type*>& dt_xyz, uint_dist_type big_number,
		uint_dist_type band_increment, int_type_t xsize, int_type_t ysize, int_type_t zsize, vector<bool*>& grid, int_type_t grid_xsize,
		int_type_t grid_ysize, int_type_t grid_zsize, int_type_t grid_resolution){

	//cout << "Probably need to check the g_y usage here .... see the DistanceTransform 0.  Right now have dt_xyz in scan 3-4. " << endl;
	//exit(1);
	// This is an implementation of the Meijster linear-time distance transform method, but for connected components instead of a grid
	// within each segment, the distances are going to be constrained to be within the object.  Need some scratch surfaces ....

	// v.1.1   The dt_xyz update is in the form of a x*y solice for each z step
	//big_number is -1
	//int_type_t big_number_int = big_number;
	uint_dist_type band_sq = band_increment*band_increment;
	double double_big_number = xsize*ysize*zsize/2;
	double_big_number = xsize*ysize;
	int_type_t n = max(max(xsize, ysize), zsize);
	int_type_t xi, yi, zi; int_type_t xj, yj, zj;
	//int_type_t en = edge_voxels.size();
	cout << "xsize, ysize, zsize " << xsize << ", " << ysize << ", " << zsize << endl;

	// we need s and t for each thread ....arrays are faster
	int_type_t number_threads = omp_get_max_threads();


	/// set up dt_xyz
#pragma omp parallel for
	for (int_type_t z_slice = 0; z_slice < zsize; z_slice++){
		for (int_type_t slice_i = 0; slice_i < xsize*ysize; slice_i++){
			dt_xyz[z_slice][slice_i] = big_number;
		}
	}

	//	// probably not much gain in parallelizing this.
	//	for (int_type_t i = 0, in = edge_voxels.size(); i <in; i++){
	//		ReturnXYZIndicesFromIndex(edge_voxels[i], xi, yi, zi, xsize, ysize);
	//		dt_xyz[zi][xi*ysize + yi] = 0;
	//	}

	int_type_t number_panels;
	int_type_t index_this_panel;
	int_type_t max_per_panel = xsize*ysize;  // index this panel is greater than, increments afterwards.
	int_type_t gx, gy, gz;

#pragma omp parallel for private(number_panels, index_this_panel, xi, yi, zi)
	for (int_type_t thread_id= 0; thread_id < number_threads; thread_id++){
		if (count_d_gamma_per_thread[thread_id] > 0){
			number_panels = count_d_gamma_per_thread[thread_id]/max_per_panel;
			for (int_type_t j = 0; j <= number_panels; j++){


				if (j == number_panels){
					// find out the index for this panel
					index_this_panel = count_d_gamma_per_thread[thread_id] % max_per_panel;

					if (index_this_panel == 0){ // cannot have 0-indexed, would be 1 because index is stopping criterion.
						index_this_panel = max_per_panel;
					}
				}	else {
					index_this_panel = max_per_panel;
				}

				//cout << "Distance transform, thread " << thread_id << " panel " << j << " and number indexes this panel " << index_this_panel << endl;
				//cout << "count " << count_d_gamma_per_thread[thread_id] << endl;
				for (int_type_t k = 0; k < index_this_panel; k++){
					ReturnXYZIndicesFromIndex(d_gamma_indices_per_thread[thread_id][j][k], xi, yi, zi, xsize, ysize);
					dt_xyz[zi][xi*ysize + yi] = 0;
				}
			}
		}
	}




	vector< int_type_t* > s(number_threads, 0);
	vector< int_type_t* > t(number_threads, 0);
	vector< bool* > occupied_line(number_threads, 0);
	//vector< int_type_t* > current_line(number_threads, 0);
	vector< int_type_t* > starts(number_threads, 0);
	vector< int_type_t* > stops(number_threads, 0);
	vector< double* > dt_line(number_threads, 0);
	vector< double* > g_xy_line(number_threads, 0);


	for (int_type_t i = 0; i < number_threads; i++){
		s[i] = new int_type_t[max(xsize, zsize)];
		t[i] = new int_type_t[max(xsize, zsize)];
		dt_line[i] = new double[max(xsize, zsize)];
		g_xy_line[i] = new double[max(xsize, zsize)];
		//current_line[i] = new int_type_t[n + 2];
		starts[i] = new int_type_t[n + 2];
		stops[i] = new int_type_t[n + 2];
		occupied_line[i] = new bool[max(xsize, zsize)];
	}
	int q; int w;

	//#pragma omp parallel for
	//	for (int_type_t z_slice = 0; z_slice < zsize; z_slice++){
	//		for (int_type_t slice_i = 0; slice_i < xsize*ysize; slice_i++){
	//			g_xy[z_slice][slice_i] = dt_xyz[z_slice][slice_i];
	//		}
	//	}


	int_type_t thread_id = 0;

	// first, go fishing.  Do scan 1 and 2 where x is static, z static, y increases and decreases.
	// scan 1

	// simplified scan 1 and 2.. Could be better, especially for sparse data ....
	cout << "Start scan 0" << endl;

	/// TODO -- instead of walking through everyone, use the edge voxels. ESPECIALLY WITH THE MASK.;;

	/// to try, use the mask, and then just mark down plus/minus within band_increment of the mask.
	int_type_t local_start, local_stop;



	int_type_t limits_y;
	int_type_t limits_top_y;
	/// new version. -- each thread does a z.
#pragma omp parallel for private (gx, gz,limits_y, limits_top_y )  /// if z should be parallelized ....
	for (int_type_t z = 0; z < zsize; z++){
		gz = z/grid_resolution;
		for (int_type_t x= 0; x < xsize; x++){
			gx = x/grid_resolution;

			/// Scan 0
			for (int_type_t gy0 = 0; gy0 < grid_ysize; gy0++){
				if (grid[gz][gx*grid_ysize + gy0] == true){
					if (gy0 == grid_ysize - 1){
						limits_y = ysize;
					}	else {
						// to next block
						limits_y = (gy0 + 1)*grid_resolution + 1;
					}

					for (int_type_t y= gy0*grid_resolution + 1; y < limits_y; y++){
						// minimum of what is there, and adding one to the previous one.
						if (dt_xyz[z][x*ysize + y - 1] < band_increment - 1){
							dt_xyz[z][x*ysize + y] = min(uint_dist_type(dt_xyz[z][x*ysize + y - 1] + 1), dt_xyz[z][x*ysize + y]);
						}
					}
				}

			}

			// opposite direction .... scan 1
			for (int_type_t gy0 = grid_ysize; gy0 > 0; gy0--){
				/// true gy index is gy0 - 1
				if (grid[gz][gx*grid_ysize + (gy0 - 1)] == true){
					if (gy0 - 1 == 0){
						limits_y = 1;
					}	else {
						// from the top of the block to the bottom of the block -- and then step over by 1..
						limits_y = (gy0 -1)*grid_resolution;
					}

					if (gy0 == grid_ysize){
						limits_top_y = ysize;

					}	else {
						limits_top_y = gy0*grid_resolution; // top of the block -- starting that the one just inside.
					}

					for (int_type_t y= limits_top_y; y > limits_y; y--){
						// minimum of what is there, and adding one to the previous one.
						if (dt_xyz[z][x*ysize + y - 1] < band_increment - 1){
							dt_xyz[z][x*ysize + y - 2] = min(uint_dist_type(dt_xyz[z][x*ysize + y - 1] + 1), dt_xyz[z][x*ysize + y - 2]);
						}
					}
				}
			}

		}
	}

	cout << "Start scan 2" << endl;

	//??? if at this step and the z step, we can only process those that are within band size.
	// TODO start here

	int current_distance; bool walk_condition;

	int_type_t line_counter = 0;


	// scan 3 and 4 -- y and z are static, x is moving.
	int_type_t j_index;
	int_type_t stop_index;


	cout << "Start scan 3" << endl;

	bool some_found = false;
	int_type_t number_starts, number_stops;
	double temporary_double;

	int_type_t relevant_gs = 0;

	int_type_t xstart, ystart, zstart, yend, xend, zend;

	for (int_type_t gz = 0; gz < grid_zsize; gz++){
		for (int_type_t gy = 0; gy < grid_ysize; gy++){
			number_starts = 0; number_stops = 0;

			if (grid[gz][0*grid_ysize + gy] == true){
				starts[0][number_starts] = 0; number_starts++;
			}

			/// This scan is  x -- so find all of the starts and stops.
			for (int_type_t gx = 0; gx < grid_xsize; gx++){

				if (grid[gz][gx*grid_ysize + gy] == true){
					if (gx  > 0){
						if (grid[gz][(gx - 1)*grid_ysize + gy] == false){
							starts[0][number_starts] = gx*grid_resolution; number_starts++;
						}
					}

					if (gx < (grid_xsize - 1)){
						if (grid[gz][(gx + 1)*grid_ysize + gy] == false){
							stops[0][number_stops] = (gx + 1)*grid_resolution - 1; number_stops++;
						}
					}
				}

			}


			/// tidy up, stop.
			if (grid[gz][(grid_xsize - 1)*grid_ysize + gy] == true){
				stops[0][number_stops] = xsize - 1; number_stops++;
			}

			/////////////////////////////////// now go through the stops and start passages. //////////////, need to do full z, y, x.
			if (number_starts != number_stops){
				cout << "gy, gz, " << gy << ", " << gz << endl;
				cout << "Start and stop size not the same " << number_starts << ", " << number_stops << endl;
				exit(1);
			}




			/// this same pattern for z, y, x.
			gy == grid_ysize - 1 ? yend = ysize : yend = (gy + 1)*grid_resolution;
			gz == grid_zsize - 1 ? zend = zsize : zend = (gz + 1)*grid_resolution;
			ystart = (gy)*grid_resolution;
			zstart = (gz)*grid_resolution;
			if (number_starts > 0){
				for (int_type_t z = zstart; z < zend; z++){
#pragma omp parallel for private(q, thread_id, w, local_start, local_stop, temporary_double)
					for (int_type_t y = ystart; y < yend; y++){
						thread_id = omp_get_thread_num();
						int_type_t line_counter = 0;
						/// do the parallel loop here, this is the true number this iteration.
						for (int_type_t ss = 0; ss < number_starts; ss++ ){
							//some_found = true;
							/// try one line
							//				local_start = starts[thread_id][ss];
							//				local_stop = stops[thread_id][ss];

							local_start = starts[0][ss];
							local_stop = stops[0][ss];

							// current strategy -- local start and local stop are the limits-- relative terms,
							// 0 to local stop - local start + 1
							// copy over dt_xyz, then at conclusion, copy back.
							//line_counter = local_stop - local_start + 1;

							line_counter = 0;
							for (int_type_t abs_u = local_start; abs_u < local_stop + 1; abs_u++, line_counter++){
								if (dt_xyz[z][abs_u*ysize + y] > band_sq){
									dt_line[thread_id][line_counter] = double_big_number;
									//cout << "dt " << dt_line[thread_id][line_counter] << endl;
								}	else {
									dt_line[thread_id][line_counter] = double(dt_xyz[z][abs_u*ysize + y]);
									//cout << "dt " << dt_line[thread_id][line_counter] << endl;
								}
							}
							//cout << "Total elements " << line_counter << endl;

							if (line_counter != local_stop - local_start + 1){
								cout << "Error line counter ..." << line_counter << endl;
								cout << local_stop - local_start + 1 << endl;
								exit(1);
							}

							// get rid of current line stuff.

							q = 0; s[thread_id][0]= 0; t[thread_id][0] = 0;

							// scan 3
							for (int_type_t u = 1; u < line_counter; u++){

								//double f(int_type_t x, int_type_t i, vector<double*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){
								//while (q >= 0 && f(t[thread_id][q], s[thread_id][q], dt_xyz, current_line[thread_id], y, z, ysize) > f(t[thread_id][q], u, dt_xyz, current_line[thread_id], y, z, ysize)){
								while (q >= 0 && f(t[thread_id][q], s[thread_id][q], dt_line[thread_id]) > f(t[thread_id][q], u, dt_line[thread_id])){
									q--;
								}

								if (q < 0){
									//cout << "q is zero " << endl;
									q = 0;  s[thread_id][0] = u;
								} else {

									//w = 1 + Sep(s[thread_id][q], u, dt_xyz, current_line[thread_id], y, z, ysize, big_number);
									w = 1 + Sep(s[thread_id][q], u, dt_line[thread_id], double_big_number);


									if (w < 0){
										cout << "EEEO! w less than zero " << endl;
										exit(1);
									}
									if (w < int(line_counter)){
										q++; s[thread_id][q] = u; t[thread_id][q] = w;
									}
								}

								//char ch; cin >> ch;
							}

							//					for (int_type_t u = 0; u < line_counter; u++){
							//						cout << "u, s[u], t[u] " << u << ", " << s[thread_id][u] << ", " << t[thread_id][u] <<endl;
							//					}

							// scan 4
							for (int_type_t u = line_counter; u > 0; u--){

								//if (current_line[thread_id][u-1] != big_number_int){
								//if (current_line[thread_id][u-1] != big_number_int){

								//g_xy[z][(current_line[thread_id][u-1])*ysize + y] = f(u - 1, s[thread_id][q], dt_xyz, current_line[thread_id], y, z, ysize);
								// copy back
								// TODO improve
								temporary_double = f(u - 1, s[thread_id][q], dt_line[thread_id]);
								//cout << "g_x " << u - 1 << " val " << temporary_double << endl;
								if (temporary_double <= band_sq){
									dt_xyz[z][(local_start + u - 1)*ysize + y] = uint_dist_type(temporary_double); /// should be rounded already
									//relevant_gs++;
								}	else {
									dt_xyz[z][(local_start + u - 1)*ysize + y] = big_number; // should be not much bigger than band size^2 + 1
								}

								//}

								if ((u-1) == t[thread_id][q]){
									q--;
								}
							}
						}

					}
				}
			}
		}
	}



//	/// can we do an abbreviated scan, only for the relevant regions?  Need to read again more on the Meijster algo, and what Sep, s, t, f represent.  Then, it will be more clear.
//	// scan 3 and 4
//	// scan 3 and 4 -- y and z are static, x is moving.
//	for (int_type_t z = 0; z < zsize; z++){
//		//for (int_type_t z = 29; z < 30; z++){
//#pragma omp parallel for private(q, thread_id, w, number_starts, number_stops, local_start, local_stop, temporary_double, some_found, gx, gy, gz)
//		for (int_type_t y= 0; y < ysize; y++){
//
//			thread_id = omp_get_thread_num();
//
//
//			gz = z/grid_resolution;
//			gy = y/grid_resolution;
//			//cout << "thread id " << thread_id << endl;
//
//			//vector<int> starts;
//			//vector<int> stops; /// redo this once it works. --- make into arrs like the current lines.
//
//			number_starts = 0; number_stops = 0;
//			some_found = false;
//
//
//			/// first, determine occupiededness, from a band size perspective.
//			for (int_type_t u = 0; u < xsize; u++){
//				occupied_line[thread_id][u] = false;
//			}
//
//
//			/// only look at this block if the grid says so ....
//			for (int_type_t u = 0; u < xsize; u++){
//				gx = u/grid_resolution;
//				if (grid[gz][gx*grid_ysize + gy] == false){
//					u = u + grid_resolution - 1;  /// end of loop will update by 1.
////					if (u % grid_resolution != (grid_resolution - 1)){
////						cout << "Some error with computation here ... " << endl;
////						cout << "u " << u << endl;
////						cout << "grid " << grid_resolution << endl;
////						exit(1);
////					}
//				}	else {
//					if (dt_xyz[z][u*ysize + y] <= band_increment){ /// currently only have the zs, so go off of that.
//						occupied_line[thread_id][u] = true;
//						some_found = true;
//
//						if (u != 0){
//							if (occupied_line[thread_id][u - 1] == false ){
//								// calculate how many back we need to write ...
//								// just write a square.
//
//								u > band_increment ? stop_index = u - band_increment : stop_index = 0;
//
//								for (int_type_t u1 = u; u1 > stop_index; u1--){
//									occupied_line[thread_id][u1 - 1] = true;
//								}
//							}
//						}
//
//
//						{
//							if (u != xsize - 1){
//								if (dt_xyz[z][(u + 1)*ysize + y] > band_increment ){
//									// calculate how many back we need to write ...
//									// just write a square.
//									u + band_increment < xsize ? stop_index = u + band_increment: stop_index = xsize;
//
//									for (int_type_t u1 = u + 1; u1 < stop_index; u1++){
//										occupied_line[thread_id][u1] = true;
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//
//			if (some_found){
//				//cout << "Some found, z " << z << " y " << y << endl;
//
//				/// now compute the intervals over the line for which we need to compute distances ....
//				if (occupied_line[thread_id][0] == true){
//					//starts.push_back(0);
//					starts[thread_id][number_starts] = 0; number_starts++;
//
//				}
//				for (int_type_t u = 1; u < xsize - 1; u++){
//
//					if (occupied_line[thread_id][u-1] == false && occupied_line[thread_id][u] == true ){
//						//starts.push_back(u);
//						starts[thread_id][number_starts] = u; number_starts++;
//					}
//
//					/// do the same with stops.
//					if (occupied_line[thread_id][u+1] == false && occupied_line[thread_id][u] == true ){
//						stops[thread_id][number_stops] = u; number_stops++;
//					}
//				}
//
//				//clean up -- stops and starts
//				if (occupied_line[thread_id][xsize - 2] == false && occupied_line[thread_id][xsize - 1] == true ){
//					//starts.push_back(xsize - 1);
//					starts[thread_id][number_starts] = xsize - 1;  number_starts++;
//				}
//
//				if (occupied_line[thread_id][xsize - 1] == true ){
//					//stops.push_back(xsize - 1);
//					stops[thread_id][number_stops] = xsize - 1;  number_stops++;
//				}
//
//				//if (starts.size() != stops.size()){
//				if (number_starts != number_stops){
//					cout << "y, z, " << y << ", " << z << endl;
//					cout << "Start and stop size not the same " << number_starts << ", " << number_stops << endl;
//					exit(1);
//				}
//
//				int_type_t line_counter = 0;
//				// scan 3 and 4 -- y and z are static, x is moving.
//
//				//			if (number_starts > 0){
//				//				cout << "Number intervals " << number_starts << endl;
//				//			}
//				//cout << "Starts size " << starts.size() << endl;
//				for (int_type_t ss = 0, sn = number_starts; ss < sn; ss++ ){
//					//some_found = true;
//					/// try one line
//					local_start = starts[thread_id][ss];
//					local_stop = stops[thread_id][ss];
//
//					// current strategy -- local start and local stop are the limits-- relative terms,
//					// 0 to local stop - local start + 1
//					// copy over dt_xyz, then at conclusion, copy back.
//					//line_counter = local_stop - local_start + 1;
//
//					line_counter = 0;
//					for (int_type_t abs_u = local_start; abs_u < local_stop + 1; abs_u++, line_counter++){
//						if (dt_xyz[z][abs_u*ysize + y] > band_sq){
//							dt_line[thread_id][line_counter] = double_big_number;
//							//cout << "dt " << dt_line[thread_id][line_counter] << endl;
//						}	else {
//							dt_line[thread_id][line_counter] = double(dt_xyz[z][abs_u*ysize + y]);
//							//cout << "dt " << dt_line[thread_id][line_counter] << endl;
//						}
//					}
//					//cout << "Total elements " << line_counter << endl;
//
//					if (line_counter != local_stop - local_start + 1){
//						cout << "Error line counter ..." << line_counter << endl;
//						cout << local_stop - local_start + 1 << endl;
//						exit(1);
//					}
//
//					// get rid of current line stuff.
//
//					q = 0; s[thread_id][0]= 0; t[thread_id][0] = 0;
//
//					// scan 3
//					for (int_type_t u = 1; u < line_counter; u++){
//
//						//double f(int_type_t x, int_type_t i, vector<double*>& g, int_type_t* x_pos_in_line, int_type_t y, int_type_t z, int_type_t ysize){
//						//while (q >= 0 && f(t[thread_id][q], s[thread_id][q], dt_xyz, current_line[thread_id], y, z, ysize) > f(t[thread_id][q], u, dt_xyz, current_line[thread_id], y, z, ysize)){
//						while (q >= 0 && f(t[thread_id][q], s[thread_id][q], dt_line[thread_id]) > f(t[thread_id][q], u, dt_line[thread_id])){
//							q--;
//						}
//
//						if (q < 0){
//							//cout << "q is zero " << endl;
//							q = 0;  s[thread_id][0] = u;
//						} else {
//
//							//w = 1 + Sep(s[thread_id][q], u, dt_xyz, current_line[thread_id], y, z, ysize, big_number);
//							w = 1 + Sep(s[thread_id][q], u, dt_line[thread_id], double_big_number);
//
//
//							if (w < 0){
//								cout << "EEEO! w less than zero " << endl;
//								exit(1);
//							}
//							if (w < int(line_counter)){
//								q++; s[thread_id][q] = u; t[thread_id][q] = w;
//							}
//						}
//
//						//char ch; cin >> ch;
//					}
//
//					//					for (int_type_t u = 0; u < line_counter; u++){
//					//						cout << "u, s[u], t[u] " << u << ", " << s[thread_id][u] << ", " << t[thread_id][u] <<endl;
//					//					}
//
//					// scan 4
//					for (int_type_t u = line_counter; u > 0; u--){
//
//						//if (current_line[thread_id][u-1] != big_number_int){
//						//if (current_line[thread_id][u-1] != big_number_int){
//
//						//g_xy[z][(current_line[thread_id][u-1])*ysize + y] = f(u - 1, s[thread_id][q], dt_xyz, current_line[thread_id], y, z, ysize);
//						// copy back
//						// TODO improve
//						temporary_double = f(u - 1, s[thread_id][q], dt_line[thread_id]);
//						//cout << "g_x " << u - 1 << " val " << temporary_double << endl;
//						if (temporary_double <= band_sq){
//							dt_xyz[z][(local_start + u - 1)*ysize + y] = uint_dist_type(temporary_double); /// should be rounded already
//							//relevant_gs++;
//						}	else {
//							dt_xyz[z][(local_start + u - 1)*ysize + y] = big_number; // should be not much bigger than band size^2 + 1
//						}
//
//						//}
//
//						if ((u-1) == t[thread_id][q]){
//							q--;
//						}
//					}
//				}
//			}
//		}
//	}

	//cout << "Number relevant marked in gmaps " << relevant_gs << endl;


	//		/// need to deal with this., too.
	//		//double band_increment_root2 = band_increment*sqrt(2);
	//		// debugging .....
	//			for (int_type_t z_slice = 0; z_slice < zsize; z_slice++){
	//				for (int_type_t slice_i = 0; slice_i < xsize*ysize; slice_i++){
	//					/// sqrts expensive, only do them if needed.
	//					//if (g_xy[z_slice][slice_i] != big_number){
	//					dt_xyz[z_slice][slice_i] = (g_xy[z_slice][slice_i]);
	//					//}
	//					//			if (g_xy[z_slice][slice_i] <= band_increment_sq){
	//					//				dt_xyz[z_slice][slice_i] = sqrt(g_xy[z_slice][slice_i]);
	//					//			}	else {
	//					//				//dt_xyz[z_slice][slice_i] = big_number; /// initialize again?
	//					//				g_xy[z_slice][slice_i] = big_number; /// initialize again for next round.
	//					//
	//					//			}
	//				}
	//			}


	cout << "End scan 4" << endl;
	if (zsize ==1){

		//		for (int_type_t z_slice = 0; z_slice < zsize; z_slice++){
		//			for (int_type_t slice_i = 0; slice_i < xsize*ysize; slice_i++){
		//				dt_xyz[z_slice][slice_i] = sqrt(g_xy[z_slice][slice_i]);
		//			}
		//		}
	}	else {


		// // scan 5 and 6 -- x and y are static, z is moving.
		for (int_type_t gx = 0; gx < grid_xsize; gx++){
			for (int_type_t gy = 0; gy < grid_ysize; gy++){
				number_starts = 0; number_stops = 0;

				if (grid[0][gx*grid_ysize + gy] == true){
					starts[0][number_starts] = 0; number_starts++;
				}

				/// This scan is  x -- so find all of the starts and stops.
				for (int_type_t gz = 0; gz < grid_zsize; gz++){

					if (grid[gz][gx*grid_ysize + gy] == true){
						if (gz  > 0){
							if (grid[gz - 1][(gx)*grid_ysize + gy] == false){
								starts[0][number_starts] = gz*grid_resolution; number_starts++;
							}
						}

						if (gz < (grid_zsize - 1)){
							if (grid[gz + 1][(gx)*grid_ysize + gy] == false){
								stops[0][number_stops] = (gz + 1)*grid_resolution - 1; number_stops++;
							}
						}
					}

				}

				/// tidy up, stop.
				if (grid[grid_zsize - 1][gx*grid_ysize + gy] == true){
					stops[0][number_stops] = zsize - 1; number_stops++;
				}

				/////////////////////////////////// now go through the stops and start passages. //////////////, need to do full z, y, x.
				if (number_starts != number_stops){
					cout << "gy, gx, " << gy << ", " << gx << endl;
					cout << "Start and stop size not the same " << number_starts << ", " << number_stops << endl;
					exit(1);
				}




				/// this same pattern for z, y, x.
				gy == grid_ysize - 1 ? yend = ysize : yend = (gy + 1)*grid_resolution;
				gx == grid_xsize - 1 ? xend = xsize : xend = (gx + 1)*grid_resolution;
				ystart = (gy)*grid_resolution;
				xstart = (gx)*grid_resolution;
				if (number_starts > 0){

					for (int_type_t x = xstart; x < xend; x++){
#pragma omp parallel for private(q, thread_id, w, local_start, local_stop, temporary_double)
						for (int_type_t y = ystart; y < yend; y++){
							thread_id = omp_get_thread_num();

							int_type_t line_counter = 0;
							// scan 5 and 6 -- x and y are static, z is moving.

							//cout << "Starts size " << starts.size() << endl;
							for (int_type_t ss = 0, sn = number_starts; ss < sn; ss++ ){
								//some_found = true;
								/// try one line
								local_start = starts[0][ss];
								local_stop = stops[0][ss];

								//g_xy[u][x*ysize + y ... g_xy_line
								line_counter = 0;
								for (int_type_t abs_u = local_start; abs_u < local_stop + 1; abs_u++, line_counter++){
									if (dt_xyz[abs_u][x*ysize + y] > band_sq){
										g_xy_line[thread_id][line_counter] = double_big_number;
									}	else {
										g_xy_line[thread_id][line_counter]  = double(dt_xyz[abs_u][x*ysize + y]);
									}
								}

								q = 0; s[thread_id][0]= 0; t[thread_id][0] = 0;

								// scan 5
								for (int_type_t u = 1; u < line_counter; u++){

									while (q >= 0 && f_with_square(t[thread_id][q], s[thread_id][q], g_xy_line[thread_id])
									> f_with_square(t[thread_id][q], u, g_xy_line[thread_id])){
										q--;
									}

									if (q < 0){
										//cout << "q is zero " << endl;
										q = 0;  s[thread_id][0] = u;
									} else {
										w = 1 + Sep_with_square(s[thread_id][q], u, g_xy_line[thread_id], double_big_number);

										if (w < 0){
											cout << "EEEO! w less than zero " << endl;
											exit(1);
										}
										if (w < int(line_counter)){
											q++; s[thread_id][q] = u; t[thread_id][q] = w;
										}
									}
								}

								// scan 6
								for (int_type_t u = line_counter; u > 0; u--){


									temporary_double = f_with_square(u - 1, s[thread_id][q], g_xy_line[thread_id]);
									if (temporary_double <= band_sq){
										dt_xyz[local_start + u - 1][x*ysize + y] = temporary_double;
									}	else {
										dt_xyz[local_start + u - 1][x*ysize + y] = big_number;
									}

									//						//dt_xyz[current_line[thread_id][u-1]][x*ysize + y] =  sqrt(f_with_square(u - 1, s[thread_id][q], x, y, g_xy, current_line[thread_id], ysize));
									//						dt_xyz[current_line[thread_id][u-1]][x*ysize + y] =  (f_with_square(u - 1, s[thread_id][q], x, y, g_xy, current_line[thread_id], ysize));
									//
									//						//}

									if ((u-1) == t[thread_id][q]){
										q--;
									}
								}
							}
						}
					}
				}
			}
		}
	}





//	cout << "Before scan 5-6" << endl;
//	for (int_type_t x = 0; x < xsize; x++){
//		// TODO
//#pragma omp parallel for private(q, thread_id, w, number_starts, number_stops, local_start, local_stop, temporary_double, some_found, gx, gy, gz)
//		//#pragma omp parallel for private(q, thread_id, w)
//		for (int_type_t y= 0; y < ysize; y++){
//			//		for (int_type_t x = 0; x < 1; x++){
//			//			for (int_type_t y= 0; y < 1; y++){
//
//			thread_id = omp_get_thread_num();
//			gx = x/grid_resolution;
//			gy = y/grid_resolution;
//
//			number_starts = 0;  number_stops = 0;
//			some_found = false;
//
//
//			/// first, determine occupiededness, from a band size perspective.
//			for (int_type_t u = 0; u < zsize; u++){
//				occupied_line[thread_id][u] = false;
//			}
//
//
//			for (int_type_t u = 0; u < zsize; u++){
//
//				gz = u/grid_resolution;
//				if (grid[gz][gx*grid_ysize + gy] == false){
//					u = u + grid_resolution - 1;  /// end of loop will update by 1.
//					//					if (u % grid_resolution != (grid_resolution - 1)){
//					//						cout << "Some error with computation here ... " << endl;
//					//						cout << "u " << u << endl;
//					//						cout << "grid " << grid_resolution << endl;
//					//						exit(1);
//					//					}
//				}	else {
//
//
//
//					if (dt_xyz[u][x*ysize + y] <= band_sq){ /// currently only have the zs, so go off of that.
//						occupied_line[thread_id][u] = true;
//						some_found = true;
//
//						if (u != 0){
//							if (occupied_line[thread_id][u - 1] == false ){
//								// calculate how many back we need to write ...
//								// just write a square.
//
//								u > band_increment ? stop_index = u - band_increment : stop_index = 0;
//
//								for (int_type_t u1 = u; u1 > stop_index; u1--){
//									occupied_line[thread_id][u1 - 1] = true;
//								}
//							}
//						}
//
//
//
//						if (u != zsize - 1){
//							if (dt_xyz[u + 1][x*ysize + y] > band_sq ){
//								// calculate how many back we need to write ...
//								// just write a square.
//								u + band_increment < zsize ? stop_index = u + band_increment: stop_index = zsize;
//
//								for (int_type_t u1 = u + 1; u1 < stop_index; u1++){
//									occupied_line[thread_id][u1] = true;
//								}
//							}
//						}
//
//					}
//				}
//			}
//
//
//			/// only do a second pass if some found
//			if (some_found){
//
//				/// now compute the intervals over the line for which we need to compute distances ....
//				if (occupied_line[thread_id][0] == true){
//					//starts.push_back(0);
//					starts[thread_id][number_starts] = 0; number_starts++;
//				}
//				for (int_type_t u = 1; u < zsize - 1; u++){
//
//					if (occupied_line[thread_id][u-1] == false && occupied_line[thread_id][u] == true ){
//						//starts.push_back(u);
//						starts[thread_id][number_starts] = u; number_starts++;
//					}
//
//					/// do the same with stops.
//					if (occupied_line[thread_id][u+1] == false && occupied_line[thread_id][u] == true ){
//						//stops.push_back(u);
//						stops[thread_id][number_stops] = u; number_stops++;
//					}
//				}
//
//				//clean up -- stops and starts
//				if (occupied_line[thread_id][zsize - 2] == false && occupied_line[thread_id][zsize - 1] == true ){
//					//starts.push_back(zsize - 1);
//					starts[thread_id][number_starts] = zsize - 1; number_starts++;
//				}
//
//				if (occupied_line[thread_id][zsize - 1] == true ){
//					//stops.push_back(zsize - 1);
//					stops[thread_id][number_stops] = zsize - 1; number_stops++;
//				}
//
//				//if (starts.size() != stops.size()){
//				if (number_starts != number_stops){
//					//cout << "y, z, " << y << ", " << z << endl;
//					cout << "Start and stop size not the same --- Z " << starts.size() << ", " << stops.size() << endl;
//					exit(1);
//				}
//			}
//
//
//
//
//			int_type_t line_counter = 0;
//			// scan 5 and 6 -- x and y are static, z is moving.
//
//			//cout << "Starts size " << starts.size() << endl;
//			for (int_type_t ss = 0, sn = number_starts; ss < sn; ss++ ){
//				//some_found = true;
//				/// try one line
//				local_start = starts[thread_id][ss];
//				local_stop = stops[thread_id][ss];
//
//				//g_xy[u][x*ysize + y ... g_xy_line
//				line_counter = 0;
//				for (int_type_t abs_u = local_start; abs_u < local_stop + 1; abs_u++, line_counter++){
//					if (dt_xyz[abs_u][x*ysize + y] > band_sq){
//						g_xy_line[thread_id][line_counter] = double_big_number;
//					}	else {
//						g_xy_line[thread_id][line_counter]  = double(dt_xyz[abs_u][x*ysize + y]);
//					}
//				}
//
//				q = 0; s[thread_id][0]= 0; t[thread_id][0] = 0;
//
//				// scan 5
//				for (int_type_t u = 1; u < line_counter; u++){
//
//					while (q >= 0 && f_with_square(t[thread_id][q], s[thread_id][q], g_xy_line[thread_id])
//					> f_with_square(t[thread_id][q], u, g_xy_line[thread_id])){
//						q--;
//					}
//
//					if (q < 0){
//						//cout << "q is zero " << endl;
//						q = 0;  s[thread_id][0] = u;
//					} else {
//						w = 1 + Sep_with_square(s[thread_id][q], u, g_xy_line[thread_id], double_big_number);
//
//						if (w < 0){
//							cout << "EEEO! w less than zero " << endl;
//							exit(1);
//						}
//						if (w < int(line_counter)){
//							q++; s[thread_id][q] = u; t[thread_id][q] = w;
//						}
//					}
//				}
//
//				// scan 6
//				for (int_type_t u = line_counter; u > 0; u--){
//
//
//					temporary_double = f_with_square(u - 1, s[thread_id][q], g_xy_line[thread_id]);
//					if (temporary_double <= band_sq){
//						dt_xyz[local_start + u - 1][x*ysize + y] = temporary_double;
//					}	else {
//						dt_xyz[local_start + u - 1][x*ysize + y] = big_number;
//					}
//
//					//						//dt_xyz[current_line[thread_id][u-1]][x*ysize + y] =  sqrt(f_with_square(u - 1, s[thread_id][q], x, y, g_xy, current_line[thread_id], ysize));
//					//						dt_xyz[current_line[thread_id][u-1]][x*ysize + y] =  (f_with_square(u - 1, s[thread_id][q], x, y, g_xy, current_line[thread_id], ysize));
//					//
//					//						//}
//
//					if ((u-1) == t[thread_id][q]){
//						q--;
//					}
//				}
//			}
//		}
//	}
//}

	cout << "Complete scan " << endl;

	for (int_type_t i = 0; i < number_threads; i++){
		delete [] s[i];
		delete [] t[i];
		delete [] dt_line[i];
		delete [] g_xy_line[i];
		//delete [] current_line[i];
		delete [] occupied_line[i];
		delete [] starts[i];
		delete [] stops[i];
	}

}




/// need to rewrite for a new purpose.  Revise on grid, then move onto more complex



