
#ifndef DISTANCETRANSFORMS_HPP_
#define DISTANCETRANSFORMS_HPP_


#include "Includes.hpp"

//void DistanceTransform1Parallel(vector<SkelGraph>& SG, double* dt_xyz, vector<int_type_t>& grid_structure, double big_number, int_type_t connectivity);

//void DistanceTransform1Parallel(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<double*>& dt_xyz, double big_number);

//void DistanceTransform0(vector<vector<vector<double> > >& distance_grid, vector<int_type_t>& grid_structure, double big_number);

//void DistanceTransform1Edges(vector<vector<vector<double> > >& distance_grid, vector<int_type_t>& grid_structure, double big_number);

//void SetTubeBoundariesBreadthFirst(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<double*>& dt_xyz, double big_number, double band_size);

//void DistanceTransform1ParallelTube(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<double*>& dt_xyz, double big_number, double band_increment);

//void DistanceTransform1ParallelTubeII(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<double*>& dt_xyz, double big_number, double band_increment);

//void DistanceTransform1ParallelTubeIII(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<double*>& dt_xyz, double big_number, double band_increment);

//void DistanceTransform1ParallelTubeIV(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<float*>& dt_xyz,
	//	vector<float*>& g_xy, float big_number, float band_increment);

//void DistanceTransform1ParallelTubeV(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<float*>& dt_xyz,
	//	vector<float*>& g_xy, float big_number, float band_increment);

//void FlipDistanceTransformInsideVolume(ReconstructionStructure& RS, vector<float*>& dt_xyz);

//void DistanceTransformSqd1ParallelTubeVI(ReconstructionStructure& RS, vector<int_type_t>& edge_voxels, vector<uint_dist_type*>& dt_xyz, vector<uint_dist_type*>& g_xy,
		//uint_dist_type big_number, uint_dist_type band_increment);


void DistanceTransformSqdParallelTubeXIV(vector<vector< int_type_t*> >& d_gamma_indices_per_thread,
		vector<int_type_t>& count_d_gamma_per_thread, vector<uint_dist_type*>& dt_xyz, uint_dist_type big_number,
		uint_dist_type band_increment, int_type_t xsize, int_type_t ysize, int_type_t zsize, vector<bool*>& grid, int_type_t grid_xsize,
		int_type_t grid_ysize, int_type_t grid_zsize, int_type_t grid_resolution);

void ReturnXYZIndicesFromIndex(int_type_t index, int_type_t& x, int_type_t& y, int_type_t& z, int_type_t xsize, int_type_t ysize);


#endif /* DISTANCETRANSFORMS_HPP_ */
