/* 
 * File:   hstruct.cpp
 * Author: Stephen Carr
 * 
 * Created on January 13, 2016, 3:16 PM
 */

#include "locality.h"

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include <matkit.h>
#include <filtlan.h>

#ifdef USE_NAMESPACE
using namespace MATKIT;
#endif



Locality::Locality(std::vector<Sdata> sdata_in,std::vector<double> heights_in,std::vector<double> angles_in) {

	// Set all of the run-specific options for matrix construction and paramters of the solver method used.

	job_name = "HSTRUCT_TEST";

	sdata = sdata_in;
	heights = heights_in;
	angles = angles_in;
	
	nShifts = 2;
	num_eigs = 5;
	interval_start = -1;
	interval_end = 1;
	
	energy_rescale = 20.0;
	energy_shift = 0.0;
	cheb_width = 0.2;
	poly_order = 3000;
	
}

Locality::Locality(const Locality& orig) {

}

Locality::~Locality() {

}

void Locality::setup(std::string name, int shifts, int eigs, int samples,double start, double end,double e_rescale, double e_shift, double c_width, int p_order, int solver) {
	
	// Edit run-specific options for matrix constructions and paramters of the solver method (to edit settings from the constructor)

	job_name = name;
	nShifts = shifts;
	num_eigs = eigs;
	num_samples = samples;
	interval_start = start;
	interval_end = end;
	solver_type = solver;
    energy_rescale = e_rescale;
	energy_shift = e_shift;
	cheb_width = c_width;
	poly_order = p_order;
}

void Locality::initMPI(int argc, char** argv){
	
	// Start  MPI on each rank
	  	
	MPI_Init(&argc,&argv);

        root = 0;
	print_rank = 1;
	
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
}

void Locality::constructGeom(){

	// The root node (rank == 0) uses the hstruct and sheet objects to figure out the geometery and tight-binding model sparse matrix structure
	// This information is then pased to all worker nodes

	time(&constructStart);
	
	if (rank == print_rank)
		printf("rank %d entering constructGeom(). \n", rank);

	int max_pairs;
	
	int* inter_pairs_i;
	int* inter_pairs_j;
	
	int* intra_pairs_i;
	int* intra_pairs_j;
	double* intra_pairs_t;
		
	int* index_to_grid_i;
	int* index_to_grid_j;
	int* index_to_grid_l;
	int* index_to_grid_s;
	
	double* index_to_pos_x;
	double* index_to_pos_y;
	double* index_to_pos_z;
	
	if (rank == root){
	
		// Build Hstruct object
		std::vector<Sheet> sheets;
		
		for (int i = 0; i < sdata.size(); ++i){
			sheets.push_back(Sheet(sdata[i]));
		}
		
		printf("rank %d building Hstruct. \n", rank);	
		Hstruct h(sheets,angles,heights);
		
		// Broadcast "index to grid" mapping information
		max_index = h.getMaxIndex();
		
		int stationary_sheet = 0;
		int stationary_x_offset = ( sheets[stationary_sheet].getShape(1,0) - sheets[stationary_sheet].getShape(0,0) ) / 2;
		int stationary_y_offset = ( sheets[stationary_sheet].getShape(1,1) - sheets[stationary_sheet].getShape(0,1) ) / 2;
		int center_grid[4] = {stationary_x_offset,stationary_y_offset,0,stationary_sheet};
		center_index = h.gridToIndex(center_grid);
		
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(&center_index, 1, MPI_INT, root);
		
		std::vector<std::vector<int> > index_vec = h.getIndexArray();
		
		index_to_grid_i = new int[max_index];
		index_to_grid_j = new int[max_index];
		index_to_grid_l = new int[max_index];
		index_to_grid_s = new int[max_index];
		
		for (int k = 0; k < max_index; ++k){
		
			index_to_grid_i[k] = index_vec[k][0];
			index_to_grid_j[k] = index_vec[k][1];
			index_to_grid_l[k] = index_vec[k][2];
			index_to_grid_s[k] = index_vec[k][3];
		
		}
		
		
		// Construct and prepare the pairs arrays for broadcasting
		
		std::vector<std::vector<int> > inter_pairs_vec;
		h.getInterPairs(inter_pairs_vec);
		
		std::vector<int> intra_pairs_vec_i;
		std::vector<int> intra_pairs_vec_j;
		std::vector<double> intra_pairs_vec_t;
		h.getIntraPairs(intra_pairs_vec_i, intra_pairs_vec_j, intra_pairs_vec_t);
	
		max_inter_pairs = static_cast<int>(inter_pairs_vec.size());
		max_intra_pairs = static_cast<int>(intra_pairs_vec_i.size());
		
		MPI::COMM_WORLD.Bcast(&max_inter_pairs, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(&max_intra_pairs, 1, MPI_INT, root);
		
		inter_pairs_i = new int[max_inter_pairs];
		inter_pairs_j = new int[max_inter_pairs];
		
		for(int x = 0; x < max_inter_pairs; ++x){
		
			inter_pairs_i[x] = inter_pairs_vec[x][0];
			inter_pairs_j[x] = inter_pairs_vec[x][1];
			
		}
		
		intra_pairs_i = new int[max_intra_pairs];
		intra_pairs_j = new int[max_intra_pairs];
		intra_pairs_t = new double[max_intra_pairs];
		
		for(int x = 0; x < max_intra_pairs; ++x){
		
			intra_pairs_i[x] = intra_pairs_vec_i[x];
			intra_pairs_j[x] = intra_pairs_vec_j[x];
			intra_pairs_t[x] = intra_pairs_vec_t[x];
			
		}
		
		// Get and prepare the index_to_pos array for broadcasting
		
		index_to_pos_x = new double[max_index];
		index_to_pos_y = new double[max_index];
		index_to_pos_z = new double[max_index];
		
		h.getIndexToPos(index_to_pos_x,0);
		h.getIndexToPos(index_to_pos_y,1);
		h.getIndexToPos(index_to_pos_z,2);
		
		
		// POSITION DEBUG PRINT
		/*
		
		for (int k = 0; k < max_index; ++k)
			printf("%lf, %lf, %lf \n",index_to_pos_x[k],index_to_pos_y[k],index_to_pos_z[k]);
		*/
		//
		
		
		
	}
	
	if (rank != root){
		
		printf("rank %d allocating memory in constructGeom(). \n", rank);
	
		// Allocate memory to receive pair and "index to grid" information
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(&center_index, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(&max_inter_pairs, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(&max_intra_pairs, 1, MPI_INT, root);
		
		index_to_grid_i = new int[max_index];
		index_to_grid_j = new int[max_index];
		index_to_grid_l = new int[max_index];
		index_to_grid_s = new int[max_index];
		
		inter_pairs_i = new int[max_inter_pairs];
		inter_pairs_j = new int[max_inter_pairs];
		
		intra_pairs_i = new int[max_intra_pairs];
		intra_pairs_j = new int[max_intra_pairs];
		intra_pairs_t = new double[max_intra_pairs];
		
		index_to_pos_x = new double[max_index];
		index_to_pos_y = new double[max_index];
		index_to_pos_z = new double[max_index];
		
	}
		
	MPI::COMM_WORLD.Bcast(index_to_grid_i, max_index, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(index_to_grid_j, max_index, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(index_to_grid_l, max_index, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(index_to_grid_s, max_index, MPI_INT, root);
	
	MPI::COMM_WORLD.Bcast(inter_pairs_i, max_inter_pairs, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(inter_pairs_j, max_inter_pairs, MPI_INT, root);
	
	MPI::COMM_WORLD.Bcast(intra_pairs_i, max_intra_pairs, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(intra_pairs_j, max_intra_pairs, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(intra_pairs_t, max_intra_pairs, MPI_DOUBLE, root);
	
	MPI::COMM_WORLD.Bcast(index_to_pos_x, max_index, MPI_DOUBLE, root);
	MPI::COMM_WORLD.Bcast(index_to_pos_y, max_index, MPI_DOUBLE, root);
	MPI::COMM_WORLD.Bcast(index_to_pos_z, max_index, MPI_DOUBLE, root);

	if (rank == print_rank)
		printf("rank %d attempting to allocate memory for its global variables. \n",rank);
	
	int* index_to_grid = new int[max_index*4];
	
	printf("rank %d has max_index = %d \n", rank, max_index);

	for (int k = 0; k < max_index; ++k){
		//if (rank == print_rank)
			//printf("attempting k = %d ... \n", k);
		index_to_grid[k*4 + 0] = index_to_grid_i[k];
		index_to_grid[k*4 + 1] = index_to_grid_j[k];
		index_to_grid[k*4 + 2] = index_to_grid_l[k];
		index_to_grid[k*4 + 3] = index_to_grid_s[k];
		//if (rank == print_rank)
			//printf("success for k = %d ! \n", k);
	}
	
	delete index_to_grid_i;
	delete index_to_grid_j;
	delete index_to_grid_l;
	delete index_to_grid_s;
	
	
	printf("rank %d finished index_to_grid creation \n", rank);
	
	int* inter_pairs = new int[2*max_inter_pairs];
	
	for (int x = 0; x < max_inter_pairs; ++x){
		inter_pairs[x*2 + 0] = inter_pairs_i[x];
		inter_pairs[x*2 + 1] = inter_pairs_j[x];
	
	}
	
	delete inter_pairs_i;
	delete inter_pairs_j;
	
	printf("inter_pairs size = %d \n", max_inter_pairs);
	printf("rank %d finished inter_pairs creation \n", rank);

	printf("allocating memory for intra pairs on rank %d \n", rank);
	printf("intra_pairs size = %d \n", max_intra_pairs);
	int* intra_pairs = new int[2*max_intra_pairs];
	printf("allocated memory for intra pairs on rank %d \n", rank);


	for (int x = 0; x < max_intra_pairs; ++x){
		intra_pairs[x*2 + 0] = intra_pairs_i[x];
		intra_pairs[x*2 + 1] = intra_pairs_j[x];
	}
	
	printf("rank %d finished intra_pairs creation \n", rank);
	
	double* index_to_pos = new double[3*max_index];
	
	for (int k = 0; k < max_index; ++k){
		index_to_pos[k*3 + 0] = index_to_pos_x[k];
		index_to_pos[k*3 + 1] = index_to_pos_y[k];
		index_to_pos[k*3 + 2] = index_to_pos_z[k];

	}
			
	if (rank == print_rank){	
		printf("Heterostructure has %d atoms. \n", max_index);
		printf("%d entries expected from intra. \n", max_intra_pairs);
		printf("%d entries expected from inter. \n", max_inter_pairs);
	}

	time(&constructEnd);
	
	// Call the next method, which will send out jobs depending on the solver type to each worker and have them construct a matrix specific to that job.
	constructMatrix(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
	
	delete index_to_grid;
	delete index_to_pos;
	delete inter_pairs;
	delete intra_pairs;
	delete intra_pairs_t;

}

void Locality::constructMatrix(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t){ 
	time(&solveStart);
	if (rank == print_rank)
		printf("rank %d entering constructMatrix(). \n", rank);

	// 0: FILTLAN local eigenvalue solver (DOS)
	if(solver_type == 0){
		if (rank == root) {
			rootEigenSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		} else {
			workerEigenSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		}
	}

	// 1: Chebyshev polynomial sampling of eigenvalue spectrum (DOS)
	if(solver_type == 1){
		if (rank == root) {
			rootChebSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		} else {
			workerChebSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		}	
	}

	time(&solveEnd);
}

void Locality::rootEigenSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t) {
	
	int maxJobs = nShifts*nShifts;
	int currentJob = 0;
	int eigs_count;
	double work[nShifts*nShifts][2];
	std::vector<std::vector<std::vector<double> > > result_array;
	
	for (int i = 0; i < nShifts; ++i){
		for (int j = 0; j < nShifts; ++j){
			double x = (1.0/(double) nShifts)*i;
			double y = (1.0/(double) nShifts)*j;
			work[i*nShifts + j][0] = x;
			work[i*nShifts + j][1] = y;

		}
	}
	
	MPI::Status status;
	
	// try to give each worker its first job
	
	for (int r = 1; r < size; ++r) {
		if (currentJob < maxJobs) {
		
			MPI::COMM_WORLD.Send(	
						work[currentJob], 	// input buffer
						2,					// size of buffer [x,y]
						MPI::DOUBLE,		// type of buffer
						r,					// worker to recieve
						WORKTAG);			// tag as work
						
			
			++currentJob;					// one more job sent out!
		}
		
	}
	
	
	printf("rank %d has sent first batch of work... \n", rank);
	
	// Receive results and dispense new work
	
	while (currentJob < maxJobs) {
	
				
		MPI::COMM_WORLD.Recv(				// get size of incoming work
					&eigs_count,
					1,
					MPI::INT,
					MPI::ANY_SOURCE,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
		
		double result[eigs_count*2];
		
		MPI::COMM_WORLD.Recv(	
					result,					// get result from worker
					eigs_count*2,
					MPI::DOUBLE,
					status.Get_source(),
					MPI::ANY_TAG);	

		std::vector<double> temp_eigs;
		std::vector<double> temp_weights;
		std::vector<std::vector<double> > temp_result;
		for (int i = 0; i < eigs_count; ++i){
			//printf("eigenvalue check: %lf \n", result[i]);
			temp_eigs.push_back(result[i*2 + 0]);
			temp_weights.push_back(result[i*2 + 1]);
			}
			
		temp_result.push_back(temp_eigs);
		temp_result.push_back(temp_weights);
			
		result_array.push_back(temp_result);
		
		MPI::COMM_WORLD.Send(	
					work[currentJob],
					2,
					MPI::DOUBLE,
					status.Get_source(),		// send to worker that just completed
					WORKTAG);
		
		++currentJob;						// one more job sent out!
		printf("rank %d has sent job to worker rank %d. (%d/%d) \n", rank, status.Get_source(),currentJob, maxJobs + 1);
	}
	
	// Receive final work
	
	for (int r = 1; r < size; ++r){
		
		MPI::COMM_WORLD.Recv(				// get size of incoming work
					&eigs_count,
					1,
					MPI::INT,
					MPI::ANY_SOURCE,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
		
		double result[eigs_count*2];
		
		MPI::COMM_WORLD.Recv(	
					result,					// get result from worker
					eigs_count*2,
					MPI::DOUBLE,
					status.Get_source(),
					MPI::ANY_TAG);
					
		std::vector<double> temp_eigs;
		std::vector<double> temp_weights;
		std::vector<std::vector<double> > temp_result;
		
		for (int i = 0; i < eigs_count; ++i){
			//printf("eigenvalue check: %lf \n", result[i]);
			temp_eigs.push_back(result[i*2 + 0]);
			temp_weights.push_back(result[i*2 + 1]);
			}
			
		temp_result.push_back(temp_eigs);
		temp_result.push_back(temp_weights);
			
		result_array.push_back(temp_result);
	
		printf("rank %d has received final work from rank %d. \n", rank, r);
	}
	
	// Tell workers to exit workerMatrixSolve()
	
	for (int r = 1; r < size; ++r){
		//printf("rank %d sending STOPTAG to rank %d. \n", rank, r);
		double temp[2];
		temp[0] = 0.0;
		temp[1] = 0.0;
		MPI::COMM_WORLD.Send(	
					temp,
					2, 
					MPI::DOUBLE, 
					r, 
					STOPTAG);
		//printf("rank %d sent STOPTAG succesfuly. \n", rank);
	}
	
	std::ofstream outFile1;
	std::string extension1 = "_eigs.out";
	outFile1.open((job_name + extension1).c_str());
	outFile1 << "Shift x, Shift y, eigenvalues \n";

	for(int i = 0; i < maxJobs; ++i){
		outFile1 << work[i][0] << ", " << work[i][1] << ", ";
		for(int j = 0; j < result_array[i][0].size(); ++j)
			outFile1 << result_array[i][0][j] << ", ";
		outFile1 << "\n";
	}
	
	outFile1.close();	
	
	std::ofstream outFile2;
	const char* extension2 = "_weights.out";
	outFile2.open((job_name + extension2).c_str());
	outFile2 << "Shift x, Shift y, weights \n";

	for(int i = 0; i < maxJobs; ++i){
		outFile2 << work[i][0] << ", " << work[i][1] << ", ";
		for(int j = 0; j < result_array[i][1].size(); ++j)
			outFile2 << result_array[i][1][j] << ", ";
		outFile2 << "\n";
	}
	
	outFile2.close();
	/*
	for (int i = 0; i < nShifts*nShifts; ++i)
		for (int j = 0; j < result_array[i][0].size(); ++j)
			printf("%lf , %lf\n", result_array[i][0][j],result_array[i][1][j]);
	*/
}

void Locality::workerEigenSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t) {
	
	double work[2];
	MPI::Status status;

	while (1) {

		if (rank == print_rank)
			printf("rank %d waiting for new job... \n",rank);
		MPI::COMM_WORLD.Recv( 
						work, 
						2, 
						MPI::DOUBLE, 
						root, 
						MPI::ANY_TAG, 
						status);
		
		if (status.Get_tag() == STOPTAG) {
			//printf("rank %d recieved STOPTAG. \n", rank);
			return;
		}

		time_t tempStart;
		time(&tempStart);
		solverTimes.push_back(tempStart);	

		printf("rank %d recieved a shift job! \n", rank);
		
		double i2pos[max_index*3];
		for (int i = 0; i < max_index; ++i) {
		
			if (index_to_grid[i*4 + 3] == 0){
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				}
			if (index_to_grid[i*4 + 3] == 1){
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + work[0];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + work[1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				}
		}
			
		
		int inter_counter = 0;
		int intra_counter = 0;
	
		printf("rank %d trying to build Sparse MATKIT Matrix! \n",rank);
		
		mkIndex max_nnz = max_intra_pairs + max_inter_pairs;
		mkIndex* row_index = new mkIndex[max_nnz];
		Real* v = new Real[max_nnz];
		mkIndex* col_pointer = new mkIndex[max_index+1];
		int input_counter = 0;

		for (int k = 0; k < max_index; ++k){
			
			col_pointer[k] = input_counter;
			
			bool same_index1 = true;
			while(same_index1) {
				if (intra_pairs[intra_counter*2 + 0] != k) {
					same_index1 = false;
				}
				else {
					row_index[input_counter] = intra_pairs[intra_counter*2 + 1];
					v[input_counter] = intra_pairs_t[intra_counter];
					//printf("[%d,%d] = %lf (intra) \n",intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
					++input_counter;
					++intra_counter;
				}
				
			}
		
			bool same_index2 = true;
			while(same_index2) {
				if (inter_pairs[inter_counter*2 + 0] != k) {
					same_index2 = false;
				}
				
				else {
					int new_k = inter_pairs[inter_counter*2 + 1];	
					row_index[input_counter] = new_k;
					double x1 = i2pos[k*3 + 0];
					double y1 = i2pos[k*3 + 1];
					double x2 = i2pos[new_k*3 + 0];
					double y2 = i2pos[new_k*3 + 1];
					int orbit1 = index_to_grid[k*4 + 2];
					int orbit2 = index_to_grid[new_k*4 + 2];
					
					//double theta = angles[index_to_grid[new_k*4 + 3]] - angles[index_to_grid[k*4 + 3]];
					
					double theta1 = angles[index_to_grid[k*4+3]];
					double theta2 = angles[index_to_grid[new_k*4 + 3]];
					
					v[input_counter] = inter_graphene(x1, y1, x2, y2, orbit1, orbit2, theta1, theta2);
					//printf("rank %d added inter_pair for index %d: [%d, %d] \n", rank, k, inter_pairs[inter_counter*2 + 0], inter_pairs[inter_counter*2+1]);
					//printf("[%d, %d] = %lf (inter) \n", k, new_k, v[input_counter]);
					++input_counter;
					++inter_counter;
				}
				
			}
		}

		col_pointer[max_index] = input_counter;
			

		SparseMatrix H(max_index, max_index, v, row_index, col_pointer, max_nnz);

		FilteredLanczosInfo info;
		FilteredLanczosOptions opts;	
		Vector eigs;
		Matrix vecs;
		Vector interval(2);
		interval(1) = interval_start;  
		interval(2) = interval_end;

		opts.neigWanted = num_eigs;
		/*
		opts.numIterForEigenRange = 30;
		opts.minIter = 10;
		opts.maxIter = 500;
		opts.extraIter = 100;
		opts.stride = 5;
		opts.tol = 0.0000000001;
		*/

		mkIndex baseDeg = 10; 
		mkIndex polyDeg = 10;
		opts.reorth = 2;
		
		// no printout
		//opts.disp = 0;
		
		// verbose printouts
		opts.disp = 3;		

		/*
		opts.intervalOpts.intervalWeights(1) = 100;
		opts.intervalOpts.intervalWeights(2) = 1;
		opts.intervalOpts.intervalWeights(3) = 1;
		opts.intervalOpts.intervalWeights(4) = 1;
		opts.intervalOpts.intervalWeights(5) = 1;
		opts.intervalOpts.numGridPoints = 10;
		opts.intervalOpts.initialShiftStep = 0.01;
		opts.intervalOpts.maxOuterIter = 50;
		opts.intervalOpts.yBottomLine = 0.001;
		opts.intervalOpts.yRippleLimit = 100;
		opts.intervalOpts.transitionIntervalRatio = 0.6;
		opts.intervalOpts.initialPlateau = 0.1;
		opts.intervalOpts.yLimitTol = 0.0001;
		opts.intervalOpts.maxInnerIter = 30;
		*/
		

		printf("rank %d entering eigensolver... \n",rank);

		info = FilteredLanczosEigenSolver(eigs, vecs, H, interval, polyDeg, baseDeg, opts);

		int eigs_count = eigs.Length();
		double result[eigs_count*2];
		for (int i = 0; i < eigs_count; ++i) {
			//printf("%lf \n", eigs(i+1));
			result[i*2 + 0] = eigs(i+1);
			result[i*2 + 1] = vecs(center_index+1,i+1);
			
		} 
		
		time_t tempEnd;
		time(&tempEnd);
		solverTimes.push_back(tempEnd);
	
		MPI::COMM_WORLD.Send(	
					&eigs_count,
					1,
					MPI::INT,
					root,
					0);

		MPI::COMM_WORLD.Send(	
					result,
					eigs_count*2,
					MPI::DOUBLE,
					root,
					0);
					
					
		if (rank == print_rank)
			printf("rank %d finished 1 job! \n", rank);
					
		}

}

void Locality::rootChebSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t) {
	
	int maxJobs = nShifts*nShifts;
	int currentJob = 0;
	int length;
	double work[nShifts*nShifts][2];
	std::vector<std::vector<double> > result_array;
	result_array.resize(maxJobs);
	
	for (int i = 0; i < nShifts; ++i){
		for (int j = 0; j < nShifts; ++j){
			double x = (1.0/(double) nShifts)*i;
			double y = (1.0/(double) nShifts)*j;
			work[i*nShifts + j][0] = x;
			work[i*nShifts + j][1] = y;

		}
	}
	
	// TESTING CODE
	/*
	int maxJobs = 1;
	int currentJob = 0;
	int num_samples;
	double work[1][2];
	work[0][0] = 0.0;
	work[0][1] = 0.5;
	std::vector<std::vector<double> > result_array;
	*/
	
	MPI::Status status;
	
	// try to give each worker its first job
	
	for (int r = 1; r < size; ++r) {
		if (currentJob < maxJobs) {
		
			MPI::COMM_WORLD.Send(	
						work[currentJob], 	// input buffer
						2,					// size of buffer [x,y]
						MPI::DOUBLE,		// type of buffer
						r,					// worker to recieve
						currentJob+1);		// work tag to label shift value
						
			
			++currentJob;					// one more job sent out!
		}
		
	}
	
	
	printf("rank %d has sent first batch of work... \n", rank);
	
	// Receive results and dispense new work
	
	while (currentJob < maxJobs) {
	
		int jobTag;
		MPI::COMM_WORLD.Recv(				// get size of incoming work
					&length,
					1,
					MPI::INT,
					MPI::ANY_SOURCE,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
		
		double results[length];
		
		MPI::COMM_WORLD.Recv(	
					results,					// get result from worker
					length,
					MPI::DOUBLE,
					status.Get_source(),
					MPI::ANY_TAG);

		jobTag = status.Get_tag();

		std::vector<double> temp_result;
		for (int i = 0; i < length; ++i){
			temp_result.push_back(results[i]);
			}
			
		result_array[jobTag-1] = (temp_result);
		
		MPI::COMM_WORLD.Send(	
					work[currentJob],
					2,
					MPI::DOUBLE,
					status.Get_source(),		// send to worker that just completed
					currentJob+1);
		
		++currentJob;						// one more job sent out!
		printf("rank %d has sent job to worker rank %d. \n", rank, status.Get_source());
	}
	
	// Receive final work
	
	for (int r = 1; r < size; ++r){
	
		int length;
		int jobTag;
		
		MPI::COMM_WORLD.Recv(				// get size of incoming work
					&length,
					1,
					MPI::INT,
					MPI::ANY_SOURCE,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
		
		double results[num_samples];
		
		MPI::COMM_WORLD.Recv(	
					results,					// get result from worker
					length,
					MPI::DOUBLE,
					status.Get_source(),
					MPI::ANY_TAG);
					
		jobTag = status.Get_tag();
		
		std::vector<double> temp_result;
		for (int i = 0; i < num_samples; ++i){
			//printf("eigenvalue check: %lf \n", result[i]);
			temp_result.push_back(results[i]);
			}
			
		result_array[jobTag-1] = (temp_result);
	
		printf("rank %d has received final work from rank %d. \n", rank, r);
		//if (result_array.size() == maxJobs)
			//break;
	}
	
	// Tell workers to exit workerMatrixSolve()
	
	for (int r = 1; r < size; ++r){
		printf("rank %d sending STOPTAG to rank %d. \n", rank, r);
		double temp[2];
		temp[0] = 0.0;
		temp[1] = 0.0;
		MPI::COMM_WORLD.Send(	
					temp,
					2, 
					MPI::DOUBLE, 
					r, 
					STOPTAG);
		printf("rank %d sent STOPTAG successfully. \n", rank);
	}
	
	double* polys = new double[poly_order];
	
	for (int i = 0; i < poly_order; ++i){
		polys[i] = i;
	}
	
	/*
	for (int i = 0; i < maxJobs; ++i){
		printf("printing results for shift %d \n", i+1);
		printf("Energy || Density \n");
		for (int j = 0; j < result_array[i].size(); ++j){
			printf("%lf || %lf \n", sample_points[j],result_array[i][j]);
		}
	}
	*/
	
	std::ofstream outFile;
	const char* extension = ".out";
	outFile.open ((job_name + extension).c_str());
	outFile << job_name << " Chebyshev T value outputs \n";
	outFile << "Shift x, Shift y, ... polynomial orders ... \n";
	outFile << "-1, -1";
	
	for(int j = 0; j < poly_order; ++j)
		outFile << ", " << polys[j];
	outFile << "\n";
	
	for(int i = 0; i < maxJobs; ++i){
		outFile << work[i][0] << ", " << work[i][1];
		for(int j = 0; j < num_samples; ++j)
			outFile << ", " << result_array[i][j];
		outFile << "\n";
	}
	
	outFile.close();
	
	delete[] polys;
}

void Locality::workerChebSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t) {

	double work[2];
	MPI::Status status;
	
	// figure out the Chebyshev polynomial coefficients for each energy sample range
	
	// p is the maximum order of the Chebyshev polynomial
	int p = poly_order;
	
	// alpha_p is a constant used in the definition of the Chebyshev polynomial
	double alpha_p = M_PI/(p+2);
	
	// sample_width is the scaled width of the non-zero portion of the polynomial
	double sample_width = cheb_width/energy_rescale;
	
	// sample_spacing is the scaled distance between the center of each sample
	double sample_spacing = (interval_end - interval_start)/(num_samples);
	
	// Save the center-point of each sample
	double* sample_points = new double[num_samples + 1];
	for (int i = 0; i < num_samples + 1; ++i){
		sample_points[i] = ( (interval_start + sample_spacing*i) + energy_shift)/(energy_rescale);
	}
	
	// Determine "dampening coefficients" (independent of sample location)
	double* damp_coeff = new double[poly_order + 1];
	damp_coeff[0] = 1.0;
	for (int j = 0; j < poly_order + 1; ++j){
		double jd = (double) j;
		damp_coeff[j] = ((1-jd/(p+2))*sin(alpha_p)*cos(jd*alpha_p)+(1/(p+2))*cos(alpha_p)*sin(jd*alpha_p))/sin(alpha_p);
	}
	
	// ------------------------
	// TESTING POLYNOMIAL BLOCK
	// Evaluates a Chebyshev "step" function at the center of the interval
	// If std::out is saved to a file, can plot this to see how well behaved your Chebyshev polynomial is
	//
	
	/*
	int x = (int) num_samples/2;
	for (int i = 0; i < num_samples; ++i){
	
		double cheb_coeff[poly_order];
	
		double a = sample_points[x] - sample_width/2;
		double b = sample_points[x] + sample_width/2;
		
		cheb_coeff[0] = (1.0/M_PI)*(acos(a) - acos(b));
		
		for (int j = 1; j < poly_order + 1; ++j) {
			double jd = (double) j;
			cheb_coeff[j] = (2/M_PI)*(sin(jd*acos(a))-sin(jd*acos(b)))/jd;
		}
		
		double T;
		double v_i = 1;
		double H = sample_points[i] + sample_width*0.5/energy_rescale;
		
		double T_prev = 1;
		
		double T_j = H*T_prev;
		double T_next;
		
		T = v_i*(T_prev*cheb_coeff[0]*damp_coeff[0] + T_j*cheb_coeff[1]*damp_coeff[1]);
		
		for (int j = 1; j < poly_order; ++j){
			T_next = 2*H*T_j - T_prev;
			T_prev = T_j;
			T_j = T_next;
			T = T + v_i*(T_next*cheb_coeff[j+1]*damp_coeff[j+1]);
		}
		
		printf("%lf \n",T);
	}
	*/
	// END POLY TESTING BLOCK
	// ----------------------
	
	
	
	// ---------------------
	// Enter MPI worker loop
	// ---------------------
	
	while (1) {
	
		int jobTag;

		//if (rank == print_rank)
			//printf("rank %d waiting for new job... \n",rank);
			
		
		// Recv to get a new job
		
		MPI::COMM_WORLD.Recv( 
						work, 			// work determines b-shift
						2, 				// length of work
						MPI::DOUBLE, 	// data-type of work
						root, 			// must come from root
						MPI::ANY_TAG,  	// either WORKTAG or STOPTAG
						status);		// keep MPI status information
		
		
		jobTag = status.Get_tag();
		// If worker gets STOPTAG it ends this method
		if (jobTag == STOPTAG) {
			//printf("rank %d received STOPTAG. \n", rank);
			return;
		}
		

		// If not STOPTAG, start timing the solver
		time_t tempStart;
		time(&tempStart);
		solverTimes.push_back(tempStart);	
		
		printf("rank %d received a shift job! [%lf,%lf] \n", rank,work[0],work[1]);
		
		// ------------------------------------------------------
		// Determine the work-specific positions of every orbital
		// ------------------------------------------------------
		
		// i2pos = "index to position"
		double i2pos[max_index*3];
		
		for (int i = 0; i < max_index; ++i) {
		
			if (index_to_grid[i*4 + 3] == 0){
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				}
			if (index_to_grid[i*4 + 3] == 1){
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + work[0];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + work[1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				}
				
		}
		
		// -----------------------
		// Build the Sparse matrix
		// -----------------------

		// Indexes how many inter terms we have entered so far
		int inter_counter = 0;
		
		// Indexes how many intra terms we have intered so far
		int intra_counter = 0;
		
		// Total number of expected non-zero matrix elements
		mkIndex max_nnz = max_intra_pairs + max_inter_pairs;
		
		// Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
		
		// row_index tells us the row of element i
		// v tells us the value of element i
		// col_pointer tells us the start of column j (at element i = col_pointer[j]) and the end of column j (at element i = col_pointer[j] - 1)
		
		mkIndex* row_index = new mkIndex[max_nnz];
		Real* v = new Real[max_nnz];
		mkIndex* col_pointer = new mkIndex[max_index+1];
		
		// Count the current element, i.e. "i = input_counter"
		int input_counter = 0;

		// Loop through every orbital
		for (int k = 0; k < max_index; ++k){
			
			// Save starting point of column k
			col_pointer[k] = input_counter;
			
			// While we are still at the correct index in our intra_pairs list:
			bool same_index1 = true;
			while(same_index1) {
			
				// if the first index of intra_pairs changes, we stop
				if (intra_pairs[intra_counter*2 + 0] != k) {
					same_index1 = false;
				}
				
				// otherwise we save that pair into our sparse matrix format
				else {
					row_index[input_counter] = intra_pairs[intra_counter*2 + 1];
					 
					// if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
					if (row_index[input_counter] == k)
						v[input_counter] = (intra_pairs_t[intra_counter] + energy_shift)/energy_rescale;
					// Otherwise we just enter the value just with scaling
					else
						v[input_counter] = intra_pairs_t[intra_counter]/energy_rescale;
					
					//printf("rank %d added intra_pair for index %d: [%d,%d] = %f \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
					
					++input_counter;
					++intra_counter;
				}
				
			}
		
			// While we are still at the correct index in our inter_pairs list:
			bool same_index2 = true;
			while(same_index2) {
			
				// if the first index of inter_pairs changes, we stop
				if (inter_pairs[inter_counter*2 + 0] != k) {
					same_index2 = false;
				}
				
				// otherwise we save that pair into our sparse matrix format
				else {
					// get the index of the other orbital in this term
					int new_k = inter_pairs[inter_counter*2 + 1];
					row_index[input_counter] = new_k;
					
					// get the position of both orbitals
					double x1 = i2pos[k*3 + 0];
					double y1 = i2pos[k*3 + 1];
					double x2 = i2pos[new_k*3 + 0];
					double y2 = i2pos[new_k*3 + 1];
					
					// and the orbit tag in their respective unit-cell
					int orbit1 = index_to_grid[k*4 + 2];
					int orbit2 = index_to_grid[new_k*4 + 2];
					
					// and the angle of the sheet each orbital is on
					double theta1 = angles[index_to_grid[k*4 + 3]];
					double theta2 = angles[index_to_grid[new_k*4 + 3]];

					// use all this information to determine coupling energy
					// !!! Currently NOT generalized for materials other than graphene, need to do a material index call for both sheets and pass to a general "inter_coupling" method !!!
					v[input_counter] = inter_graphene(x1, y1, x2, y2, orbit1, orbit2, theta1, theta2)/energy_rescale;
					
					//printf("rank %d added inter_pair for index %d: [%d, %d] = %f \n", rank, k, inter_pairs[inter_counter*2 + 0], inter_pairs[inter_counter*2+1], v[input_counter]);
					++input_counter;
					++inter_counter;
				}
				
			}
		}

		// Save the end point + 1 of the last column
		col_pointer[max_index] = input_counter;
		
		// Construct the Sparse Tight-binding Hamiltonian matrix
		SparseMatrix H(max_index, max_index, v, row_index, col_pointer, max_nnz);
		
		// ------------------------------
		// Follwoing saves Matrix to file
		//
		// Should only be uncommented for 1-job processes, otherwise they will overwrite each other!
		
		/*
		std::ofstream outFile;
		const char* extension = "_matrix.dat";
		outFile.open ((job_name + extension).c_str());
		
		for(int i = 0; i < max_index; ++i){
			int start_index = col_pointer[i];
			int stop_index = col_pointer[i+1];
				for(int j = start_index; j < stop_index; ++j){
					outFile << row_index[j] + 1 << ", " << i + 1 << ", " << v[j] << "\n";
				}
		}
		
		outFile.close();
		*/
		
		// End Matrix Save
		// ---------------
		
		// Starting vector for Chebyshev method is a unit-vector at the center-orbital
		Vector v_i(max_index);
		v_i(center_index) = 1.0;
		
		if(rank == print_rank)
			printf("rank %d starting Chebychev solver work... \n", rank);
			
		// Chebyshev values
		double T;
		
		// Saves T values
		double* T_array = new double[poly_order];
			
		// Temporary vector for algorithm ("previous" vector T_j-1)
		Vector T_prev(max_index);
		T_prev(center_index) = 1.0;
		
		// Temporary vector for algorithm ("current" vector T_j)
		Vector T_j = H*T_prev;
		
		// Temporary vector for algorithm ("next" vector T_j+1)
		Vector T_next;
		
		// first T value is always 1
		T_array[0] = 1;	
		
		// Next one is calculated simply
		T_array[1] = T_j(center_index);

		// Now loop algorithm up to poly_order to find all T values
		for (int j = 2; j < poly_order; ++j){
		
			T_next = 2*H*T_j - T_prev;
			T_prev = T_j;
			T_j = T_next;
			T_array[j] = T_j(center_index);
			
			// print every 100 steps on print rank
			//if (rank == print_rank)
				//if (j%100 == 0)
					//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
		}
		
		/* COMMENTED OUT (we save T values now, NOT densities)
		// Find the density of eigenvalue spectrum
		double* densities = new double[num_samples];
		
		// At each sample_point, we center a sample of sample_width
		for (int i =0; i < num_samples; ++i){
		
			// Now we determine the Chebyshev coefficients for a specific sample range in this loop
			double cheb_coeff[num_samples];
			double a = sample_points[i] - sample_width/2;
			double b = sample_points[i] + sample_width/2;
		
			cheb_coeff[0] = (1.0/M_PI)*(acos(a) - acos(b));
			
			for (int j = 1; j < poly_order + 1; ++j) {
				double jd = (double) j;
				cheb_coeff[j] = (2/M_PI)*(sin(jd*acos(a))-sin(jd*acos(b)))/jd;
			}
			
			// Finally we evaluate our polynomial using the T array
			double T = 0;
			for (int j = 0; j < poly_order; ++j){
				T = T + T_array[j]*cheb_coeff[j]*damp_coeff[j];
				}
			
			densities[i] = T;
		}
		
		*/
		
		// Save time at which solver finished
		time_t tempEnd;
		time(&tempEnd);
		solverTimes.push_back(tempEnd);	
		
		int length = poly_order;
		
		// Notify root about incoming data size
		MPI::COMM_WORLD.Send(	
					&length,
					1,
					MPI::INT,
					root,
					jobTag);

		// Send root the density information
		MPI::COMM_WORLD.Send(	
					T_array,
					length,
					MPI::DOUBLE,
					root,
					jobTag);
					
		//if (rank == print_rank)
			//printf("rank %d finished 1 job! \n", rank);
		
		// Cleanup C++ allocated memory
		delete[] T_array;
		//delete[] densities;
		
		// End of while(1) means we wait for another instruction from root
		}

}

// Not implemented (plotting is done via output files through i.e. MATLAB utility scripts)
void Locality::plot(){ 

}

// Prints timing information
void Locality::save(){

	if (rank != root){
		int jobCount = solverTimes.size() / 2;	
		double avgTime = 0;
		double maxTime = 0;
		double minTime = difftime(solverTimes[1],solverTimes[0]);	
	
		for (int i = 0; i < jobCount; ++i){
			double tempTime = difftime(solverTimes[2*i + 1],solverTimes[2*i]);
			if (tempTime > maxTime)
				maxTime = tempTime;
			if (tempTime < minTime)
				minTime = tempTime;
			avgTime += tempTime/( (double) jobCount);
		}

		printf(	"=======- RANK %d TIMING -======= \n"
			"Number of Jobs: %d \n"
			"Avg Solve Time: %lf sec \n"
			"Max Solve Time: %lf sec \n"
			"Min Solve Time: %lf sec \n"
			"================================ \n",
			rank, jobCount, avgTime, maxTime, minTime);		
	}

	if (rank == print_rank){
		printf(	"=======- TIMING RESULTS -======= \n"
		        "Construct Geom: %lf sec \n"
			"Solver        : %lf sec \n"
			"================================ \n",
			difftime(constructEnd, constructStart), 
			difftime(solveEnd, solveStart));
	}
}

void Locality::finMPI(){ 
	//printf("rank %d finalizing MPI. \n", rank);
		
	MPI_Finalize();

}
