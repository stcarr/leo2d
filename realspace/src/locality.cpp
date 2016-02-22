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
#include <matkit.h>
#include <filtlan.h>
//#include <spmatrix.h>

//#include <petscksp.h>
//#include <slepceps.h>
//#include <petscmatlab.h>

#ifdef USE_NAMESPACE
using namespace MATKIT;
#endif



Locality::Locality(std::vector<Sdata> sdata_in,std::vector<double> heights_in,std::vector<double> angles_in) {

	sdata = sdata_in;
	heights = heights_in;
	angles = angles_in;
	num_eigs = 15;
	max_num_local_jobs = 4;
    
}

Locality::Locality(const Locality& orig) {

}

Locality::~Locality() {

}

void Locality::setup() {

}

void Locality::initMPI(int argc, char** argv){
	
//	char help[] = "what this program does in brief can go here.";

//	PetscInitialize(&argc,&argv,(char*)0,help);
//	SlepcInitialize(&argc,&argv,(char*)0,help);
 	
	MPI_Init(&argc,&argv);

        root = 0;
	print_rank = 1;
	
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
}

void Locality::constructGeom(){
	
	if (rank == print_rank)
		printf("rank %d entering constructGeom(). \n", rank);

	int max_pairs;
	
	int* nnz;
	
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
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
		
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
		
		
	}
	
	if (rank != root){
		
		printf("rank %d allocating memory in constructGeom(). \n", rank);
	
		// Allocate memory to receive pair and "index to grid" information
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
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

	
	// Some C code follows to allocate memory for our completed arrays
	// Perhaps one can get MPI to take std::vector as a valid data type instead?
	
	if (rank == print_rank)
		printf("rank %d attempting to allocate memory for its global variables. \n",rank);
	
	printf("rank %d has max index = %d (pre allocation) \n", rank, max_index);
	int index_to_grid[max_index*4];
	
	printf("rank %d has max_index = %d \n", rank, max_index);

	for (int k = 0; k < max_index; ++k){
		//if (rank == print_rank)
		//	printf("attempting k = %d ... \n", k);
		index_to_grid[k*4 + 0] = index_to_grid_i[k];
		index_to_grid[k*4 + 1] = index_to_grid_j[k];
		index_to_grid[k*4 + 2] = index_to_grid_l[k];
		index_to_grid[k*4 + 3] = index_to_grid_s[k];
		//if (rank == print_rank)
		//	printf("success for k = %d ! \n", k);
	}
	
	int inter_pairs[2*max_inter_pairs];
	
	for (int x = 0; x < max_inter_pairs; ++x){
		inter_pairs[x*2 + 0] = inter_pairs_i[x];
		inter_pairs[x*2 + 1] = inter_pairs_j[x];
	
	}
	
	int intra_pairs[2*max_intra_pairs];

	for (int x = 0; x < max_intra_pairs; ++x){
		intra_pairs[x*2 + 0] = intra_pairs_i[x];
		intra_pairs[x*2 + 1] = intra_pairs_j[x];
	}
	
	double index_to_pos[3*max_index];
	
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
	
	constructMatrix(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t,nnz);

}

void Locality::constructMatrix(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, int* nnz){ 
	if (rank == print_rank)
		printf("rank %d entering constructMatrix(). \n", rank);
		
	if (rank == root) {
		rootMatrixSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t,nnz);
	} else {
		workerMatrixSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t,nnz);
	}

}

void Locality::rootMatrixSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, int* nnz) {

	int nShifts = 2;	
	int maxJobs = nShifts*nShifts;
	int currentJob = 0;
	double work[nShifts*nShifts][2];
	double result[num_eigs];
	
	for (int i = 0; i < nShifts; ++i){
		for (int j = 0; j < nShifts; ++j){
			double x = 1.0/(double) i;
			double y = 1.0/(double) j;
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
		
		MPI::COMM_WORLD.Recv(	
					result,					// get result from worker
					num_eigs,
					MPI::DOUBLE,
					MPI::ANY_SOURCE,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
		
		MPI::COMM_WORLD.Send(	
					work[currentJob],
					2,
					MPI::DOUBLE,
					status.Get_source(),		// send to worker that just completed
					WORKTAG);
		
		++currentJob;						// one more job sent out!
		printf("rank %d has sent job to worker rank %d. \n", rank, status.Get_source());
	}
	
	// Receive final work
	
	for (int r = 1; r < size; ++r){
		
		MPI::COMM_WORLD.Recv(	
					result,					// get result from worker
					num_eigs,
					MPI::DOUBLE,
					r,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
	
		printf("rank %d has received final work from rank %d. \n", rank, r);
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
		printf("rank %d sent STOPTAG succesfuly. \n", rank);
	}
	
}

void Locality::workerMatrixSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, int* nnz) {
	
	double result[num_eigs];
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
			printf("rank %d recieved STOPTAG. \n", rank);
			return;
		}

		printf("rank %d recieved a shift job! \n", rank);		
		for (int i = 0; i < num_eigs; ++i){
			result[i] = 12;
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
					//printf("rank %d added intra_pair for index %d: [%d,%d] = %f \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
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
					double x = index_to_pos[new_k*3 + 0] - index_to_pos[k*3 + 0];
					double y = index_to_pos[new_k*3 + 1] - index_to_pos[k*3 + 1];
					int orbit1 = index_to_grid[k*4 + 2];
					int orbit2 = index_to_grid[new_k*4 + 2];
					double theta = angles[index_to_grid[new_k*4 + 3]] - angles[index_to_grid[k*4 + 3]];					

					v[input_counter] = inter_graphene(x, y, orbit1, orbit2, theta);
					//printf("rank %d added inter_pair for index %d: [%d, %d] \n", rank, k, inter_pairs[inter_counter*2 + 0], inter_pairs[inter_counter*2+1]);
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
		interval(1) = -1;  
		interval(2) = 1;

		opts.neigWanted = 5;
		opts.numIterForEigenRange = 30;
		opts.minIter = 10;
		opts.maxIter = 500;
		opts.extraIter = 100;
		opts.stride = 5;
		opts.tol = 0.0000000001;
		mkIndex baseDeg = 10; 
		mkIndex polyDeg = 10;
		opts.reorth = 0;
		opts.disp = 3;		


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

		printf("rank %d entering eigensolver... \n",rank);

		info = FilteredLanczosEigenSolver(eigs, vecs, H, interval, polyDeg, baseDeg, opts);

		MPI::COMM_WORLD.Send(	
					result,
					num_eigs,
					MPI::DOUBLE,
					root,
					0);
					
					
		if (rank == print_rank)
			printf("rank %d finished 1 job! \n", rank);
					
		}

}

void Locality::plot(){ 

}

void Locality::save(){ 

}

void Locality::finMPI(){ 
	printf("rank %d finalizing MPI. \n", rank);
	if (rank == print_rank)
		printf("rank %d finalizing MPI. \n", rank);

	//PetscErrorCode ierr;
	//ierr = PetscFinalize();CHKERRV(ierr);

	MPI_Finalize();

}

