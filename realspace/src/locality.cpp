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

Locality::Locality(std::vector<Sdata> sdata_in,std::vector<double> heights_in,std::vector<double> angles_in) {

	sdata = sdata_in;
	heights = heights_in;
	angles = angles_in;
	num_eigs = 15;
    
}

Locality::Locality(const Locality& orig) {

}

Locality::~Locality() {

}

void Locality::setup() {

}

void Locality::initMPI(int argc, char** argv){

	MPI::Init(argc, argv);
	
	root = 0;
	print_rank = 1;
	
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
}

void Locality::constructGeom(){
	
	if (rank == print_rank)
		printf("rank %d entering constructGeom(). \n", rank);

	int max_pairs;
	
	int* pairs_i;
	int* pairs_j;
		
	int* index_to_grid_i;
	int* index_to_grid_j;
	int* index_to_grid_l;
	int* index_to_grid_s;
	
	if (rank == root){
	
		// Build Hstruct object
		std::vector<Sheet> sheets;
		std::vector<int> vec_i;
		std::vector<int> vec_j;
		
		for (int i = 0; i < sdata.size(); ++i){
			sheets.push_back(Sheet(sdata[i]));
		}
		
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
		
		
		// Construct and prepare the pairs array for broadcasting
		// !!!!! MISSING INTRALAYER PAIRS IN h.getPairs() !!!!!
		
		std::vector<std::vector<int> > pairs_vec = h.getPairs();
		max_pairs = static_cast<int>(pairs_vec.size());
		
		MPI::COMM_WORLD.Bcast(&max_pairs, 1, MPI_INT, root);
		
		pairs_i = new int[max_pairs];
		pairs_j = new int[max_pairs];
		
		for(int x = 0; x < max_pairs; ++x){
			pairs_i[x] = pairs_vec[x][0];
			pairs_j[x] = pairs_vec[x][1];
		}
		
	}
	
	if (rank != root){
	
		// Allocate memory to receive pair and "index to grid" information
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(&max_pairs, 1, MPI_INT, root);
		
		pairs_i = new int[max_pairs];
		pairs_j = new int[max_pairs];
		
		index_to_grid_i = new int[max_index];
		index_to_grid_j = new int[max_index];
		index_to_grid_l = new int[max_index];
		index_to_grid_s = new int[max_index];
	}
		
	
	MPI::COMM_WORLD.Bcast(pairs_i, max_pairs, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(pairs_j, max_pairs, MPI_INT, root);
	
	MPI::COMM_WORLD.Bcast(index_to_grid_i, max_index, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(index_to_grid_j, max_index, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(index_to_grid_l, max_index, MPI_INT, root);
	MPI::COMM_WORLD.Bcast(index_to_grid_s, max_index, MPI_INT, root);
	
	// Some C code follows to allocate memory for our completed arrays
	// Perhaps one can get MPI to take std::vector as a valid data type instead?
	
	index_to_grid = (int **) malloc(max_index * sizeof(int *));
	
	for (int k = 0; k < max_index; ++k){
		index_to_grid[k] = (int *) malloc(4 * sizeof(int));
		index_to_grid[k][0] = index_to_grid_i[k];
		index_to_grid[k][1] = index_to_grid_j[k];
		index_to_grid[k][2] = index_to_grid_l[k];
		index_to_grid[k][3] = index_to_grid_s[k];
	}
	
	pairs = (int **) malloc(max_pairs * sizeof(int *));
	
	for (int x = 0; x < max_pairs; ++x){
		pairs[x] = (int *) malloc(2 * sizeof(int));
		pairs[x][0] = pairs_i[x];
		pairs[x][1] = pairs_j[x];
	
	}
			
	if (rank == print_rank){	
		printf("Heterostructure has %d q points. \n", max_index);
		printf("Sparse matrix with %d entries expected. \n", max_pairs);
	}
	

}

void Locality::constructMatrix(){ 
	if (rank == print_rank)
		printf("rank %d entering constructMatrix(). \n", rank);
		
	if (rank == root) {
		rootMatrixSolve();
	} else {
		workerMatrixSolve();
	}

}

void Locality::rootMatrixSolve() {
	
	int nShifts = 5;
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
						&work[currentJob], 	// input buffer
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
					1,
					MPI::DOUBLE,
					MPI::ANY_SOURCE,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
		
		MPI::COMM_WORLD.Send(	
					&work[currentJob],
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
					1,
					MPI::DOUBLE,
					r,
					MPI::ANY_TAG,
					status);				// keeps tag and source information
	
		printf("rank %d has received final work from rank %d. \n", rank, r);
	}
	
	// Tell workers to exit workerMatrixSolve()
	
	for (int r = 1; r < size; ++r){
		double temp[2];
		temp[0] = 0.0;
		temp[1] = 0.0;
		MPI::COMM_WORLD.Send(	
					&temp,
					2, 
					MPI::DOUBLE, 
					r, 
					STOPTAG);
	}
	
}

void Locality::workerMatrixSolve() {

	double result[num_eigs];
	double work[2];
	MPI::Status status;
	
	while (1) {
		MPI::COMM_WORLD.Recv( 
						&work, 
						2, 
						MPI::DOUBLE, 
						root, 
						MPI::ANY_TAG, 
						status);
		
		if (status.Get_tag() == STOPTAG) {
			return;
		}
		
		for (int i = 0; i < num_eigs; ++i){
			result[i] = 12;
		}
		
		MPI::COMM_WORLD.Send(	
					result,
					1,
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
	if (rank == print_rank)
		printf("rank %d finalizing MPI. \n", rank);


	MPI::Finalize();

}

