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

#include <petscksp.h>
#include <slepceps.h>

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

	char help[] = "what this program does in brief can go here.";

	PetscInitialize(&argc,&argv,(char*)0,help);
 

        root = 0;
	print_rank = 1;
	
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
}

void Locality::constructGeom(){
	
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
		
		std::vector<std::vector<int> > inter_pairs_vec = h.getInterPairs();
		
		std::vector<std::vector<int> > intra_pairs_vec_i;
		std::vector<std::vector<int> > intra_pairs_vec_j;
		std::vector<std::vector<double> > intra_pairs_vec_t;
		h.getIntraPairs(intra_pairs_vec_i, intra_pairs_vec_j, intra_pairs_vec_t);
		
		max_inter_pairs = static_cast<int>(inter_pairs_vec.size());
		max_intra_pairs = static_cast<int>(intra_pairs_vec_x.size());
		
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
		
		// Get the index_to_pos array
		
		index_to_pos_x = new double[max_index];
		index_to_pos_y = new double[max_index];
		index_to_pos_z = new double[max_index];
		
		h.getIndexToPos(index_to_pos_x,0);
		h.getIndexToPos(index_to_pos_y,1);
		h.getIndexToPos(index_to_pos_z,2);
		
		
	}
	
	if (rank != root){
	
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
	
	index_to_grid = (int **) malloc(max_index * sizeof(int *));
	
	for (int k = 0; k < max_index; ++k){
		index_to_grid[k] = (int *) malloc(4 * sizeof(int));
		index_to_grid[k][0] = index_to_grid_i[k];
		index_to_grid[k][1] = index_to_grid_j[k];
		index_to_grid[k][2] = index_to_grid_l[k];
		index_to_grid[k][3] = index_to_grid_s[k];
	}
	
	inter_pairs = (int **) malloc(max_inter_pairs * sizeof(int *));
	
	for (int x = 0; x < max_inter_pairs; ++x){
		inter_pairs[x] = (int *) malloc(2 * sizeof(int));
		inter_pairs[x][0] = inter_pairs_i[x];
		inter_pairs[x][1] = inter_pairs_j[x];
	
	}
	
	intra_pairs = (int **) malloc(max_intra_pairs * sizeof(int *));
	intra_pairs_t = (double *) malloc(max_intra_pairs *sizeof(double));
	
	for (int x = 0; x < max_intra_pairs; ++x){
		intra_pairs[x] = (int *) malloc(2 * sizeof(int));
		intra_pairs[x][0] = intra_pairs_i[x];
		intra_pairs[x][1] = intra_pairs_j[x];
		intra_pairs_t[x] = intra_pairs_t[x];
	
	}
	
	index_to_pos = (double **) malloc((max_index) * sizeof(double *));
	
	for (int k = 0; k < max_index; ++k){
		index_to_pos[k] = (double *) malloc(3 * sizeof(double));
		index_to_pos[k][0] = index_to_pos_x[k];
		index_to_pos[k][1] = index_to_pos_y[k];
		index_to_pos[k][2] = index_to_pos_z[k];

	}
			
	if (rank == print_rank){	
		printf("Heterostructure has %d atoms. \n", max_index);
		printf("%d entries expected from intra. \n", max_intra_pairs);
		printf("%d entries expected from inter. \n", max_inter_pairs);
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
					num_eigs,
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
					num_eigs,
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
		
		

		// Build and solve TBH matrix
			
	        printf("Entering PETSc code area. \n");
		
		PetscErrorCode ierr;
		PetscInt N = max_index;
		Mat H;
		PetscInt nnz[N];


		MatCreate(PETSC_COMM_SELF,&H);
		MatSetType(H,MATSEQAIJ);
		MatSetSizes(H,N,N,N,N);
		MatSeqAIJSetPreallocation(H,NULL,nnz);


		// ******

		// some loop here . . .

		PetscInt m = 1; // number of rows being added
		PetscInt idxm = 1; // row index values
	
		for (int k = 0; k < max_index; ++k){
		
			PetscInt n; // number of cols being added
			PetscInt idxn[n]; // col index values
			PetscScalar v[n]; // entry values

			// typically you want m = 1, because v corresponds to the columns
			// from what I understand. So if you want the same row on several rows
			// you can use m > 1, but otherwise just m = 1

			//  ierr = MatSetValues(H, m, &idxm, n, idxn, v, INSERT_VALUES); 
			
			// loop ends

			// *******
		}
		
		MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);
		PetscLogStagePop();
		
		// slepc begins

		
		// slepc ends
		
		printf("5~~ \n");

		MatDestroy(&H);
		
		// delete H matrix
			

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
	if (rank == print_rank)
		printf("rank %d finalizing MPI. \n", rank);


	MPI::Finalize();

}

