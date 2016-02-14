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
 	
	// MPI_Init(&argc,&argv);

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
		
			
		// figure out sparse matrix form (nnz per row)
		
		//nnz = (int *) malloc(max_index * sizeof(int));
		nnz = new int[max_index];		

		int intra_counter = 0;
		int inter_counter = 0;
		
		for (int k = 0; k < max_index; ++k) {
			int nonzeros = 0;
			
			bool same_index1 = true;
			while(same_index1) {
				if (intra_pairs_i[intra_counter] != k) {
					same_index1 = false;
				}
				
				else {
					++nonzeros;
					++intra_counter;
				}
				
			}
			
			bool same_index2 = true;
			while(same_index2) {
				if (inter_pairs_i[inter_counter] != k) {
					same_index2 = false;
				}
				
				else {
					++nonzeros;
					++inter_counter;
				}
				
			}
			
			nnz[k] = nonzeros;
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
		
		nnz = new int[max_index];
		
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
	
	MPI::COMM_WORLD.Bcast(nnz, max_index, MPI_DOUBLE, root);
	
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
	
	PetscErrorCode ierr;
	Mat H,A;
	// int curr_mat = 0;
        EPS 	eps;
	ST 	st;
	KSP	ksp;
	PC 	pc;

	ierr = EPSCreate(PETSC_COMM_SELF,&eps);CHKERRV(ierr);
	ierr = MatCreate(PETSC_COMM_SELF,&H);CHKERRV(ierr);

	// Build and solve TBH matrix
		
	PetscInt N = max_index;
	
	MatSetType(H,MATSEQAIJ);
	MatSetSizes(H,N,N,N,N);
	
	PetscInt petsc_nnz[N];
	
	for (int k = 0; k < N; ++k) {
		//if (rank == print_rank)
		//	printf("rank %d with nnz[%d] = %d. \n",rank,k,nnz[k]);
		petsc_nnz[k] = nnz[k];
	}

	
	MatSeqAIJSetPreallocation(H,NULL,petsc_nnz);


	// ******
	// loop to build our sparse H matrix row-by-row

	// typically you want m = 1, because v corresponds to the columns
	// from what I understand. So if you want the same row on several rows
	// you can use m > 1, but otherwise just m = 1


	PetscInt m = 1; // number of rows being added
	
	int intra_counter = 0;
	int inter_counter = 0;

	printf("rank %d trying to build PETSc Matrix! \n",rank);

	for (int k = 0; k < max_index; ++k){
	
		PetscInt idxm = k;
		int n = nnz[k]; // number of cols being added
		PetscInt idxn[n]; // col index values
		PetscScalar v[n]; // entry values

		// printf("rank %d entering construction loops. \n", rank);			

		int input_counter = 0;
		
		bool same_index1 = true;
		while(same_index1) {
			if (intra_pairs[intra_counter*2 + 0] != k) {
				same_index1 = false;
			}
			else {
				idxn[input_counter] = intra_pairs[intra_counter*2 + 1];
				v[input_counter] = intra_pairs_t[intra_counter];
if (k == idxn[input_counter])
	printf("intra val: (%d,%d,%lf)\n",k,idxn[input_counter],intra_pairs_t[intra_counter]);
				//printf("rank %d added intra_pair for index %d: [%d,%d] \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1]);
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
				idxn[input_counter] = inter_pairs[inter_counter*2 + 1];
				v[input_counter] = -123.0; // NEED TO ADD IN INTERLAYER INTERACTION!!
				//printf("rank %d added inter_pair for index %d: [%d, %d] \n", rank, k, inter_pairs[inter_counter*2 + 0], inter_pairs[inter_counter*2+1]);
				++input_counter;
				++inter_counter;
			}
			
		}
		
		//printf("rank %d attempting to add row to H matrix \n", rank);
		if (m == NULL || &idxm == NULL || n == NULL || idxn == NULL || v == NULL)
			printf("NULL PTR FOUND!! \n \n \n \n ----------------- \n \n \n -------------- \n ");
		if (H == NULL)
			printf("H is null!! \n \n \n");
		ierr = MatSetValues(H, m, &idxm, n, idxn, v, INSERT_VALUES);CHKERRV(ierr); 
		
	}
	
	ierr = MatAssemblyBegin(H,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	ierr = MatAssemblyEnd(H,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
//	ierr = PetscLogStagePop();CHKERRV(ierr);

	while (1) {
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
		
		

		// Build and solve TBH matrix
			
	    	printf("Entering PETSc code area. \n");
	
	
		// PetscErrorCode ierr;
		// Mat H;

		ierr = MatDuplicate(H,MAT_DO_NOT_COPY_VALUES,&A);CHKERRV(ierr);

                PetscBool assembly_check;
                ierr = MatAssembled(H,&assembly_check);CHKERRV(ierr);
		printf("H is assembled: %d\n",assembly_check);
		ierr = MatAssembled(A,&assembly_check);CHKERRV(ierr);
		printf("A is assembled: %d\n",assembly_check);
		
		// loop to build our sparse H matrix row-by-row

		// typically you want m = 1, because v corresponds to the columns
		// from what I understand. So if you want the same row on several rows
		// you can use m > 1, but otherwise just m = 1
	
		m = 1; // number of rows being added
		
		inter_counter = 0;
		intra_counter = 0;
	
		printf("rank %d trying to build PETSc Matrix! \n",rank);
	
		for (int k = 0; k < max_index; ++k){
			
			PetscInt idxm = k;
                	int n = nnz[k]; // number of cols being added
                	PetscInt idxn[n]; // col index values
                	PetscScalar v[n]; // entry values
			int input_counter = 0;
			
			bool same_index1 = true;
			while(same_index1) {
				if (intra_pairs[intra_counter*2 + 0] != k) {
					same_index1 = false;
				}
				else {
					idxn[input_counter] = intra_pairs[intra_counter*2 + 1];
					v[input_counter] = intra_pairs_t[intra_counter];
					// printf("rank %d added intra_pair for index %d: [%d,%d] \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1]);
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
					idxn[input_counter] = inter_pairs[inter_counter*2 + 1];
					v[input_counter] = 69.0; // NEED TO ADD IN INTERLAYER INTERACTION!!
					// printf("rank %d added inter_pair for index %d: [%d, %d] \n", rank, k, inter_pairs[inter_counter*2 + 0], inter_pairs[inter_counter*2+1]);
					++input_counter;
					++inter_counter;
				}
				
			}
			
			// printf("rank %d attempting to add row to H matrix \n", rank);
			if (m == NULL || &idxm == NULL || n == NULL || idxn == NULL || v == NULL)
				printf("NULL PTR FOUND!! \n \n \n \n ----------------- \n \n \n -------------- \n ");
			if (H == NULL)
				printf("H is null :( \n \n \n");
			ierr = MatSetValues(A, m, &idxm, n, idxn, v, INSERT_VALUES);CHKERRV(ierr); 
			
		}
		
		// slepc begins

                ierr = MatAssembled(A,&assembly_check);CHKERRV(ierr);
                printf("A is assembled after insert: %d\n",assembly_check);


        	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
	        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRV(ierr);
		
                ierr = MatAssembled(A,&assembly_check);CHKERRV(ierr);
		printf("A is assembled after Assembly: %d\n", assembly_check);

 		// Old Debug Code
 		//
 		/*
		for (int i = 0; i < max_index; i++)
		{
			PetscScalar value_check;
			PetscInt index_diag = i;
			MatGetValues(A,1,&index_diag,1,&index_diag,&value_check);
			printf("checking (%d,%d): %lf\n", i,i,value_check);
		}
		//
*/

/*		ierr = EPSSetOperators(eps,A,NULL);CHKERRV(ierr);




		PetscReal E1 = -1;
		PetscReal E2 =  1; // energy range of interest

		ierr = EPSSetInterval(eps,E1,E2);CHKERRV(ierr);
	        ierr = EPSSetWhichEigenpairs(eps, EPS_ALL);CHKERRV(ierr);
		// ierr = EPSSetType(eps,EPSKRYLOVSCHUR);CHKERRV(ierr);
        	ierr = EPSGetST(eps,&st);CHKERRV(ierr);
        	ierr = STSetType(st,STSINVERT);CHKERRV(ierr);
        	ierr = STGetKSP(st,&ksp);CHKERRV(ierr);
        	ierr = KSPGetPC(ksp,&pc);CHKERRV(ierr);
        	ierr = KSPSetType(ksp,  KSPPREONLY);CHKERRV(ierr);
        	ierr = PCSetType(pc,PCCHOLESKY);CHKERRV(ierr);
        	ierr = EPSSetFromOptions(eps);CHKERRV(ierr);
        	ierr = EPSSetProblemType(eps,EPS_HEP);CHKERRV(ierr); // sets as hermitian
		PetscInt nconv;

		ierr = EPSSolve(eps);CHKERRV(ierr);

		ierr = EPSGetConverged(eps,&nconv);CHKERRV(ierr); // how many converged eigenpairs

		Vec xr,xi; // real part and imaginary part of eigenvector
		PetscScalar *ki; // real part, imaginary part of eigenvalue
		PetscScalar *kr;
		ki = new PetscScalar[nconv];
		kr = new PetscScalar[nconv];
		PetscInt i; // which eigenpair

		for (int i = 0; i < nconv; i++)
			ierr = EPSGetEigenpair(eps,i,&kr[i],&ki[i],xr,xi);CHKERRV(ierr);

		printf("Beginning eigenvalue printing: \n");
		for (int i = 0; i < nconv; i++)
			printf("%lf + i %lf\n", kr[i],ki[i]);
*/
		// slepc ends
		
		// ierr = MatDestroy(&H);CHKERRV(ierr);
		
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
	ierr = EPSDestroy(&eps);CHKERRV(ierr);
	ierr = MatDestroy(&H);CHKERRV(ierr);
	ierr = MatDestroy(&A);CHKERRV(ierr);

}

void Locality::plot(){ 

}

void Locality::save(){ 

}

void Locality::finMPI(){ 
	printf("rank %d finalizing MPI. \n", rank);
	if (rank == print_rank)
		printf("rank %d finalizing MPI. \n", rank);

	PetscErrorCode ierr;
	ierr = SlepcFinalize();CHKERRV(ierr);

}

