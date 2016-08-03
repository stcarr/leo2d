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
#include <sstream>
	
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

void Locality::setup(std::string name, int shifts, int eigs, int samples,double start, double end,double e_rescale, double e_shift, double c_width, int p_order, int solver, int intra_search, int inter_search, int magOn_in, double B_in, double vc_in, int nts_in, std::vector<int> ts_in) {
	// Edits run-specific options for matrix constructions and parameters of the solver method (to edit settings from the constructor)

	job_name = name;
	nShifts = shifts;
	
	// Following are no longer used for calculation, done in post-processing
	num_eigs = eigs;
	num_samples = samples;
	interval_start = start;
	interval_end = end;

	// H construction and Chebyshev options
	solver_type = solver;
    energy_rescale = e_rescale;
	energy_shift = e_shift;
	cheb_width = c_width;
	poly_order = p_order;
	intra_searchsize = intra_search;
	inter_searchsize = inter_search;
	
	// Magnetic field options
	magOn = magOn_in;
	B = B_in;
	
	// Vacancy options
	vacancy_chance = vc_in;
	
	// Local DOS targets
	num_target_sheets = nts_in;
	target_sheets = ts_in;
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
		
		center_index = new int[sdata.size()];
		for (int target_sheet = 0; target_sheet < sdata.size(); target_sheet++){
			int target_x_offset = ( sheets[target_sheet].getShape(1,0) - sheets[target_sheet].getShape(0,0) ) / 2;
			int target_y_offset = ( sheets[target_sheet].getShape(1,1) - sheets[target_sheet].getShape(0,1) ) / 2;
			int center_grid[4] = {target_x_offset,target_y_offset,0,target_sheet};
			center_index[target_sheet] = h.gridToIndex(center_grid);
		}
		
		// Get Vacancies if solver_type == 3
		if (solver_type == 3){
			v_work = h.getVacancyList(center_index[0],nShifts*nShifts);
		}
		
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(center_index, sdata.size(), MPI_INT, root);
		
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
	
		printf("Building inter and intra pairs. \n");	
		std::vector<std::vector<int> > inter_pairs_vec;
		h.getInterPairs(inter_pairs_vec,inter_searchsize);
		
		std::vector<int> intra_pairs_vec_i;
		std::vector<int> intra_pairs_vec_j;
		std::vector<double> intra_pairs_vec_t;
		h.getIntraPairs(intra_pairs_vec_i, intra_pairs_vec_j, intra_pairs_vec_t, intra_searchsize);
	
		printf("Inter and intra pair construction complete. \n");
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
		printf("Getting indexToPos arrays. \n");
		index_to_pos_x = new double[max_index];
		index_to_pos_y = new double[max_index];
		index_to_pos_z = new double[max_index];
		
		h.getIndexToPos(index_to_pos_x,0);
		h.getIndexToPos(index_to_pos_y,1);
		h.getIndexToPos(index_to_pos_z,2);
		printf("indexToPos complete. \n");

		// POSITION DEBUG PRINT
		/*
		
		for (int k = 0; k < max_index; ++k)
			printf("%lf, %lf, %lf \n",index_to_pos_x[k],index_to_pos_y[k],index_to_pos_z[k]);
		*/
		//
		
		
		
	}
	
	if (rank != root){
		
		printf("rank %d allocating memory in constructGeom(). \n", rank);
		
		center_index = new int[sdata.size()];
	
		// Allocate memory to receive pair and "index to grid" information
		MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, root);
		MPI::COMM_WORLD.Bcast(center_index, sdata.size(), MPI_INT, root);
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
		index_to_grid[k*4 + 0] = index_to_grid_i[k];
		index_to_grid[k*4 + 1] = index_to_grid_j[k];
		index_to_grid[k*4 + 2] = index_to_grid_l[k];
		index_to_grid[k*4 + 3] = index_to_grid_s[k];
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

	// 1: Chebyshev polynomial sampling of eigenvalue spectrum (DOS)
	if(solver_type == 1 or solver_type == 2 or solver_type == 3){
		if (rank == root) {
			rootChebSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		} else {
			workerChebSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		}	
	}

	time(&solveEnd);
}


void Locality::sendRootWork(int type, int jobIndex, int target_r, std::vector< std::vector<double> > work, std::vector< std::vector<int> > v){

	if (type == 1 or type == 2){
		if(jobIndex == -1){
		
			double temp[2];
			temp[0] = 0.0;
			temp[1] = 0.0;
			MPI::COMM_WORLD.Send(	
					temp,
					2, 
					MPI::DOUBLE, 
					target_r, 
					STOPTAG);
		
		}
		else{
			double work_out[2];
			work_out[0] = work[jobIndex][0];
			work_out[1] = work[jobIndex][1];
	
			MPI::COMM_WORLD.Send(	
						work_out, 	// input buffer
						2,					// size of buffer [x,y]
						MPI::DOUBLE,		// type of buffer
						target_r,			// worker to receive
						jobIndex+1);		// work tag to label shift value
		}		
	} else if (type == 3){
	
		if(jobIndex == -1){
			int temp_length = 2;
			int temp[2];
			temp[0] = 0;
			temp[1] = 0;
			MPI::COMM_WORLD.Send(	
					&temp_length,
					1, 
					MPI::INT, 
					target_r, 
					STOPTAG);
					
			MPI::COMM_WORLD.Send(	
					temp,
					2, 
					MPI::INT, 
					target_r, 
					STOPTAG);
		
		}
		else{
	
			int vac_length = v[jobIndex].size();
			int vac_list[vac_length];
			for (int i = 0; i < vac_length; ++i){
				vac_list[i] = v[jobIndex][i];
			}
			
			MPI::COMM_WORLD.Send(	
						&vac_length, 		// input buffer
						1,					// size of buffer 
						MPI::INT,			// type of buffer
						target_r,			// worker to receive
						jobIndex+1);		// work tag to label work
						
			MPI::COMM_WORLD.Send(	
						vac_list,		 	// input buffer
						vac_length,			// size of buffer
						MPI::INT,			// type of buffer
						target_r,			// worker to receive
						jobIndex+1);		// work tag to label work
		}
	}

}

void Locality::rootChebSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t) {
	
	int maxJobs = nShifts*nShifts;
	int currentJob = 0;
	int length;
	std::vector<std::vector<double> > work (maxJobs, std::vector<double> ( 2, 0 ) );
	std::vector<std::vector<double> > result_array;
	result_array.resize(maxJobs);
	
	// Uniform sample over a grid	
	if (solver_type == 1){
		for (int i = 0; i < nShifts; ++i){
			for (int j = 0; j < nShifts; ++j){
				double x = (1.0/(double) (nShifts))*i;
				double y = (1.0/(double) (nShifts))*j;
				work[i*nShifts + j][0] = x;
				work[i*nShifts + j][1] = y;

			}
		}
	}
	
	// Cut through the unit cell
	if (solver_type == 2){
		for (int i = 0; i < maxJobs; ++i){
			double x = (1.0/((double) maxJobs))*i;
			work[i][0] = x;
			work[i][1] = x;
		}
	}
	
	// No Shifts, used for vacancy in monolayer
	if (solver_type == 3){
		for (int i = 0; i < maxJobs; ++i){
			work[i][0] = 0;
			work[i][1] = 0;
		}
	}
	
	/*	
	// TESTING CODE
	
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
	
	v_work.push_back(std::vector<int> (5,12));
	
	for (int r = 1; r < size; ++r) {
		if (currentJob < maxJobs) {

			sendRootWork(solver_type, currentJob, r, work, v_work);
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
		
		sendRootWork(solver_type, currentJob, status.Get_source(), work, v_work);

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
		
		sendRootWork(solver_type, -1, r, work, v_work);
		
		printf("rank %d sent STOPTAG successfully. \n", rank);
	}
	
	double* polys = new double[poly_order];
	
	for (int i = 0; i < poly_order; ++i){
		polys[i] = i;
	}
	
	
	int num_orbitals[num_target_sheets];
	int total_orbs = 0;
	
	for (int i = 0; i < num_target_sheets; i++){
		num_orbitals[i] = sdata[target_sheets[i]].atom_types.size();
		total_orbs = total_orbs + num_orbitals[i];
	}
	
	int sIndex = 0;
	for (int s_i = 0; s_i < num_target_sheets; s_i++){
		int sheet = target_sheets[s_i];
		for (int orb = 0; orb < num_orbitals[s_i]; ++orb) {

			std::ostringstream sheet_st;
			sheet_st << sheet+1;
			std::string sheet_string = sheet_st.str();
			if ( sheet + 1 < 10){
				sheet_string = "0" + sheet_string;
			}
		
			std::ostringstream orb_st;
			orb_st << orb+1;
			std::string orbital_string = orb_st.str();
			if ( orb + 1 < 10){
				orbital_string = "0" + orbital_string;
			}
			
			std::ofstream outFile;
			const char* extension =".cheb";
			outFile.open ((job_name + "_sheet" + sheet_string + "_orbital" + orbital_string + extension).c_str());
			std::cout << "Saving " << job_name << "_sheet" << sheet_string << "_orbital" << orbital_string << " to file. \n";
			outFile << job_name << " Chebyshev T value outputs \n";
			outFile << "Shift x, Shift y, ... polynomial orders ... \n";
			outFile << "-1, -1";
			
			for(int j = 0; j < poly_order; ++j)
				outFile << ", " << polys[j];
			outFile << "\n";
			
			for(int i = 0; i < maxJobs; ++i){
				outFile << work[i][0] << ", " << work[i][1];
				for(int j = 0; j < poly_order; ++j)
					outFile << ", " << result_array[i][j + poly_order*orb + sIndex];
				outFile << "\n";
			}
			
			outFile.close();
		
		}
		sIndex = sIndex + poly_order*num_orbitals[s_i];
	}
	
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
	
	// Sheet solving information
	int num_orbitals[num_target_sheets];
	int total_orbs = 0;
	
	for (int i = 0; i < num_target_sheets; i++){
		num_orbitals[i] = sdata[target_sheets[i]].atom_types.size();
		total_orbs = total_orbs + num_orbitals[i];
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
		int vac_length;
		int* vac_list;

		//if (rank == print_rank)
			//printf("rank %d waiting for new job... \n",rank);
			
		
		// Recv to get a new job
		
		if (solver_type == 1 || solver_type == 2){
		
			MPI::COMM_WORLD.Recv( 
							work, 			// work determines b-shift
							2, 				// length of work
							MPI::DOUBLE, 	// data-type of work
							root, 			// must come from root
							MPI::ANY_TAG,  	// either WORKTAG or STOPTAG
							status);		// keep MPI status information
							
		} else if (solver_type == 3){
		
			MPI::COMM_WORLD.Recv(	
						&vac_length, 		// input buffer
						1,					// size of buffer 
						MPI::INT,			// type of buffer
						root,				// must come from root
						MPI::ANY_TAG,		// either WORKTAG or STOPTAG
						status);			// keep MPI status information
			
			vac_list = new int[vac_length];
			
			MPI::COMM_WORLD.Recv(	
						vac_list,		 	// input buffer
						vac_length,			// size of buffer
						MPI::INT,			// type of buffer
						root,				// must come from root
						status.Get_tag(),	// make sure to get the right WORKTAG
						status);			// work tag to label work
		}
		
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

		double s1_a[2][2];
		for (int i = 0; i < 2; ++i){
			for (int j = 0; j < 2; ++j){
				s1_a[i][j] = sdata[0].a[i][j];
			}
		}		

		for (int i = 0; i < max_index; ++i) {
		
			if (index_to_grid[i*4 + 3] == 0){
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				}
			if (index_to_grid[i*4 + 3] == 1){
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + work[0]*s1_a[0][0] + work[1]*s1_a[0][1];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + work[0]*s1_a[1][0] + work[1]*s1_a[1][1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				}
				
		}
		
		// -----------------------
		// Build the Sparse matrix
		// -----------------------
		
		printf("rank %d building Sparse Matrix H. \n",rank);
		
		// Indexes how many inter terms we have entered so far
		int inter_counter = 0;
		
		// Indexes how many intra terms we have intered so far
		int intra_counter = 0;
		
		// Total number of expected non-zero matrix elements
		int max_nnz = max_intra_pairs + max_inter_pairs;
		
		// Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
		
		// col_index tells us the col of element i
		// v tells us the value of element i
		// row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
		
		int* col_index = new int[max_nnz];
		
		double* v;
		std::complex<double>* v_c;
		
		if (magOn == 0){
			v = new double[max_nnz];
		}
		else if (magOn == 1){
			v_c = new std::complex<double>[max_nnz];
		}
		
		// Keep track of variables for vacancy simulation
		
		int local_max_index = max_index;
		
		std::vector<int> vacancies;									// indices of the orbitals to remove from the tight-binding model
		std::vector<int> current_index_reduction(max_index+1,0);	// tells us how to relabel our indices in the reduced size matrix
		int num_vacancies = 0;										// total number of vacancies
		
		// Vacancy creation loop
		if (solver_type == 3){
			for (int i = 0; i < max_index; ++i){
			
				current_index_reduction[i] = num_vacancies;

				for (int j = 0; j < vac_length; ++j){
					if(i == vac_list[j]){
					
						vacancies.push_back(i);
						num_vacancies++;
						break;
						
					}
				}
		
			}
			
			current_index_reduction[max_index] = num_vacancies;
			local_max_index = max_index - num_vacancies;
			
			printf("Number of vacancies = %d \n",num_vacancies);
			printf("First few vacancies = [");
			for (int j = 0; j < 6; ++j){
				if (vac_length > j){
					printf("%d, ",vac_list[j]);
				}
			}
			if (vac_length > 5){
				printf("...] \n");
			} else{
				printf("] \n");
			}
		}
		/*
		
		
		if (solver_type == 3){
		
			srand((unsigned)time(NULL));
			
			for (int i = 0; i < max_index; ++i){
			
				current_index_reduction[i] = num_vacancies;
				int protected_orbital = 0;
				
				for (int s = 0; s < num_target_sheets; ++s){
				
					int sheet = target_sheets[s];
					
					for (int o = 0; o < num_orbitals[s]; ++o){
					
						if (i == center_index[sheet] + o){
							protected_orbital = 1;
						}
						
						
					}
					
				}
				
				if (protected_orbital == 1){
					continue;
				}
					
				if (((double)rand()/(double)RAND_MAX) < vacancy_chance){
					vacancies.push_back(i);
					num_vacancies++;
				}
			}
			
			current_index_reduction[max_index] = num_vacancies;
			local_max_index = max_index - num_vacancies;
			
			printf("Number of vacancies = %d \n",num_vacancies);
		
		} // End of vacancy creation loop (i.e. solver_type == 3)
		*/
		
		int* row_pointer = new int[local_max_index+1];

		// Count the current element, i.e. "i = input_counter"
		int input_counter = 0;
	
		// Loop through every orbital (i.e. rows of H)
		for (int k_i = 0; k_i < max_index; ++k_i){
			
			int skip_here1 = 0;
			
			if (solver_type == 3){
				if (current_index_reduction[k_i] + 1 == current_index_reduction[k_i + 1]){
					skip_here1 = 1;
				}
			}
			
			int k = k_i - current_index_reduction[k_i];
					
			// Save starting point of row k
			row_pointer[k] = input_counter;
			
			// While we are still at the correct index in our intra_pairs list:
			bool same_index1 = true;
			while(same_index1) {
			
				int skip_here2 = 0;
			
				// if the first index of intra_pairs changes, we stop
				if (intra_pairs[intra_counter*2 + 0] != k_i) {
					same_index1 = false;
					continue; // go to "while(same_index1)" which will end this loop
				}
				
				int new_k = intra_pairs[intra_counter*2 + 1];

				
				// if we are accounting for defects, we check if the other half of the pair has been removed
				if (solver_type == 3){
					if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
						skip_here2 = 1;
					}
				}
				
				// we save this pair into our sparse matrix format
				if (skip_here1 == 0 && skip_here2 == 0){
					col_index[input_counter] = new_k - current_index_reduction[new_k];
					
					double t;
					
					// if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
					if (new_k == k_i)
						t = (intra_pairs_t[intra_counter] + energy_shift)/energy_rescale;
					// Otherwise we just enter the value just with scaling
					else
						t = intra_pairs_t[intra_counter]/energy_rescale;
					
					//printf("rank %d added intra_pair for index %d: [%d,%d] = %f \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
					
					if (magOn == 0){
						v[input_counter] = t;
					}
					else if (magOn == 1) {
						
						double x1 = i2pos[k*3 + 0];
						double y1 = i2pos[k*3 + 1];
						double x2 = i2pos[new_k*3 + 0];
						double y2 = i2pos[new_k*3 + 1];
						
						double pPhase = peierlsPhase(x1, x2, y1, y2, B);
						v_c[input_counter] = std::polar(t, pPhase);
					}
					++input_counter;
				}
				++intra_counter;
				
			}
		
			// While we are still at the correct index in our inter_pairs list:
			bool same_index2 = true;
			while(same_index2) {
			
				int skip_here2 = 0;
			
				// if the first index of intra_pairs changes, we stop
				if (inter_pairs[inter_counter*2 + 0] != k_i) {
					same_index2 = false;
					continue; // go to "while(same_index2)" which will end this loop
				}
				
				int new_k = inter_pairs[inter_counter*2 + 1];
				
				// if we are accounting for defects, we check if the other half of the pair has been removed
				if (solver_type == 3){
					if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
						skip_here2 = 1;
					}
				}
				
				// we save this pair into our sparse matrix format
				if (skip_here1 == 0 && skip_here2 == 0){
					// get the index of the other orbital in this term
					
					col_index[input_counter] = new_k - current_index_reduction[new_k];
					
					// get the position of both orbitals
					double x1 = i2pos[k_i*3 + 0];
					double y1 = i2pos[k_i*3 + 1];
					double z1 = i2pos[k_i*3 + 2];
					double x2 = i2pos[new_k*3 + 0];
					double y2 = i2pos[new_k*3 + 1];
					double z2 = i2pos[new_k*3 + 2];
					
					// and the orbit tag in their respective unit-cell
					int orbit1 = index_to_grid[k_i*4 + 2];
					int orbit2 = index_to_grid[new_k*4 + 2];
					
					int mat1 = sdata[index_to_grid[k_i*4 + 3]].mat;
					int mat2 = sdata[index_to_grid[new_k*4 + 3]].mat;
					
					// and the angle of the sheet each orbital is on
					double theta1 = angles[index_to_grid[k_i*4 + 3]];
					double theta2 = angles[index_to_grid[new_k*4 + 3]];

					// use all this information to determine coupling energy
					// !!! Currently NOT generalized for materials other than graphene, need to do a material index call for both sheets and pass to a general "inter_coupling" method !!!
					
					double t = interlayer_term(x1, y1, z1, x2, y2, z2, orbit1, orbit2, theta1, theta2, mat1, mat2)/energy_rescale;
					if (t != 0 ){
						if (magOn == 0){
							v[input_counter] = t;
						}
						else if (magOn == 1){
							double pPhase = peierlsPhase(x1, x2, y1, y2, B);
						    v_c[input_counter] = std::polar(t, pPhase);
						}
						//printf("inter [%d,%d] = %lf \n",k,new_k,t);
						++input_counter;
					}
				}
				++inter_counter;

			}
		}
		
		// Save the end point + 1 of the last row
		row_pointer[local_max_index] = input_counter;
		
		// Construct the Sparse Tight-binding Hamiltonian matrix
		//!!	SparseMatrix H(max_index, max_index, v, row_index, col_pointer, max_nnz);
		
		// ------------------------------
		// Following saves Matrix to file
		//
		// Should only be uncommented for 1-job processes, otherwise they will overwrite each other!
		
		/*
		std::ofstream outFile;
		const char* extension = "_matrix.dat";
		outFile.open ((job_name + extension).c_str());
		
		for(int i = 0; i < max_index; ++i){
			int start_index = row_pointer[i];
			int stop_index = row_pointer[i+1];
				for(int j = start_index; j < stop_index; ++j){
					outFile << col_index[j] + 1 << ", " << i + 1 << ", " << v[j] << ", " << i2pos[i*3 + 0] << ", " << i2pos[i*3 + 1] << ", " << i2pos[i*3 + 2] << "\n";
				}
		}
		
		outFile.close();
		*/
		
		// End Matrix Save
		// ---------------
				
		//if(rank == print_rank)
			printf("rank %d starting Chebychev solver work... \n", rank);
		/*
		char mv_type = 'N';
        char matdescra[6] = {'G',' ',' ','C',' ',' '};
        double alpha = 1;
        double beta = 0;
		*/
		
		SpMatrix H = SpMatrix();
		if (magOn == 0){
			H.setup(local_max_index, local_max_index, v,   col_index, row_pointer, max_nnz); 
		}
		else if (magOn == 1){
			H.setup(local_max_index, local_max_index, v_c, col_index, row_pointer, max_nnz);
		}
		
		// Chebyshev values
		
		// Saves T values
		
		double* T_array = new double[poly_order*total_orbs];
		
		if (magOn == 0){

			// Temporary vector for algorithm ("previous" vector T_j-1)
			// Starting vector for Chebyshev method is a unit-vector at the center-orbital
			int sIndex = 0;
			for (int t_s = 0; t_s < num_target_sheets; t_s++) {
			
				int sheet = target_sheets[t_s];
			
				for (int orb = 0; orb < num_orbitals[t_s]; ++orb) {
				
					int target_index = center_index[sheet] + orb - current_index_reduction[center_index[sheet] + orb];
					
					double T_prev[local_max_index];
					
					for (int i = 0; i < local_max_index; ++i){
						T_prev[i] = 0.0;
					}
					
					T_prev[target_index] = 1.0;
				
					// Temporary vector for algorithm ("current" vector T_j)
					double T_j[local_max_index];	
					for (int i = 0; i < local_max_index; ++i){
						T_j[i] = 0.0;
					}		
					
					H.vectorMultiply(T_prev, T_j, 1, 0);
					
					// Temporary vector for algorithm ("next" vector T_j+1)
					double T_next[local_max_index];
					for (int i = 0; i < local_max_index; ++i){
						T_next[i] = 0.0;
					}
					
					// first T value is always 1
					T_array[0 + orb*poly_order + sIndex] = 1;
					
					// Next one is calculated simply
					T_array[1 + orb*poly_order + sIndex] = T_j[target_index];
					
					// Now loop algorithm up to poly_order to find all T values
					// double alpha2 = 2;
					
					// want to do: T_next = 2*H*T_j - T_prev;
					H.vectorMultiply(T_j, T_next, 2, 0);
					for (int c = 0; c < local_max_index; ++c){
						T_next[c] = T_next[c] - T_prev[c];
					}
					
					for (int j = 2; j < poly_order/2; ++j){
					
						// reassign values from previous iteration
						for (int c = 0; c < local_max_index; ++c){
							T_prev[c] = T_j[c];	//T_prev = T_j;
							T_j[c] = T_next[c];	//T_j = T_next;
						}
						
						// get the jth entry
						T_array[j + orb*poly_order + sIndex] = T_j[target_index];
						
						// compute the next vector
						H.vectorMultiply(T_j, T_next, 2, 0);
						for (int c = 0; c < local_max_index; ++c){
							T_next[c] = T_next[c] - T_prev[c];
						}
						
						// use Chebyshev recursion relations to populate the {2j,2j+1} entries.
						if (j >= poly_order/4){
						
							double an_an = 0;
							double anp_an = 0;
							for (int c = 0; c < local_max_index; ++c){
								an_an += T_j[c]*T_j[c];
								anp_an += T_next[c]*T_j[c];
							}
							
							// u_{2n} 	= 2*<a_n|a_n>		- u_0;
							// u_{2n+1} = 2*<a_{n+1}|a_n> 	- u_1; 
							T_array[2*j + orb*poly_order + sIndex] = 2*an_an - T_array[0 + orb*poly_order + sIndex];
							T_array[2*j + 1 + orb*poly_order + sIndex] = 2*anp_an - T_array[1 + orb*poly_order] + sIndex;

						}

						
						// print every 100 steps on print rank
						//if (rank == print_rank)
							//if (j%100 == 0)
								//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
					}
					
				}
				
				sIndex = sIndex + num_orbitals[sheet]*poly_order;
			}
			
		}	// end magOn == 0 block
		else if (magOn == 1) {
		
			// Temporary vector for algorithm ("previous" vector T_j-1)
			// Starting vector for Chebyshev method is a unit-vector at the center-orbital
			int sIndex = 0;
			for (int t_s = 0; t_s < num_target_sheets; t_s++) {
			
				int sheet = target_sheets[t_s];
			
				for (int orb = 0; orb < num_orbitals[t_s]; ++orb) {
				
					int target_index = center_index[sheet] + orb - current_index_reduction[center_index[sheet] + orb];
				
					std::complex<double> T_prev[local_max_index];
					
					for (int i = 0; i < local_max_index; ++i){
						T_prev[i] = 0.0;
					}
					
					T_prev[target_index] = 1.0;
				
					// Temporary vector for algorithm ("current" vector T_j)
					std::complex<double> T_j[local_max_index];	
					for (int i = 0; i < local_max_index; ++i){
						T_j[i] = 0.0;
					}		
					
					H.vectorMultiply(T_prev, T_j, 1, 0);

					// Temporary vector for algorithm ("next" vector T_j+1)
					std::complex<double> T_next[local_max_index];
					for (int i = 0; i < local_max_index; ++i){
						T_next[i] = 0.0;
					}
					
					// first T value is always 1
					T_array[0 + orb*poly_order + sIndex] = 1;
					
					// Next one is calculated simply
					T_array[1 + orb*poly_order + sIndex] = T_j[target_index].real();
					
					// Now loop algorithm up to poly_order to find all T values
					// double alpha2 = 2;
					
					// want to do: T_next = 2*H*T_j - T_prev;
					H.vectorMultiply(T_j, T_next, 2, 0);
					for (int c = 0; c < local_max_index; ++c){
						T_next[c] = T_next[c] - T_prev[c];
					}
					
					for (int j = 2; j < poly_order/2; ++j){
					
						// reassign values from previous iteration
						for (int c = 0; c < local_max_index; ++c){
							T_prev[c] = T_j[c];	//T_prev = T_j;
							T_j[c] = T_next[c];	//T_j = T_next;
						}
						
						// get the jth entry
						T_array[j + orb*poly_order + sIndex] = T_j[target_index].real();
						
						// compute the next vector
						H.vectorMultiply(T_j, T_next, 2, 0);
						for (int c = 0; c < local_max_index; ++c){
							T_next[c] = T_next[c] - T_prev[c];
						}
						
						// use Chebyshev recursion relations to populate the {2j,2j+1} entries.
						if (j >= poly_order/4){
						
							double an_an = 0;
							double anp_an = 0;
							for (int c = 0; c < local_max_index; ++c){
								an_an += (std::conj(T_j[c])*T_j[c]).real();
								anp_an += (std::conj(T_next[c])*T_j[c]).real();
							}
							
							// u_{2n} 	= 2*<a_n|a_n>		- u_0;
							// u_{2n+1} = 2*<a_{n+1}|a_n> 	- u_1; 
							T_array[2*j + orb*poly_order + sIndex] = 2*an_an - T_array[0 + orb*poly_order + sIndex];
							T_array[2*j + 1 + orb*poly_order + sIndex] = 2*anp_an - T_array[1 + orb*poly_order + sIndex];

						}
						
						// print every 100 steps on print rank
						//if (rank == print_rank)
							//if (j%100 == 0)
								//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
					}
					
				}
				
				sIndex = sIndex + num_orbitals[sheet]*poly_order;
			}
			
		} // end magOn == 1 block
		
		// Save time at which solver finished
		time_t tempEnd;
		time(&tempEnd);
		solverTimes.push_back(tempEnd);	
		
		int length = poly_order*total_orbs;
		
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
		
		// End of while(1) means we wait for another instruction from root
		}

}

double Locality::peierlsPhase(double x1, double x2, double y1, double y2, double B_in){
	double preFactor = 3.14e-5; // This should be changed to be eqaul to (2 Pi)/(Flux Quanta), with Flux Quanta = h/e, in units of (T*Ang^2)^-1
	double phase = preFactor*B_in*(x2 - x1)*(0.5)*(y1 + y2);
	return phase;
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
