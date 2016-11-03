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
#include <string>
#include <string>
#include <math.h>
	
	
Locality::Locality(std::vector<Sdata> sdata_in,std::vector<double> heights_in,std::vector<double> angles_in) {

	// Set all of the run-specific options for matrix construction and paramters of the solver method used.

	job_name = "HSTRUCT_TEST";

	sdata = sdata_in;
	heights = heights_in;
	angles = angles_in;
	
}

Locality::Locality(const Locality& orig) {

}

Locality::~Locality() {

}

void Locality::setup(Loc_params opts_in){
// Edits run-specific options for matrix constructions and parameters of the solver method (to edit settings from the constructor)
	
	opts = opts_in;
	job_name = opts.getString("job_name");
	
	if (rank == root){
		opts.printParams();
	}
	
	/*
	nShifts = opts.getInt("nShifts");

	// H construction and Chebyshev options
	solver_type = opts.getInt("solver_type");
    energy_rescale = opts.getDouble("energy_rescale");
	energy_shift = opts.getDouble("energy_shift");
	poly_order = opts.getInt("poly_order");
	
	// Magnetic field options
	magOn = opts.getInt("magOn");
	B = opts.getDouble("B");
	
	// Electric field options
	elecOn = opts.getInt("elecOn");
	E = opts.getDouble("E");
	
	// Vacancy options
	vacancy_chance = opts.getDouble("vacancy_chance");
	
	// Local DOS targets
	num_target_sheets = opts.getInt("num_target_sheets");
	target_sheets = opts.getVecInt("target_sheets");
	*/
}

void Locality::getVacanciesFromFile(std::vector<std::vector<int> > &v, std::vector<std::vector<int> > &t){
	
	std::string line;
	std::ifstream in_file;
	in_file.open("vacancies.dat");
	if (in_file.is_open())
	{
		while ( getline(in_file,line) ) {
			
			std::vector<int> temp_v;
			std::vector<int> temp_t;
			
			std::istringstream in_line(line);
		
			std::string in_string;
			
			while ( getline(in_line, in_string, '=') )	{
				
				if (in_string == "JOB"){
				
					std::string v_line;
					std::string t_line;
					
					getline(in_file,v_line);
					getline(in_file,t_line);
					
					std::istringstream in_v_line(v_line);
					std::istringstream in_t_line(t_line);
					
					while ( getline(in_v_line, in_string, ',') )	{
						temp_v.push_back(atoi(in_string.c_str()) - 1);
					}
					while ( getline(in_t_line, in_string, ',') )	{
						temp_t.push_back(atoi(in_string.c_str()) - 1);
					}
					
					v.push_back(temp_v);
					t.push_back(temp_t);
					
				}
				
			}
		}
	}
}

int Locality::initMPI(int argc, char** argv){
	
	// Start  MPI on each rank
	  	
	MPI_Init(&argc,&argv);

    root = 0;
	print_rank = 1;
	
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
	
	if (size == 1)
		return -1;
	else
		return 0;
}

void Locality::constructGeom(){

	// The root node (rank == 0) uses the hstruct and sheet objects to figure out the geometery and tight-binding model sparse matrix structure
	// This information is then pased to all worker nodes

	time(&constructStart);
	
	int solver_type = opts.getInt("solver_type");
	int solver_space = opts.getInt("solver_space");
	int fft_from_file = opts.getInt("fft_from_file");
	int intra_searchsize = opts.getInt("intra_searchsize");
	int inter_searchsize = opts.getInt("inter_searchsize");
	int nShifts = opts.getInt("nShifts");
	
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
	
	// vacancy and target indices
	
	std::vector<std::vector<int> > v_work;
	std::vector<std::vector<int> > target_indices;
		
	
	if (rank == root){
	
		// Build Hstruct object
		std::vector<Sheet> sheets;
		
		for (int i = 0; i < sdata.size(); ++i){
			sheets.push_back(Sheet(sdata[i]));
		}
		
		printf("rank %d building Hstruct. \n", rank);	
		Hstruct h(sheets,angles,heights,solver_space);
		
		// Broadcast "index to grid" mapping information
		max_index = h.getMaxIndex();
		
		// Get Vacancies if solver_type == 3 or 4
		
		if (solver_type == 1 || solver_type == 2){
			target_indices = h.getTargetList(opts);
		}
		
		if (solver_type == 3){
			//v_work = h.getVacancyList(center_index[0],nShifts);
		}
		
		if (solver_type == 4){
			getVacanciesFromFile(v_work, target_indices);
		}
		
		
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
	
		printf("Building inter and intra pairs. \n");	
		std::vector<std::vector<int> > inter_pairs_vec;
		
		// !! MOMENTUM SPACE TO DO !!
		// Probably don't need to change this, just double check!
		h.getInterPairs(inter_pairs_vec,inter_searchsize);
		
		std::vector<int> intra_pairs_vec_i;
		std::vector<int> intra_pairs_vec_j;
		std::vector<double> intra_pairs_vec_t;
		
		h.getIntraPairs(intra_pairs_vec_i, intra_pairs_vec_j, intra_pairs_vec_t, intra_searchsize);
	
		printf("Inter and intra pair construction complete. \n");
		max_inter_pairs = static_cast<int>(inter_pairs_vec.size());
		max_intra_pairs = static_cast<int>(intra_pairs_vec_i.size());

		printf("Heterostructure has %d orbitals. \n", max_index);
		printf("%d entries expected from intra. \n", max_intra_pairs);
		printf("%d entries expected from inter. \n", max_inter_pairs);
		
		if (solver_space == 1){
			if (fft_from_file == 0){
			
				// The following 7 variables should eventually be taken as input parameters in Loc_params.cpp from hstruct.in
					int n_x = 20;
					int n_y = 20;
					int L_x = 20;
					int L_y = 20;
					int length_x = 2;
					int length_y = 2;
					std::string fft_file = "interlayer_fft.dat";
				//
				
				printf("Making *fft.dat file. \n");
				h.makeInterFFTFile(n_x, n_y, L_x, L_y, length_x, length_y, fft_file);
				
				// fft_data is accessed via fft_data[orbital_1][orbital_2][position][real/cpx]
				
				// !! also when putting data[o1][o2] into H, will need to multiply every interlayer term by sqrt(Area(RL_1)*Area(RL_2))
				// !! Using real-space area: This term is ~ 4*pi^2/sqrt(Area1*Area2) ~ 7.4308 for tBLG
		


			} else {
				printf("Warning: Not creating FFT file, assuming a *fft.dat file exists in current directory. \n");
			}
		}
		
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
		printf("Building indexToPos arrays. \n");
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
	int* index_to_grid = new int[max_index*4];
	
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
	
	int* inter_pairs = new int[2*max_inter_pairs];
	
	for (int x = 0; x < max_inter_pairs; ++x){
		inter_pairs[x*2 + 0] = inter_pairs_i[x];
		inter_pairs[x*2 + 1] = inter_pairs_j[x];
	
	}
	
	delete inter_pairs_i;
	delete inter_pairs_j;
	
	int* intra_pairs = new int[2*max_intra_pairs];
	
	for (int x = 0; x < max_intra_pairs; ++x){
		intra_pairs[x*2 + 0] = intra_pairs_i[x];
		intra_pairs[x*2 + 1] = intra_pairs_j[x];
	}
	
	double* index_to_pos = new double[3*max_index];
	
	for (int k = 0; k < max_index; ++k){
		index_to_pos[k*3 + 0] = index_to_pos_x[k];
		index_to_pos[k*3 + 1] = index_to_pos_y[k];
		index_to_pos[k*3 + 2] = index_to_pos_z[k];

	}

	time(&constructEnd);
	
	// Call the next method, which will send out jobs depending on the solver type to each worker and have them construct a matrix specific to that job.
	constructMatrix(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t,v_work,target_indices);
	
	delete index_to_grid;
	delete index_to_pos;
	delete inter_pairs;
	delete intra_pairs;
	delete intra_pairs_t;

}

void Locality::constructMatrix(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > v_work, std::vector< std::vector<int> > target_indices){ 
	time(&solveStart);
	
	int solver_type = opts.getInt("solver_type");
	
	//if (rank == print_rank)
		//printf("rank %d entering constructMatrix(). \n", rank);

	// 1: Chebyshev polynomial sampling of eigenvalue spectrum (DOS)
	if(solver_type < 5){
		if (rank == root) {
			rootChebSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t,v_work,target_indices);
		} else {
			workerChebSolve(index_to_grid,index_to_pos,inter_pairs,intra_pairs,intra_pairs_t);
		}	
	}

	time(&solveEnd);
}

void Locality::sendRootWork(Mpi_job_params jobIn, int target_r){
	
	jobIn.sendParams(target_r,1);
	
}

void Locality::recursiveShiftCalc(std::vector<Mpi_job_params>& jobArray, double* shifts, int solver_type, int nShifts, int maxJobs, int num_sheets, int num_shift_sheets, int* shift_sheets, std::vector< std::vector<int> > target_indices) {

	// Uniform sample over a grid	
	if (solver_type == 1){
	
		if (num_shift_sheets == 0){
			
			int jobID = (int)jobArray.size() + 1;
			
			int n_targets = (int)target_indices[0].size();
			int targets[n_targets];
			
			for (int t = 0; t < n_targets; ++t){
				targets[t] = target_indices[0][t];
			}

			Mpi_job_params tempJob;
			tempJob.loadLocParams(opts);
			tempJob.setParam("shifts",shifts,num_sheets,3);
			tempJob.setParam("jobID",jobID);
			tempJob.setParam("max_jobs",maxJobs);
			tempJob.setParam("target_list",targets,n_targets);
			jobArray.push_back(tempJob);
		
		} else {
		
			int tar_sheet = shift_sheets[num_shift_sheets-1];
		
			for (int i = 0; i < nShifts; ++i){
				for (int j = 0; j < nShifts; ++j){
					
					double x = (1.0/(double) (nShifts))*i;
					double y = (1.0/(double) (nShifts))*j;
					
					shifts[tar_sheet*3 + 0] = x;
					shifts[tar_sheet*3 + 1] = y;
					shifts[tar_sheet*3 + 2] = 0;
					
					recursiveShiftCalc(jobArray, shifts, solver_type, nShifts, maxJobs, num_sheets, num_shift_sheets-1, shift_sheets, target_indices);


				}
			}
		}
		
	}

}

void Locality::rootChebSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > v_work, std::vector< std::vector<int> > target_indices) {
		
	int solver_type = opts.getInt("solver_type");
	int observable_type = opts.getInt("observable_type");
	int solver_space = opts.getInt("solver_space");
	
	int maxJobs;
	int nShifts;
	
	int currentJob = 0;
	int length;
	
	std::vector<Mpi_job_params> jobArray;
	
	int num_sheets = (int)sdata.size();

	if (solver_space == 0) {
	
		// Uniform sample over a grid
		if (solver_type == 1) {
			nShifts = opts.getInt("nShifts");
			int num_shift_sheets = opts.getInt("num_shift_sheets");
			int* shift_sheets = opts.getIntVec("shift_sheets");
			
			maxJobs = pow(nShifts*nShifts,num_shift_sheets);
			
			double shifts[num_sheets*3];
			
			for (int s = 0; s < num_sheets; ++s){
				shifts[s*3 + 0] = 0;
				shifts[s*3 + 1] = 0;
				shifts[s*3 + 2] = 0;
				
			}
			
			recursiveShiftCalc(jobArray,shifts, solver_type, nShifts, maxJobs, num_sheets, num_shift_sheets, shift_sheets, target_indices);
			
		}
		
		// Cut through the unit cell
		if (solver_type == 2){
		
			nShifts = opts.getInt("nShifts");
			maxJobs = nShifts*nShifts;
			
			for (int i = 0; i < maxJobs; ++i){
			
				double x = (1.0/((double) maxJobs))*i;
				
				double shifts[num_sheets*3];
				Mpi_job_params tempJob;
				tempJob = Mpi_job_params();
				
				for(int s = 0; s < num_sheets - 1; ++s){
					shifts[s*3 + 0] = 0;
					shifts[s*3 + 1] = 0;
					shifts[s*3 + 2] = 0;
				}
				
				shifts[(num_sheets-1)*3 + 0] = x;
				shifts[(num_sheets-1)*3 + 1] = x;
				shifts[(num_sheets-1)*3 + 2] = 0;
							
				int n_targets = (int)target_indices[0].size();
				int targets[n_targets];
				
				for (int t = 0; t < n_targets; ++t){
					targets[t] = target_indices[0][t];
				}
			
								
				
				tempJob.loadLocParams(opts);
				tempJob.setParam("shifts",shifts,num_sheets,3);
				tempJob.setParam("jobID",i+1);
				tempJob.setParam("max_jobs",maxJobs);
				tempJob.setParam("target_list",targets,n_targets);
				jobArray.push_back(tempJob);

			}
		}
		
		// Vacancy solvers
		
		if (solver_type == 3 || solver_type == 4){
		
			maxJobs = (int)v_work.size();
			
			for (int i = 0; i < maxJobs; ++i){
				
				Mpi_job_params tempJob;
				
				double shifts[num_sheets*3];

				for(int s = 0; s < num_sheets ; ++s){
					shifts[s*3 + 0] = 0;
					shifts[s*3 + 1] = 0;
					shifts[s*3 + 2] = 0;
				}
					
				int n_vac = (int)v_work[i].size();
				int vacancies[n_vac];
				
				for (int v = 0; v < n_vac; ++v){
					vacancies[v] = v_work[i][v];
				}
				
				int n_targets = (int)target_indices[i].size();
				int targets[n_targets];
				
				for (int t = 0; t < n_targets; ++t){
					targets[t] = target_indices[i][t];
				}
				
							
				tempJob.loadLocParams(opts);
				tempJob.setParam("shifts",shifts,num_sheets,3);
				tempJob.setParam("target_list",targets,n_targets);
				tempJob.setParam("vacancy_list",vacancies,n_vac);
				tempJob.setParam("jobID",i+1);
				tempJob.setParam("max_jobs",maxJobs);
				jobArray.push_back(tempJob);
			}
		}
	} else if (solver_space == 1) {
	// !! MOMENTUM SPACE TO DO !!
	// if Momentum-space, want all sheets to always have the same "shift" (q value)
	// Can still have a LC and grid sampling, and maybe even a "high-symm LC" which does Gamma->M->K->Gamma
		
		// Square sampling
		if (solver_type == 1){
			nShifts = opts.getInt("nShifts");
			maxJobs = nShifts*nShifts;
			
			for (int i = 0; i < nShifts; ++i){
				double x = (1.0/((double) nShifts))*i;
				for (int j = 0; j < nShifts; ++j){
					double y = (1.0/((double) nShifts))*j;
					
					double shifts[num_sheets*3];
					Mpi_job_params tempJob;
					tempJob = Mpi_job_params();
					
					for(int s = 0; s < num_sheets; ++s){
						shifts[s*3 + 0] = x;
						shifts[s*3 + 1] = y;
						shifts[s*3 + 2] = 0;
					}
								
					int n_targets = (int)target_indices[0].size();
					int targets[n_targets];
					
					for (int t = 0; t < n_targets; ++t){
						targets[t] = target_indices[0][t];
					}
				
									
					
					tempJob.loadLocParams(opts);
					tempJob.setParam("shifts",shifts,num_sheets,3);
					tempJob.setParam("jobID",i*nShifts + j + 1);
					tempJob.setParam("max_jobs",maxJobs);
					tempJob.setParam("target_list",targets,n_targets);
					jobArray.push_back(tempJob);
				}

			}
			
		}
		
		// Linecut sampling
		if (solver_type == 2){
			nShifts = opts.getInt("nShifts");
			maxJobs = nShifts*nShifts;
			
			for (int i = 0; i < maxJobs; ++i){
			
				double x = (1.0/((double) maxJobs))*i;
				
				double shifts[num_sheets*3];
				Mpi_job_params tempJob;
				tempJob = Mpi_job_params();
				
				for(int s = 0; s < num_sheets; ++s){
					shifts[s*3 + 0] = x;
					shifts[s*3 + 1] = x;
					shifts[s*3 + 2] = 0;
				}
							
				int n_targets = (int)target_indices[0].size();
				int targets[n_targets];
				
				for (int t = 0; t < n_targets; ++t){
					targets[t] = target_indices[0][t];
				}
			
								
				
				tempJob.loadLocParams(opts);
				tempJob.setParam("shifts",shifts,num_sheets,3);
				tempJob.setParam("jobID",i+1);
				tempJob.setParam("max_jobs",maxJobs);
				tempJob.setParam("target_list",targets,n_targets);
				jobArray.push_back(tempJob);

			}
		}
		
		if (solver_type == 3 || solver_type == 4){
			printf("!!WARNING!!: Momentum-space (solver_space = M) is NOT compatible with vacancy solvers. \n");
		}
	
	}
	
	
	MPI::Status status;
	
	std::vector<std::vector<double> > result_array;
	result_array.resize(maxJobs);
	
	// try to give each worker its first job
	
	for (int r = 1; r < size; ++r) {
		if (currentJob < maxJobs) {
		
			sendRootWork(jobArray[currentJob], r);
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
		
		sendRootWork(jobArray[currentJob], status.Get_source());
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
			
		result_array[jobTag-1] = temp_result;
	
		printf("rank %d has received final work from rank %d. \n", rank, r);
		//if (result_array.size() == maxJobs)
			//break;
	}
	
	// Tell workers to exit workerMatrixSolve()
	
	for (int r = 1; r < size; ++r){
		printf("rank %d sending STOPTAG to rank %d. \n", rank, r);
		
		Mpi_job_params tempJob;
		
		// solver_type == -1 is the STOPTAG for workers
		tempJob.setParam("solver_type", -1);
		
		sendRootWork(tempJob, r);
		
		printf("rank %d sent STOPTAG successfully. \n", rank);
	}
	
	int poly_order = opts.getInt("poly_order");
	double* polys = new double[poly_order];
	for (int i = 0; i < poly_order; ++i){
		polys[i] = i;
	}
	
	int num_target_sheets = opts.getInt("num_target_sheets");
	std::vector<int> target_sheets = opts.getVecInt("target_sheets");
	
	std::cout << "Saving " << job_name << ".cheb to disk. \n";

	
	std::ofstream outFile;
	const char* extension =".cheb";
	outFile.open ((job_name + extension).c_str());
	
	for (int job = 0; job < maxJobs; ++job){
	
		int num_targets = jobArray[job].getInt("num_targets");
		int* target_list = jobArray[job].getIntVec("target_list");

		jobArray[job].printCheb(outFile);
		
		//outFile << job_name << " Chebyshev T value outputs for job " << job+1 << ":\n";
		
		if (observable_type == 0){
			outFile << "T:   ";
			
			for(int j = 0; j < poly_order-1; ++j){
				outFile << polys[j] << ", ";
			}
			outFile << polys[poly_order-1] << "\n";
			
			for(int t = 0; t < num_targets; ++t){
				outFile << target_list[t] << ": ";
				for(int j = 0; j < poly_order-1; ++j){
					outFile << result_array[job][j + t*poly_order] << ", ";
				}
				outFile << result_array[job][poly_order-1 + t*poly_order] << "\n";
			}
			
			outFile << "\n";
			
		} else if (observable_type == 1){
			
			outFile << "T: \n";
			
			for(int t = 0; t < num_targets; ++t){
				outFile << target_list[t] << ": \n";
				for (int r = 0; r < poly_order; ++r){
					for(int j = 0; j < poly_order-1; ++j){
						outFile << result_array[job][j + r*poly_order + t*poly_order*poly_order] << ", ";
					}
					outFile << result_array[job][poly_order-1 + r*poly_order + t*poly_order*poly_order] << "\n";
				}
			}	
			
			outFile << "\n";
		
		}
	}
	
	outFile.close();
	
	delete polys;
}

void Locality::workerChebSolve(int* index_to_grid, double* index_to_pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t) {
	
	MPI::Status status;
	/*
	double work[2];
	work[0] = 0;
	work[1] = 0;
	
	// figure out the Chebyshev polynomial coefficients for each energy sample range
	
	// p is the maximum order of the Chebyshev polynomial
	
	int p = poly_order;
	
	// Sheet solving information
	*/
	
	// ---------------------
	// Enter MPI worker loop
	// ---------------------
	
	while (1) {
	
		Mpi_job_params jobIn;
		
		jobIn.recvParams(root);
		
		int jobID = jobIn.getInt("jobID");
		int max_jobs = jobIn.getInt("max_jobs");
		int solver_type = jobIn.getInt("solver_type");
		int observable_type = jobIn.getInt("observable_type");
		int solver_space = jobIn.getInt("solver_space");
		
		// If worker gets STOPTAG it ends this method
		if (solver_type == -1) {
			//printf("rank %d received STOPTAG. \n", rank);
			return;
		}
		

		// If not STOPTAG, start timing the solver
		time_t tempStart;
		time(&tempStart);
		solverTimes.push_back(tempStart);

		double vacancy_chance = jobIn.getDouble("vacancy_chance");
		
		int num_targets = jobIn.getInt("num_targets");
		int poly_order = jobIn.getInt("poly_order");
		
		double* shifts = jobIn.getDoubleMat("shifts");
		
		printf("rank %d received a job (%d/%d) \n", rank, jobID,max_jobs);
		
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
		
		if (solver_space == 0){
			for (int i = 0; i < max_index; ++i) {
			
				int s = index_to_grid[i*4 + 3];
				
				i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + shifts[s*3 + 0]*s1_a[0][0] + shifts[s*3 + 1]*s1_a[0][1];
				i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + shifts[s*3 + 0]*s1_a[1][0] + shifts[s*3 + 1]*s1_a[1][1];
				i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
			
			}
		} else if (solver_space == 1){
		
			// Here we define Lattice 1 as " A = q + K "
			// and we define Lattice 2 as " B = q - K' "
			// This is done to take care of the the fact that the code looks for K - K' close to zero for inter-pairs (inherited from real space)
			// but in Momentum space the pairs that are "close" should K + K' close to zero!
			// However the monolayer bloch states depend on the sign of K, so we need to make this substitution so the monolayer terms are correct
			
			std::vector< std::vector<double> > b1 = getReciprocal(sdata[0].a);
			
			for (int i = 0; i < max_index; ++i) {
			
				int s = index_to_grid[i*4 + 3];
				if (s == 0){
					i2pos[i*3 + 0] = index_to_pos[i*3 + 0] + shifts[s*3 + 0]*b1[0][0] + shifts[s*3 + 1]*b1[0][1];
					i2pos[i*3 + 1] = index_to_pos[i*3 + 1] + shifts[s*3 + 0]*b1[1][0] + shifts[s*3 + 1]*b1[1][1];
					i2pos[i*3 + 2] = index_to_pos[i*3 + 2];
				} if (s == 1){
					i2pos[i*3 + 0] = -index_to_pos[i*3 + 0] + shifts[s*3 + 0]*b1[0][0] + shifts[s*3 + 1]*b1[0][1];
					i2pos[i*3 + 1] = -index_to_pos[i*3 + 1] + shifts[s*3 + 0]*b1[1][0] + shifts[s*3 + 1]*b1[1][1];
					i2pos[i*3 + 2] = -index_to_pos[i*3 + 2];
				}
			
			}
		
		
		}
		
		// Keep track of variables for vacancy simulation
		
		int local_max_index = max_index;
		
		std::vector<int> current_index_reduction(max_index+1,0);	// tells us how to relabel our indices in the reduced size matrix
		int vacancy_counter = 0;									// total number of vacancies
		
		// Vacancy creation loop
		
		int num_vacancies = jobIn.getInt("num_vacancies");
		int* vac_list = jobIn.getIntVec("vacancy_list");
		
		
		if (solver_type == 3 || solver_type == 4){
			for (int i = 0; i < max_index; ++i){
			
				current_index_reduction[i] = vacancy_counter;

				for (int j = 0; j < num_vacancies; ++j){
					if(i == vac_list[j]){
						vacancy_counter++;
						break;
						
					}
				}
		
			}
			
			current_index_reduction[max_index] = vacancy_counter;
			local_max_index = max_index - vacancy_counter;
			/*
			printf("First few vacancies = [");
			for (int j = 0; j < 6; ++j){
				if (num_vacancies > j){
					printf("%d, ",vac_list[j]);
				}
			}
			if (num_vacancies > 5){
				printf("...] \n");
			} else{
				printf("] \n");
			}
			*/
		}
		
		
		
		// -----------------------
		// Build the Sparse matrix
		// -----------------------
		
		SpMatrix H;
		SpMatrix dxH;
		double* alpha_0_arr = new double[num_targets*local_max_index];
		
		if (solver_space == 0){
			if (observable_type == 0){ // DOS
				generateH(H, jobIn, index_to_grid, i2pos, inter_pairs, intra_pairs, intra_pairs_t, current_index_reduction, local_max_index);
			} else if (observable_type == 1) { // Conductivity
				generateCondH(H, dxH, alpha_0_arr, jobIn, index_to_grid, i2pos, inter_pairs, intra_pairs, intra_pairs_t, current_index_reduction, local_max_index);
			}
		} else if (solver_space == 1){
			generateMomH(H, jobIn, index_to_grid, i2pos, inter_pairs, intra_pairs, intra_pairs_t, current_index_reduction, local_max_index);
		}
		
		// !! MOMENTUM SPACE TO DO !!
		// if solver_space == 1: Instead call generateMomH
		// generateMomH: 
		// 			needs to always be complex
		//			cannot deal with E or B fields yet
		//			should check that there are no vacancies
		//			needs to do bloch-theory on the intralayer blocks
		//			needs to call the interpolate FFT file for interlayer and sample at r2 - r1 + q (i.e. the "shift")
		//
		// Has same input format as generateH(...), because shift is stored in jobIn!
		
		// ---------------------
		// Chebyshev Computation
		// ---------------------
		
		// Saves T values
		double* T_array;
		if (observable_type == 0){
			T_array = new double[poly_order*num_targets];
			computeDosKPM(T_array,H,jobIn,current_index_reduction,local_max_index);
		} else if (observable_type == 1){
			T_array = new double[poly_order*poly_order*num_targets];
			computeCondKPM(T_array, H, dxH, jobIn, current_index_reduction, local_max_index, alpha_0_arr);
		}
		
		// Save time at which solver finished
		time_t tempEnd;
		time(&tempEnd);
		solverTimes.push_back(tempEnd);	
		int length;
		if (observable_type == 0){
			length = poly_order*num_targets;
		} else if (observable_type == 1){
			length = poly_order*poly_order*num_targets;
		}
		
		// Notify root about incoming data size
		MPI::COMM_WORLD.Send(	
					&length,
					1,
					MPI::INT,
					root,
					jobID);

		// Send root the density information
		MPI::COMM_WORLD.Send(	
					T_array,
					length,
					MPI::DOUBLE,
					root,
					jobID);
					
		//if (rank == print_rank)
			//printf("rank %d finished 1 job! \n", rank);
			
		// Cleanup C++ allocated memory
		delete T_array;

		
		// End of while(1) means we wait for another instruction from root
		}

}

void Locality::generateH(SpMatrix &H, Mpi_job_params jobIn, int* index_to_grid, double* i2pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, std::vector<int> current_index_reduction, int local_max_index){

	int solver_type = jobIn.getInt("solver_type");
	int magOn = jobIn.getInt("magOn");
	int elecOn = jobIn.getInt("elecOn");
	double B = jobIn.getDouble("B");
	double E = jobIn.getDouble("E");
	double energy_rescale = jobIn.getDouble("energy_rescale");
	double energy_shift = jobIn.getDouble("energy_shift");
	
	//jobIn.printParams();
	
	
	// Indexes how many inter terms we have entered so far
	int inter_counter = 0;
	
	// Indexes how many intra terms we have entered so far
	int intra_counter = 0;
	
	// Total number of expected non-zero matrix elements
	int max_nnz = max_intra_pairs + max_inter_pairs;
	
	// Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
	
	// col_index tells us the col of element i
	// v tells us the value of element i
	// row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
	
	H.setup(max_nnz, local_max_index, local_max_index);
	
	int* col_index = H.allocColIndx();
	int* row_pointer = H.allocRowPtr();
	
	double* v;
	std::complex<double>* v_c;
	
	if (magOn == 0){
		v = H.allocRealVal();
	}
	else if (magOn == 1){
		v_c = H.allocCpxVal();
	}

	// Count the current element, i.e. "i = input_counter"
	int input_counter = 0;

	// Loop through every orbital (i.e. rows of H)
	for (int k_i = 0; k_i < max_index; ++k_i){
	
		int skip_here1 = 0;
		
		if (solver_type == 3 || solver_type == 4){
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
			if (solver_type == 3 || solver_type == 4){
				if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
					skip_here2 = 1;
				}
			}
			
			// we save this pair into our sparse matrix format
			if (skip_here1 == 0 && skip_here2 == 0){
				col_index[input_counter] = new_k - current_index_reduction[new_k];
				
				double t;
				
				// if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
				// Also, if electric field is included (elecOn == 1) we add in an on-site offset due to this gate voltage.
				if (new_k == k_i){
					if (elecOn == 1){
						t = (intra_pairs_t[intra_counter] + energy_shift + onSiteE(i2pos[k_i*3 + 0],i2pos[k_i*3 + 1],i2pos[k_i*3 + 2],E))/energy_rescale;
					} else if (elecOn == 0)
						t = (intra_pairs_t[intra_counter] + energy_shift)/energy_rescale;
				// Otherwise we enter the value just with rescaling
				}
				else
					t = intra_pairs_t[intra_counter]/energy_rescale;
				
				//printf("rank %d added intra_pair for index %d: [%d,%d] = %f \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
				
				if (magOn == 0){
					v[input_counter] = t;
				}
				else if (magOn == 1) {
					
					double x1 = i2pos[k_i*3 + 0];
					double y1 = i2pos[k_i*3 + 1];
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
			if (solver_type == 3 || solver_type == 4){
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
	
	for(int i = 0; i < local_max_index; ++i){
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
	
}

void Locality::generateCondH(SpMatrix &H, SpMatrix &dxH, double* alpha_0_arr, Mpi_job_params jobIn, int* index_to_grid, double* i2pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, std::vector<int> current_index_reduction, int local_max_index){
	
	int solver_type = jobIn.getInt("solver_type");
	int magOn = jobIn.getInt("magOn");
	int elecOn = jobIn.getInt("elecOn");
	double B = jobIn.getDouble("B");
	double E = jobIn.getDouble("E");
	double energy_rescale = jobIn.getDouble("energy_rescale");
	double energy_shift = jobIn.getDouble("energy_shift");
	
	//jobIn.printParams();
	
	
	// Indexes how many inter terms we have entered so far
	int inter_counter = 0;
	
	// Indexes how many intra terms we have entered so far
	int intra_counter = 0;
	
	// Total number of expected non-zero matrix elements
	int max_nnz = max_intra_pairs + max_inter_pairs;
	
	// Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
	
	// col_index tells us the col of element i
	// v tells us the value of element i
	// row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
	
	H.setup(max_nnz, local_max_index, local_max_index);
	dxH.setup(max_nnz, local_max_index, local_max_index);

	
	int* col_index = H.allocColIndx();
	int* row_pointer = H.allocRowPtr();
	
	int* col_index_dx = dxH.allocColIndx();
	int* row_pointer_dx = dxH.allocRowPtr();

	
	double* v;
	double* v_dx;
	std::complex<double>* v_c;
	std::complex<double>* v_c_dx;

	
	if (magOn == 0){
		v = H.allocRealVal();
		v_dx = dxH.allocRealVal();
	}
	else if (magOn == 1){
		v_c = H.allocCpxVal();
		v_c_dx = dxH.allocCpxVal();
	}

	// Count the current element, i.e. "i = input_counter"
	int input_counter = 0;

	// Loop through every orbital (i.e. rows of H)
	for (int k_i = 0; k_i < max_index; ++k_i){
	
		int skip_here1 = 0;
		
		if (solver_type == 3 || solver_type == 4){
			if (current_index_reduction[k_i] + 1 == current_index_reduction[k_i + 1]){
				skip_here1 = 1;
			}
		}
		
		int k = k_i - current_index_reduction[k_i];
		
		// Save starting point of row k
		row_pointer[k] = input_counter;
		row_pointer_dx[k] = input_counter;
		
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
			if (solver_type == 3 || solver_type == 4){
				if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
					skip_here2 = 1;
				}
			}
			
			// we save this pair into our sparse matrix format	
			if (skip_here1 == 0 && skip_here2 == 0){
			
				col_index[input_counter] = new_k - current_index_reduction[new_k];
				col_index_dx[input_counter] =  new_k - current_index_reduction[new_k];
				
				double t;
				
				// if it is the diagonal element, we "shift" the matrix up or down in energy scale (to make sure the spectrum fits in [-1,1] for the Chebyshev method)
				// Also, if electric field is included (elecOn == 1) we add in an on-site offset due to this gate voltage.
				if (new_k == k_i){
					if (elecOn == 1){
						t = (intra_pairs_t[intra_counter] + energy_shift + onSiteE(i2pos[k_i*3 + 0],i2pos[k_i*3 + 1],i2pos[k_i*3 + 2],E))/energy_rescale;
					} else if (elecOn == 0)
						t = (intra_pairs_t[intra_counter] + energy_shift)/energy_rescale;
				// Otherwise we enter the value just with rescaling
				}
				else
					t = intra_pairs_t[intra_counter]/energy_rescale;
				
				//printf("rank %d added intra_pair for index %d: [%d,%d] = %f \n", rank, k, intra_pairs[intra_counter*2 + 0], intra_pairs[intra_counter*2 + 1],v[input_counter]);
				
				if (magOn == 0){
					v[input_counter] = t;
					v_dx[input_counter] = (i2pos[k_i*3 + 0] - i2pos[new_k*3 + 0])*t;
				}
				else if (magOn == 1) {
					
					double x1 = i2pos[k_i*3 + 0];
					double y1 = i2pos[k_i*3 + 1];
					double x2 = i2pos[new_k*3 + 0];
					double y2 = i2pos[new_k*3 + 1];
					
					double pPhase = peierlsPhase(x1, x2, y1, y2, B);
					v_c[input_counter] = std::polar(t, pPhase);
					v_c_dx[input_counter] = std::polar((i2pos[k_i*3 + 0] - i2pos[new_k*3 + 0])*t,pPhase);
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
			if (solver_type == 3 || solver_type == 4){
				if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
					skip_here2 = 1;
				}
			}
			
			// we save this pair into our sparse matrix format
			if (skip_here1 == 0 && skip_here2 == 0){
				// get the index of the other orbital in this term
				
				col_index[input_counter] = new_k - current_index_reduction[new_k];
				col_index_dx[input_counter] = new_k - current_index_reduction[new_k];
				
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
						v_dx[input_counter] = (x1 - x2)*t;
					}
					else if (magOn == 1){
						double pPhase = peierlsPhase(x1, x2, y1, y2, B);
						v_c[input_counter] = std::polar(t, pPhase);
						v_c_dx[input_counter] = std::polar((x1 - x2)*t, pPhase);
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
	row_pointer_dx[local_max_index] = input_counter;
	
	int num_targets = jobIn.getInt("num_targets");
	int* target_list = jobIn.getIntVec("target_list");
	
	for (int t = 0; t < num_targets; ++t){
		
		int target_here = target_list[t];
		double target_vec[local_max_index];
		for (int i = 0; i < local_max_index; ++i){
			target_vec[i] = 0;
		}
		target_vec[target_here] = 1;
		
		double temp_vec[local_max_index];
		for (int i = 0; i < local_max_index; ++i){
			temp_vec[i] = 0;
		}
		
		dxH.vectorMultiply(target_vec,temp_vec,1,0);
		
		// want < 0 | dxH, not dxH | 0 >, so need to use: Transpose( < 0 | dxH ) = - dxH | 0 >
		for (int i = 0; i < local_max_index; ++ i){
			alpha_0_arr[t*local_max_index + i] = -temp_vec[i];
		}
	}
	
	/*
			// ------------------------------
			// Following saves Matrix to file
			//
			// Should only be uncommented for 1-job processes, otherwise they will overwrite each other!
			
			//*
			std::ofstream outFile;
			const char* extension = "_matrix.dat";
			outFile.open ((job_name + extension).c_str());
			
			for(int i = 0; i < local_max_index; ++i){
				int start_index = row_pointer[i];
				int stop_index = row_pointer[i+1];
					for(int j = start_index; j < stop_index; ++j){
						outFile << col_index[j] + 1 << ", " << i + 1 << ", " << v[j] << ", " << i2pos[i*3 + 0] << ", " << i2pos[i*3 + 1] << ", " << i2pos[i*3 + 2] << "\n";
					}
			}
			
			outFile.close();
			//
			
			// End Matrix Save
			// ---------------
			
			
			// ------------------------------
			// Following saves Matrix to file
			//
			// Should only be uncommented for 1-job processes, otherwise they will overwrite each other!
			
			///*
			std::ofstream outFile2;
			const char* extension2 = "_dxH_matrix.dat";
			outFile2.open ((job_name + extension2).c_str());
			
			for(int i = 0; i < local_max_index; ++i){
				int start_index = row_pointer_dx[i];
				int stop_index = row_pointer_dx[i+1];
					for(int j = start_index; j < stop_index; ++j){
						outFile2 << col_index_dx[j] + 1 << ", " << i + 1 << ", " << v_dx[j] << ", " << i2pos[i*3 + 0] << ", " << i2pos[i*3 + 1] << ", " << i2pos[i*3 + 2] << "\n";
					}
			}
			
			outFile2.close();
			//
			
			// End Matrix Save
			// ---------------
			
	*/
	
}

void Locality::generateMomH(SpMatrix &H, Mpi_job_params jobIn, int* index_to_grid, double* i2pos, int* inter_pairs, int* intra_pairs, double* intra_pairs_t, std::vector<int> current_index_reduction, int local_max_index){

	int solver_type = jobIn.getInt("solver_type");
	int magOn = jobIn.getInt("magOn");
	int elecOn = jobIn.getInt("elecOn");
	double B = jobIn.getDouble("B");
	double E = jobIn.getDouble("E");
	double energy_rescale = jobIn.getDouble("energy_rescale");
	double energy_shift = jobIn.getDouble("energy_shift");
	
	double* shifts = jobIn.getDoubleMat("shifts");
	
	double shift_x = shifts[0];
	double shift_y = shifts[1];
	
	// Indexes how many inter terms we have entered so far
	int inter_counter = 0;
	
	// Indexes how many intra terms we have entered so far
	int intra_counter = 0;
	
	// Total number of expected non-zero matrix elements
	int max_nnz = max_intra_pairs + max_inter_pairs;
	
	// Do monolayer bloch theory for each sheet:
	
	int num_sheets = jobIn.getInt("num_sheets");
	
	int L = 5; // L is the half-side length of a 2L+1 x 2L+1 search grid for the monolayer bloch terms
	
	std::vector<int> N_R_array; // num_sheets vector of the number of real-space positions in each sheet
	std::vector< std::vector< std::vector<double> > > R_array; // (num_sheets) x (N_R) x (3) array of real-space positions
	std::vector< std::vector< std::vector< std::vector<double> > > > bloch_t_array; // (num_sheets) x (N_R) x (num_orbitals) x (num_orbitals) of real-space t_ij coupling values
	
	for (int s = 0; s < num_sheets; ++s){
	
		int N_R = 0;
		
		std::vector< std::vector<double> > temp_R_array;
		std::vector< std::vector< std::vector<double> > > temp_bloch_t_array;
		
		std::vector<std::vector<double> > a = sdata[s].a;
		double theta = angles[s];
		int num_orbs = (int) sdata[s].atom_types.size();
		int mat = sdata[s].mat;
		
		std::vector<std::vector<double> > orb_pos = sdata[s].atom_pos;
		
		for (int i = -L; i < L+1; ++i){
			for (int j = -L; j < L+1; ++j){
			
			
				// unrotated variables
				
				double temp_x = i*a[0][0] + j*a[0][1];
				double temp_y = i*a[1][0] + j*a[1][1];
				double temp_z = 0;
				
				// rotated variables	
				double x = temp_x*cos(theta) - temp_y*sin(theta);
				double y = temp_y*cos(theta) + temp_x*sin(theta);
				double z = temp_z;
				
				
				std::vector< std::vector<double> > temp_o1_bloch_t_array;
				for (int o1 = 0; o1 < num_orbs; ++o1){
					std::vector<double> temp_o2_bloch_t_array;
					for (int o2 = 0; o2 < num_orbs; ++o2){
					
						// compute R(i,j) based on <a> and <theta>
						// find monolayer t_ij(R(i,j),o1,o2)
						// if t_ij != 0, add R and t_ij to respective vectors, and ++N_R;
						
						// !! Include orbital displacement in x,y,z*E
						
						double o1_x = orb_pos[o1][0];
						double o1_y = orb_pos[o1][1];
						double o1_z = orb_pos[o1][2];
						
						double o2_x = temp_x + orb_pos[o2][0];
						double o2_y = temp_y + orb_pos[o2][1];
						double o2_z = temp_z + orb_pos[o2][2];
						
						double t = intralayer_term(o1_x, o1_y, o1_z, o2_x, o2_y, o2_z, o1, o2, mat);
						temp_o2_bloch_t_array.push_back(t);
					}
					temp_o1_bloch_t_array.push_back(temp_o2_bloch_t_array);
				}
				
				std::vector<double> temp_r_vec;
				temp_r_vec.push_back(x);
				temp_r_vec.push_back(y);
				temp_r_vec.push_back(z);
				
				temp_R_array.push_back(temp_r_vec);
				temp_bloch_t_array.push_back(temp_o1_bloch_t_array);
			}
		}
		
		R_array.push_back(temp_R_array);
		bloch_t_array.push_back(temp_bloch_t_array);
		N_R_array.push_back( (int) temp_R_array.size());
	}
	
	
	
	// Sparse matrix format is 2 arrays with length = nnz, and 1 array with length = max_index + 1
	
	// col_index tells us the col of element i
	// v tells us the value of element i
	// row_pointer tells us the start of row j (at element i = row_pointer[j]) and the end of row j (at element i = row_pointer[j] - 1)
	
	H.setup(max_nnz, local_max_index, local_max_index);
	
	int* col_index = H.allocColIndx();
	int* row_pointer = H.allocRowPtr();
	
	std::complex<double>* v_c;
	
	v_c = H.allocCpxVal();

	// Count the current element, i.e. "i = input_counter"
	int input_counter = 0;

	// Load FFTW data into Interlayer_coupling object
	
	// The following 7 variables should eventually be taken as input parameters in Loc_params.cpp from hstruct.in
		int n_x = 20;
		int n_y = 20;
		int L_x = 20;
		int L_y = 20;
		int length_x = 2;
		int length_y = 2;
		std::string fft_file = "interlayer_fft.dat";
	//
	
	Interlayer_coupling fftw_inter;
	fftw_inter.fft_setup(L_x,L_y,length_x,length_y,fft_file);
	
	//printf("Momentum Space coupling term [0][0] at k = (0,0) = %lf + i*%lf \n",7.4308*fftw_inter.interp_fft(0.0,0.0,0,0,0),7.4308*fftw_inter.interp_fft(0.0,0.0,0,0,1));
	
	
	// fft_data is accessed via fft_data[orbital_1][orbital_2][position][real/cpx]
	
	// !! also when putting data[o1][o2] into H, will need to multiply every interlayer term by sqrt(Area(RL_1)*Area(RL_2))
	// !! Using real-space area: This term is ~ 4*pi^2/sqrt(Area1*Area2) ~ 7.4308 for tBLG
		
	
	// We keep the vacancy methods here, although momentum space in general has no vacancies...
	
	// Loop through every orbital (i.e. rows of H)
	for (int k_i = 0; k_i < max_index; ++k_i){
	
		int skip_here1 = 0;
		
		if (solver_type == 3 || solver_type == 4){
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
			if (solver_type == 3 || solver_type == 4){
				if(current_index_reduction[new_k] + 1 == current_index_reduction[new_k + 1]){
					skip_here2 = 1;
				}
			}
			
			// we save this pair into our sparse matrix format
			if (skip_here1 == 0 && skip_here2 == 0){
				col_index[input_counter] = new_k - current_index_reduction[new_k];
				
				std::complex<double> t(0.0,0.0);
				
				// Need to do monolayer bloch theory here!
				
				// finding pairing between k_i and new_k
				

				// and the orbit tag in their respective unit-cell
				int orbit1 = index_to_grid[k_i*4 + 2];
				int orbit2 = index_to_grid[new_k*4 + 2];
				int s = index_to_grid[k_i*4 + 3];
				
				int N_R = N_R_array[s];
				
				for (int n = 0; n < N_R; ++n){
					
					double q1 = i2pos[k_i*3 + 0];
					double q2 = i2pos[k_i*3 + 1];
					double q3 = i2pos[k_i*3 + 2];
					
					double r1 = R_array[s][n][0];
					double r2 = R_array[s][n][1];
					double r3 = R_array[s][n][2];
					
					double phase = q1*r1 + q2*r2 + q3*r3;
					
					double t_here = bloch_t_array[s][n][orbit1][orbit2];
					
					
					t = t + std::polar(t_here,phase);
				}
				
				
				v_c[input_counter] = t/energy_rescale;
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
			if (solver_type == 3 || solver_type == 4){
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
				
				
				// Momentum-space Interlayer coupling is evaluated at k = k2 + k1 - q
				// This is because we defined k1 = q + K orig, and k2 = q - K' original, and THAT was done because we did not want to change the code for inter-pairs in hstruct
				
				double dx = x2 + x1 - shift_x;
				double dy = y2 + y1 - shift_y;
				
				std::complex<double> t;
				
				// FFT file is based on sheet1 -> sheet2, so we need to keep track of which sheet k_i and new_k are on (i.e. k_i < new_k -> k_i on 1st sheet, k_i > new_k -> k_i on 2nd sheet)
				// Not checking for this is the equivalent of twisting the layer "theta" for 1->2 coupling and "-theta" for 2->1 coupling, which makes H non-Hermitian.
				
				
				if (k_i <= new_k){
					t = std::complex<double>(7.4308*fftw_inter.interp_fft(dx,dy,orbit1,orbit2,0),7.4308*fftw_inter.interp_fft(dx,dy,orbit1,orbit2,1));
				} else if (k_i > new_k) {
					t = std::complex<double>(7.4308*fftw_inter.interp_fft(dx,dy,orbit2,orbit1,0),7.4308*fftw_inter.interp_fft(dx,dy,orbit2,orbit1,1));
				}
				
				// double t = interlayer_term(x1, y1, z1, x2, y2, z2, orbit1, orbit2, theta1, theta2, mat1, mat2)/energy_rescale;
				if (std::abs(t) != 0){
					
					if (k_i <= new_k){
						v_c[input_counter] = t/energy_rescale;
					} else if (k_i > new_k) {
						v_c[input_counter] = std::conj(t)/energy_rescale;
					}
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
	const char* extension = "_matrix_real.dat";
	outFile.open ((job_name + extension).c_str());
	
	for(int i = 0; i < local_max_index; ++i){
		int start_index = row_pointer[i];
		int stop_index = row_pointer[i+1];
			for(int j = start_index; j < stop_index; ++j){
				outFile << col_index[j] + 1 << ", " << i + 1 << ", " << std::real(v_c[j]) << ", " << i2pos[i*3 + 0] << ", " << i2pos[i*3 + 1] << ", " << i2pos[i*3 + 2] << "\n";
			}
	}
	
	outFile.close();
	
	std::ofstream outFile2;
	const char* extension2 = "_matrix_cpx.dat";
	outFile2.open ((job_name + extension2).c_str());
	
	for(int i = 0; i < local_max_index; ++i){
		int start_index = row_pointer[i];
		int stop_index = row_pointer[i+1];
			for(int j = start_index; j < stop_index; ++j){
				outFile2 << col_index[j] + 1 << ", " << i + 1 << ", " << std::imag(v_c[j]) << ", " << i2pos[i*3 + 0] << ", " << i2pos[i*3 + 1] << ", " << i2pos[i*3 + 2] << "\n";
			}
	}
	
	outFile.close();
	*/
	
	// End Matrix Save
	// ---------------

	
}

void Locality::computeDosKPM(double* T_array, SpMatrix &H, Mpi_job_params jobIn, std::vector<int> current_index_reduction, int local_max_index){

	int magOn = jobIn.getInt("magOn");
	int poly_order = jobIn.getInt("poly_order");
	int solver_space = jobIn.getInt("solver_space");
	
	int num_targets = jobIn.getInt("num_targets");
	int* target_list = jobIn.getIntVec("target_list");

	
	if (magOn == 0 && solver_space == 0){

		// Starting vector for Chebyshev method is a unit-vector at the target orbital
		
		for (int t_count = 0; t_count < num_targets; ++t_count){
	
			int target_index = target_list[t_count] - current_index_reduction[target_list[t_count]];
			
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
			T_array[0 + t_count*poly_order] = 1;
			
			// Next one is calculated simply
			T_array[1 + t_count*poly_order] = T_j[target_index];
			
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
				T_array[j + t_count*poly_order] = T_j[target_index];
				
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
					T_array[2*j + t_count*poly_order ] = 2*an_an - T_array[0 + t_count*poly_order];
					T_array[2*j + 1 + t_count*poly_order] = 2*anp_an - T_array[1 + t_count*poly_order ];

				}

				
				// print every 100 steps on print rank
				//if (rank == print_rank)
					//if (j%100 == 0)
						//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
			}
		
		}

	}	// end magOn == 0 block
	else if (magOn == 1 || solver_space == 1) {
	
		// Starting vector for Chebyshev method is a unit-vector at the target rbital
		for (int t_count = 0; t_count < num_targets; ++t_count){
	
			int target_index = target_list[t_count] - current_index_reduction[target_list[t_count]];			
			
			std::complex<double>* T_prev = new std::complex<double>[local_max_index];
			
			for (int i = 0; i < local_max_index; ++i){
				T_prev[i] = 0.0;
			}
			
			T_prev[target_index] = 1.0;
		
			// Temporary vector for algorithm ("current" vector T_j)
			std::complex<double>* T_j = new std::complex<double>[local_max_index];	
			for (int i = 0; i < local_max_index; ++i){
				T_j[i] = 0.0;
			}		
			
			H.vectorMultiply(T_prev, T_j, 1, 0);

			// Temporary vector for algorithm ("next" vector T_j+1)
			std::complex<double>* T_next = new std::complex<double>[local_max_index];
			for (int i = 0; i < local_max_index; ++i){
				T_next[i] = 0.0;
			}
			
			// first T value is always 1
			T_array[0 + t_count*poly_order] = 1;
			
			// Next one is calculated simply
			T_array[1 + t_count*poly_order] = T_j[target_index].real();
			
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
				T_array[j + t_count*poly_order] = T_j[target_index].real();
				
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
					T_array[2*j + t_count*poly_order] = 2*an_an - T_array[0 + t_count*poly_order];
					T_array[2*j + 1 + t_count*poly_order] = 2*anp_an - T_array[1 + t_count*poly_order];

				}
				
				// print every 100 steps on print rank
				//if (rank == print_rank)
					//if (j%100 == 0)
						//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
			}
			
			delete T_prev;
			delete T_j;
			delete T_next;
			
		}
		
	} // end magOn == 1 block
}

void Locality::computeCondKPM(double* T_array, SpMatrix &H, SpMatrix &dxH, Mpi_job_params jobIn, std::vector<int> current_index_reduction, int local_max_index, double* alpha_0){

	
	int magOn = jobIn.getInt("magOn");
	int poly_order = jobIn.getInt("poly_order");
	
	int num_targets = jobIn.getInt("num_targets");
	int* target_list = jobIn.getIntVec("target_list");

	if (magOn == 0){

		// Starting vector for Chebyshev method is a unit-vector at the target orbital
		for (int t_count = 0; t_count < num_targets; ++t_count){
			
			double* beta_array;
			beta_array = new double[local_max_index*poly_order];
			
			int target_index = target_list[t_count] - current_index_reduction[target_list[t_count]];
			
			double T_prev[local_max_index];
			
			for (int i = 0; i < local_max_index; ++i){
				T_prev[i] = 0.0;
			}
			
			T_prev[target_index] = 1.0;
			
			for (int i = 0; i < local_max_index; ++i){
				beta_array[i] = T_prev[i];
			}
			
			// Temporary vector for algorithm ("current" vector T_j)
			double T_j[local_max_index];	
			for (int i = 0; i < local_max_index; ++i){
				T_j[i] = 0.0;
			}		
			
			H.vectorMultiply(T_prev, T_j, 1, 0);
			
			for (int i = 0; i < local_max_index; ++i){
				beta_array[local_max_index + i] = T_j[i];
			}
			
			// Temporary vector for algorithm ("next" vector T_j+1)
			double T_next[local_max_index];
			for (int i = 0; i < local_max_index; ++i){
				T_next[i] = 0.0;
			}
			
			// Now loop algorithm along T_m*|beta>, up to poly_order
			// double alpha2 = 2;
			
			// want to do: T_next = 2*H*T_j - T_prev;
			H.vectorMultiply(T_j, T_next, 2, 0);
			for (int c = 0; c < local_max_index; ++c){
				T_next[c] = T_next[c] - T_prev[c];
			}
			
			for (int j = 2; j < poly_order; ++j){
			
				// reassign values from previous iteration
				for (int c = 0; c < local_max_index; ++c){
					T_prev[c] = T_j[c];	//T_prev = T_j;
					T_j[c] = T_next[c];	//T_j = T_next;
				}
				
				for (int i = 0; i < local_max_index; ++i){
					beta_array[j*local_max_index + i] = T_j[i];
				}
				
				// compute the next vector
				H.vectorMultiply(T_j, T_next, 2, 0);
				for (int c = 0; c < local_max_index; ++c){
					T_next[c] = T_next[c] - T_prev[c];
				}
				
				// print every 100 steps on print rank
				//if (rank == print_rank)
					//if (j%100 == 0)
						//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
			}	// end beta_array construction
			
			// start loop over <alpha|*T_n for C_nm construction
			
			for (int i = 0; i < local_max_index; ++i){
				T_prev[i] = alpha_0[t_count*local_max_index + i];
			}
			
			for (int p = 0; p < poly_order; ++p){
			
				// temp_vec is = - <alpha|dxH ... because dxH.vectorMultiply defines dxH|alpha> and Transpose(dxH) = -dxH
				double temp_vec[local_max_index];
				dxH.vectorMultiply(T_prev, temp_vec,1,0);

				double temp_sum = 0;
				for (int i = 0; i < local_max_index; ++i){
					temp_sum = temp_sum + temp_vec[i]*beta_array[p*local_max_index + i];
				}

				T_array[p + 0*poly_order + t_count*poly_order*poly_order] = temp_sum;
				
			}
			
			// Temporary vector for algorithm ("current" vector T_j)
			for (int i = 0; i < local_max_index; ++i){
				T_j[i] = 0.0;
			}	

			// This is OK because H is Hermitian, Transpose(H) = H, so <alpha|H can be computed with H|alpha>
			H.vectorMultiply(T_prev, T_j, 1, 0);
			for (int p = 0; p < poly_order; ++p){
			
				// temp_vec is = - <alpha|dxH ... because dxH.vectorMultiply defines dxH|alpha> and Transpose(dxH) = -dxH
				double temp_vec[local_max_index];
				dxH.vectorMultiply(T_j, temp_vec,1,0);

				double temp_sum = 0;
				for (int i = 0; i < local_max_index; ++i){
					temp_sum = temp_sum + temp_vec[i]*beta_array[p*local_max_index + i];
				}

				T_array[p + 1*poly_order + t_count*poly_order*poly_order] = temp_sum;
				
			}
			
			// Temporary vector for algorithm ("next" vector T_j+1)
			for (int i = 0; i < local_max_index; ++i){
				T_next[i] = 0.0;
			}
		
			// Now loop algorithm up to poly_order to find all T values
			// double alpha2 = 2;
			
			// want to do: T_next = 2*H*T_j - T_prev;
			H.vectorMultiply(T_j, T_next, 2, 0);
			for (int c = 0; c < local_max_index; ++c){
				T_next[c] = T_next[c] - T_prev[c];
			}
			
			for (int j = 2; j < poly_order; ++j){
			
				// reassign values from previous iteration
				for (int c = 0; c < local_max_index; ++c){
					T_prev[c] = T_j[c];	//T_prev = T_j;
					T_j[c] = T_next[c];	//T_j = T_next;
				}
				
				
				for (int p = 0; p < poly_order; ++p){
				
					// temp_vec is = - <alpha|dxH ... because dxH.vectorMultiply defines dxH|alpha> and Transpose(dxH) = -dxH
					double temp_vec[local_max_index];
					dxH.vectorMultiply(T_j, temp_vec,1,0);
					
					double temp_sum = 0;
					for (int i = 0; i < local_max_index; ++i){
						temp_sum = temp_sum + temp_vec[i]*beta_array[p*local_max_index + i];
					}
					T_array[p + j*poly_order + t_count*poly_order*poly_order] = temp_sum;
					
				}
				
				// compute the next vector
				H.vectorMultiply(T_j, T_next, 2, 0);
				for (int c = 0; c < local_max_index; ++c){
					T_next[c] = T_next[c] - T_prev[c];
				}

				
				// print every 100 steps on print rank
				//if (rank == print_rank)
					//if (j%100 == 0)
						//printf("Chebyshev iteration (%d/%d) complete. \n",j,poly_order);
			}
			
			delete beta_array;
			
		}
		
		

	}	// end magOn == 0 block
	


}

double Locality::peierlsPhase(double x1, double x2, double y1, double y2, double B_in){
	double preFactor = 3.14e-5; // This should be changed to be equal to (2 Pi)/(Flux Quanta), with Flux Quanta = h/e, in units of (T*Ang^2)^-1
	double phase = preFactor*B_in*(x2 - x1)*(0.5)*(y1 + y2);
	return phase;
}

double Locality::onSiteE(double x, double y, double z, double E_in){
	// E a constant in the Z direction, return z*E (physical prefactors missing!!!)
	return z*E_in;
}

std::vector< std::vector<double> > Locality::getReciprocal(std::vector< std::vector<double> > a){

	std::vector< std::vector<double> > b;
	b.resize(3);
	for (int i = 0; i < 3; ++i)
		b[i].resize(3);
	
	double denom1 = 0.0;
	double denom2 = 0.0;
	double denom3 = 0.0;
	
	for (int j = 0; j < 3; ++j) {
		denom1 += a[0][j]*crossProd(a[1],a[2],j);
		denom2 += a[1][j]*crossProd(a[2],a[0],j);
		denom3 += a[2][j]*crossProd(a[0],a[1],j);
	}
	
	for (int k = 0; k < 3; ++k){
		b[0][k] = 2*M_PI*crossProd(a[1],a[2],k)/denom1;
		b[1][k] = 2*M_PI*crossProd(a[2],a[0],k)/denom2;
		b[2][k] = 2*M_PI*crossProd(a[0],a[1],k)/denom3;
	}
	
	/*
	for (int l = 0; l < 3; ++l){
		printf("b[%d] = [%lf, %lf, %lf] \n", l, b[l][0], b[l][1], b[l][2]);
	}
	*/

	return b;
}

double Locality::crossProd(std::vector<double> x, std::vector<double> y, int dim){

	if (dim == 0){
		return ( x[1]*y[2] - x[2]*y[1] );
	} else if (dim == 1) {
		return ( x[2]*y[0] - x[0]*y[2] );
	} else if (dim == 2) {
		return ( x[0]*y[1] - x[1]*y[0] );
	} else {
		return 0;
	}

}

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
