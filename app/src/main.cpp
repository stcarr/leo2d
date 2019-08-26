/*
 * File:   main.cpp
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:45 PM
 */

#include "locality.h"
#include "materials/materials.h"

#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>



using namespace std;

int main(int argc, char** argv) {

	// ------------------------------
	// Generate input for simulation.
	// ------------------------------

	// -----------------------------------------
	// First set basic information about the job
	// -----------------------------------------

	// Determines the grid size (from min_size to max_size) which the simulation attempts to populate using a geometric condition (currently checks r < max_size)
	int min_size = -50;
	int max_size = 50;
	int boundary_condition = 0;
	std::vector<int> min;
	std::vector<int> max;
	double max_R;

	// Number of sheets in simulation, s_data,heights,angles determines their properties
	int num_sheets = 0;
	int current_sheet = -1;
	vector<Sdata> s_data;
	vector<double> heights, angles;

	// Supercell parameters (for periodic  boundary conditions)
	double sc_a = 1.0;
	vector< vector<double> > supercell;
	supercell.resize(2);
	supercell[0].resize(2);
	supercell[1].resize(2);

	// ---------------------------------------------------------
	// Next three Categories define a single sheet's information
	// ---------------------------------------------------------

	// Height (in angstroms) and twist angle (CCW, in radians)
	double height = 0.;
	double angle = 0.;

	// File name for the strained position or configuration data
	std::string strain_file;

	// ------------------------------------------------------
	// Solver methods and b-shift options in Job_params class
	// ------------------------------------------------------

	Job_params opts;

	// -----------------------------------------------------------
	// Now we parse the command-line input file for these settings
	// -----------------------------------------------------------

	string line;
	ifstream in_file;
	in_file.open(argv[1]);
	if (in_file.is_open())
	{
		while ( getline(in_file,line) )
		{

			istringstream in_line(line);
			string in_string;

			while ( getline(in_line, in_string, ' ') )	{

				if (in_string == "JOB_NAME"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("job_name", in_string);
					}


				if (in_string == "MAXSIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					max_size = atoi(in_string.c_str());

					max.push_back(max_size);
					max.push_back(max_size);
					max.push_back(max_size);
					max_R = (double) max_size;
					}

				if (in_string == "MINSIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					min_size = atoi(in_string.c_str());

					min.push_back(min_size);
					min.push_back(min_size);
					min.push_back(min_size);
					}

					if (in_string == "MAX_R"){
						getline(in_line,in_string,' ');
						getline(in_line,in_string,' ');
						max_R = atof(in_string.c_str());
						}


				if (in_string == "BOUNDARY_CONDITION"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					boundary_condition = atoi(in_string.c_str());
					opts.setParam("boundary_condition",boundary_condition);
				}


				if (in_string == "SUPERCELL_GROUPS"){
						int sheet_count = 0;
						std::vector<int> sc_groups;
						//int num_mom_groups = opts.getInt("num_mom_groups");
						int num_sc_groups = 2;
						sc_groups.resize(num_sheets);
						getline(in_line,in_string,' ');
						for (int g_count = 0; g_count < num_sc_groups; ++g_count){
							int same_group = 1;
							while (same_group && (sheet_count < num_sheets) ){
								getline(in_line,in_string,' ');
								std::string temp_string = in_string;
								if (in_string == "/"){
									same_group = 0;
								} else {
									// subtract one to put sheets into 0-based indexing for rest of code
									printf("adding sheet %d to group %d \n",stoi(temp_string),g_count+1);
									sc_groups[(stoi(temp_string) - 1)] = g_count;
									sheet_count = sheet_count + 1;
								}
							}

						}
						opts.setParam("sc_groups",sc_groups);
				}

				if (in_string == "SUPERCELL_GRID"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					int m = atoi(in_string.c_str());
					getline(in_line,in_string,' ');
					int n = atoi(in_string.c_str());
					opts.setParam("x_supercell",m);
					opts.setParam("y_supercell",n);
					int two = 2;
					opts.setParam("supercell_type",two);
				}

				if (in_string == "SUPERCELL_SEARCH_SIZE") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("sc_search_size",atoi(in_string.c_str()));
				}

				if (in_string == "K_SAMPLING") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("k_sampling",atoi(in_string.c_str()));
				}

				if (in_string == "K_TYPE") {
					// 0 -> Grid sampling of layer 1's Hexagonal BZ
					// 1 -> LC sampling of twisted supercell Hexagonal BZ
					// 2 -> LC sampling of layer 1's Hexagonal BZ ()
					// 3 -> LC sampling for sandwich project (twisted supercell, both K,K')
					// 4 -> LC sampling of trilayer graphene supercell (3 incomm twist angles)
					// 5 -> LC sampling of trilayer graphene through each K point (3 incomm twist angles)

					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("k_type",atoi(in_string.c_str()));
				}

				if (in_string == "K_GRID") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("num_k1",atoi(in_string.c_str()));
					getline(in_line,in_string,' ');
					opts.setParam("num_k2",atoi(in_string.c_str()));
				}

				if (in_string == "MAT_FROM_FILE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mat_from_file", atoi(in_string.c_str()));
				}

				if (in_string == "MAT_FILENAME"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mat_filename", in_string);
				}

				if (in_string == "NUM_SHEETS") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					num_sheets = atoi(in_string.c_str());
					s_data.resize(num_sheets);
					heights.resize(num_sheets);
					angles.resize(num_sheets);
					opts.setParam("num_sheets",num_sheets);
					}

				if (in_string == "SUPERCELL_M_N"){
				    getline(in_line,in_string,' ');
				    getline(in_line,in_string,' ');
				    int m = atoi(in_string.c_str());
				    getline(in_line,in_string,' ');
				    int n = atoi(in_string.c_str());
				    opts.setParam("m_supercell",m);
				    opts.setParam("n_supercell",n);
		        int one = 1; // explicit type cast, saftey for setParam
		        opts.setParam("supercell_type",one);


				}

				// set to 1 for TMDCs (different a2 relative to graphene definition)
				if (in_string == "HEX_SUPERCELL_MODIFY"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("hex_supercell_modify", atoi(in_string.c_str()));
				}

				if (in_string == "TRILAYER_SUPERCELL"){
				    getline(in_line,in_string,' ');
				    getline(in_line,in_string,' ');
				    int trilayer_on = atoi(in_string.c_str());
						if (trilayer_on == 1) {
			        int three = 3;
			        opts.setParam("supercell_type",three);
						}
				}

				if (in_string == "START_SHEET") {
					getline(in_line,in_string,' ');
					current_sheet = atoi(in_string.c_str()) - 1;
					if (current_sheet+1 > num_sheets){
						printf("ERROR: Cannot take in sheet %d if NUM_SHEETS = %d ! Quitting... \n",current_sheet+1,num_sheets);
						return -1;
					}
				}

				if (in_string == "END_SHEET") {
					getline(in_line,in_string,' ');
					if (current_sheet == atoi(in_string.c_str()) - 1) {
						// last two entries are 0, for solver_type and strain_type. They may be set later in the input file
						heights[current_sheet] = height;
						angles[current_sheet] = angle;
						current_sheet = -1;
					}
				}

				if (in_string == "MATERIAL") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					//if (in_string == "FILE"){
					if (opts.getInt("mat_from_file") == 1){
						//opts.setParam("mat_from_file", 1);
						Materials::Mat blank_mat;
						s_data[current_sheet] = Sdata(blank_mat, min,max,max_R,boundary_condition,0,0,strain_file,1);
						s_data[current_sheet].lmat_name = in_string;
					} else {
						Materials::Mat mat = Materials::string_to_mat(in_string);
						s_data[current_sheet] = Sdata(mat,min,max,max_R,boundary_condition,0,0,strain_file,0);
					}
				}

				if (in_string == "HEIGHT"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					height = atof(in_string.c_str());
				}

				if (in_string == "ANGLE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					// In Degrees:
					angle = (2*M_PI*atof(in_string.c_str()))/(360.0);
					// In Radians:
					//angle = atof(in_string.c_str());
				}


				if (in_string == "STRAIN_FILE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					strain_file = in_string;
					opts.setParam("strain_file",strain_file);
				}

				// for the Fourier decomposed coefficients of relaxation (on a supercell)
				if (in_string == "STRAIN_THETAS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					strain_file = in_string;
					opts.setParam("strain_thetas",strain_file);
				}

				if (in_string == "STRAIN_X_COEFFS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					strain_file = in_string;
					opts.setParam("strain_x_coeffs",strain_file);
				}

				if (in_string == "STRAIN_Y_COEFFS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					strain_file = in_string;
					opts.setParam("strain_y_coeffs",strain_file);
				}

				if (in_string == "STRAIN_Z_COEFFS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					strain_file = in_string;
					opts.setParam("strain_z_coeffs",strain_file);
				}

				if (in_string == "NSHIFTS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("nShifts",atoi(in_string.c_str()));
				}

				if (in_string == "ENERGY_RESCALE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("energy_rescale",atof(in_string.c_str()));
				}

				if (in_string == "ENERGY_SHIFT"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("energy_shift",atof(in_string.c_str()));
				}

				if (in_string == "POLY_ORDER"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("poly_order",atoi(in_string.c_str()));
				}

				if (in_string == "SOLVER_TYPE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					if (in_string == "SQ"){
						opts.setParam("solver_type",1);
					} else if (in_string == "LC"){
						opts.setParam("solver_type",2);
					} else if (in_string == "MLMC_VAC"){
						opts.setParam("solver_type",3);
					} else if (in_string == "VD_FILE"){
						opts.setParam("solver_type",4);
					} else if (in_string == "STRAIN"){
						opts.setParam("solver_type",5);
					} else if (in_string == "FIXED_SHIFT"){
						opts.setParam("solver_type",6);
					}
				}

				if (in_string == "NUM_LC_POINTS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("num_lc_points",atoi(in_string.c_str()));
				}
				if (in_string == "LC_POINTS"){

					std::vector< std::vector<double> > lc_points;
					int num_lc_points = opts.getInt("num_lc_points");
					lc_points.resize(num_lc_points);
					getline(in_line,in_string,' ');
					for (int lc_count = 0; lc_count < num_lc_points; ++lc_count){
						lc_points[lc_count].resize(2);
						getline(in_line,in_string,' ');
						lc_points[lc_count][0] = atof(in_string.c_str());
						getline(in_line,in_string,' ');
						lc_points[lc_count][1] = atof(in_string.c_str());
					}

					opts.setParam("lc_points",lc_points);

				}
				if (in_string == "OBSERVABLE_TYPE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					if (in_string == "DOS"){
						opts.setParam("observable_type",0);
					} else if (in_string == "COND"){
						opts.setParam("observable_type",1);
					}
				}

				if (in_string == "DOS_TRANSFORM"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
				        opts.setParam("dos_transform",atoi(in_string.c_str()));
				}

				if (in_string == "COND_TRANSFORM"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
				        opts.setParam("cond_transform",atoi(in_string.c_str()));
				}

				if (in_string == "COND_POLY_PAR"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
				        opts.setParam("cond_poly_par",atoi(in_string.c_str()));
				}

				if (in_string == "COND_POLY_PAR_SCALE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
				        opts.setParam("cond_poly_par_scale",atoi(in_string.c_str()));
				}

				if (in_string == "SOLVER_SPACE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					if (in_string[0] == 'R'){
						opts.setParam("solver_space",0);
					} else if (in_string[0] == 'M'){
						opts.setParam("solver_space",1);
					}
				}

				if (in_string == "MOM_VF_ONLY"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("mom_vf_only",atoi(in_string.c_str()));
				}

				if (in_string == "MOM_VF_TYPE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("mom_vf_type",atoi(in_string.c_str()));
				}

				if (in_string == "NUM_MOM_GROUPS"){
						getline(in_line,in_string,' ');
						getline(in_line,in_string,' ');
						opts.setParam("num_mom_groups",atoi(in_string.c_str()));
				}

				if (in_string == "MOM_GROUPS"){
						int sheet_count = 0;
						std::vector< std::vector<int> > mom_groups;
						int num_mom_groups = opts.getInt("num_mom_groups");
						mom_groups.resize(num_mom_groups);
						getline(in_line,in_string,' ');
						for (int mom_count = 0; mom_count < num_mom_groups; ++mom_count){
							std::vector<int> curr_mom_group;
							int same_group = 1;
							while (same_group && (sheet_count < num_sheets) ){
								getline(in_line,in_string,' ');
								std::string temp_string = in_string;
								if (in_string == "/"){
									same_group = 0;
								} else {
									// subtract one to put sheets into 0-based indexing for rest of code
									printf("adding sheet %d to group %d \n",stoi(temp_string),mom_count+1);
									curr_mom_group.push_back(stoi(temp_string) - 1);
									sheet_count = sheet_count + 1;
								}
							}

							mom_groups[mom_count] = curr_mom_group;

						}
						opts.setParam("mom_groups",mom_groups);
				}



				if (in_string == "FFT_FROM_FILE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_from_file",atoi(in_string.c_str()));
				}

				if (in_string == "FFT_FILENAME"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("fft_file", in_string);
				}

				if (in_string == "FFT_N_X"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_n_x",atoi(in_string.c_str()));
				}

				if (in_string == "FFT_N_Y"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_n_y",atoi(in_string.c_str()));
				}

				if (in_string == "FFT_L_X"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_L_x",atoi(in_string.c_str()));
				}

				if (in_string == "FFT_L_Y"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_L_y",atoi(in_string.c_str()));
				}

				if (in_string == "FFT_LENGTH_X"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_length_x",atoi(in_string.c_str()));
				}

				if (in_string == "FFT_LENGTH_Y"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
	        opts.setParam("fft_length_y",atoi(in_string.c_str()));
				}

				if (in_string == "STRAIN_TYPE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');

					if (in_string[0] == 'N'){ // STRAIN_TYPE = NONE
						opts.setParam("strain_type",0);
					} else if (in_string[0] == 'S'){ // STRAIN_TYPE = SUPERCELL
						opts.setParam("strain_type",1);
					} else if (in_string[0] == 'C'){ // STRAIN_TYPE = CONFIGURATION
						opts.setParam("strain_type",2);
					} else if (in_string[0] == 'R'){ // STRAIN_TYPE = REALSPACE
						opts.setParam("strain_type",3);
					} else if (in_string[0] == 'F'){ // STRAIN_TYPE = POSITIONS FROM FILE
						opts.setParam("strain_type",4);
					} else if (in_string[0] == 'P'){ // STRAIN_TYPE = PLANE WAVES (Fourier comp. in Config space)
						opts.setParam("strain_type",5);
					} else if (in_string[0] == 'D'){ // STRAIN_TYPE = DIRECT PLANE WAVES (Fourier comp. on supercell)
						opts.setParam("strain_type",6);
					}
				}

				if (in_string == "STRAIN_LAMBDA"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("strain_lambda",atof(in_string.c_str()));
				}

				if (in_string == "UNIFORM_STRAIN"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("uniform_strain",atoi(in_string.c_str()));
				}

				if (in_string == "UNIFORM_STRAIN_MAT"){
					// expects:
					// UNIFORM_STRAIN_MAT = u_xx u_xy u_yx u_yy
					getline(in_line,in_string,' '); // skip "="
					std::vector< std::vector<double> > u_ij;
					u_ij.resize(2);
					u_ij[0].resize(2);
					u_ij[1].resize(2);
					for (int i = 0; i < 2; ++i){
						for (int j = 0; j < 2; ++j){
							getline(in_line,in_string,' ');
							u_ij[i][j] = atof(in_string.c_str());
						}
					}
					opts.setParam("u_ij",u_ij);
				}

				if (in_string == "GSFE_Z_ON"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("gsfe_z_on",atoi(in_string.c_str()));
				}

				if (in_string == "GSFE_R_SCALE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("gsfe_r_scale",atof(in_string.c_str()));
				}

				if (in_string == "MATRIX_SAVE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("matrix_save",atoi(in_string.c_str()));
				}

				if (in_string == "MATRIX_POS_SAVE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("matrix_pos_save",atoi(in_string.c_str()));
				}

				if (in_string == "VERBOSE_SAVE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("verbose_save",atoi(in_string.c_str()));
				}

				if (in_string == "DIAGONALIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("diagonalize",atoi(in_string.c_str()));
				}

				if (in_string == "D_TYPE"){ // 0-> full eigenvalue and eigenvec solve. || 1 -> only eigenvalues || 2-> only the 8 eigenvalues near center
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("d_type",atoi(in_string.c_str()));

				}

				if (in_string == "D_KPM_DOS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("d_kpm_dos",atoi(in_string.c_str()));
				}

				if (in_string == "KPM_TRACE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("kpm_trace",atoi(in_string.c_str()));
				}

				if (in_string == "KPM_TRACE_SAMPS") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("kpm_trace_samps",atoi(in_string.c_str()));
				}

				if (in_string == "D_WEIGHTS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("d_weights",atoi(in_string.c_str()));
				}

				if (in_string == "D_VECS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("d_vecs",atoi(in_string.c_str()));
				}

				if (in_string == "D_COND"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("d_cond",atoi(in_string.c_str()));
				}

				if (in_string == "CHIRAL_ON"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("chiral_on",atoi(in_string.c_str()));
				}

				// MLMC parameters

				if (in_string == "MLMC"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc",atoi(in_string.c_str()));
				}

				if (in_string == "MLMC_MAX_LEVEL"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc_max_level",atoi(in_string.c_str()));
				}
				if (in_string == "MLMC_LEVEL"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc_level",atoi(in_string.c_str()));
				}
				if (in_string == "MLMC_NUM_CLUSTERS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc_num_clusters",atoi(in_string.c_str()));
				}
				if (in_string == "MLMC_CLUSTER_SIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc_cluster_size",atoi(in_string.c_str()));
				}


				if (in_string == "MLMC_OUT_ROOT"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc_out_root", in_string);
					}

				if (in_string == "MLMC_TEMP_ROOT"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("mlmc_temp_root", in_string);
					}

				if (in_string == "VACANCY_FILE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("vacancy_file", in_string);
					}

				// Wannierzation parameters

				if (in_string == "WAN_SAVE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("wan_save",atoi(in_string.c_str()));
				}

				if (in_string == "WAN_NUM_BANDS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("wan_num_bands",atoi(in_string.c_str()));
				}

				// Ballsitic Transport parameters

				if (in_string == "BALLISTIC_TRANSPORT"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("ballistic_transport",atoi(in_string.c_str()));
				}

				if (in_string == "BALLISTIC_SIGMA"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("ballistic_sigma",atof(in_string.c_str()));
				}

				if (in_string == "BALLISTIC_TIME_STEP"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("ballistic_time_step",atof(in_string.c_str()));
				}

				if (in_string == "BALLISTIC_MAX_STEPS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("ballistic_max_steps",atoi(in_string.c_str()));
				}

				// E,B field parameters

				if (in_string == "USE_B_FIELD"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("magOn",atoi(in_string.c_str()));
				}

				if (in_string == "B_FIELD"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("B",atof(in_string.c_str()));
				}

				if (in_string == "USE_E_FIELD"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("elecOn",atoi(in_string.c_str()));
				}

				if (in_string == "E_FIELD"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("E",atof(in_string.c_str()));
				}

				if (in_string == "VACANCY_CHANCE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("vacancy_chance",atof(in_string.c_str()));
				}

				if (in_string == "NUM_TARGET_SHEETS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("num_target_sheets",atoi(in_string.c_str()));
				}

				if (in_string == "TARGET_SHEETS"){
					std::vector<int> temp_sheets;
					temp_sheets.resize(opts.getInt("num_target_sheets"));
					getline(in_line,in_string,' ');
					for (int i = 0; i < opts.getInt("num_target_sheets"); ++i){
						getline(in_line,in_string,' ');
						int temp_target_sheet = atoi(in_string.c_str()) - 1;
						temp_sheets[i] = temp_target_sheet;
					}
					opts.setParam("target_sheets",temp_sheets);
				}

				if (in_string == "NUM_SHIFT_SHEETS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("num_shift_sheets",atoi(in_string.c_str()));
				}

				if (in_string == "SHIFT_SHEETS"){
					std::vector<int> temp_sheets;
          temp_sheets.resize(opts.getInt("num_shift_sheets"));
					getline(in_line,in_string,' ');
					for (int i = 0; i < opts.getInt("num_shift_sheets"); ++i){
						getline(in_line,in_string,' ');
						int temp_target_sheet = atoi(in_string.c_str()) - 1;
						temp_sheets[i] = temp_target_sheet;
					}
					opts.setParam("shift_sheets",temp_sheets);
				}

				if (in_string == "GLOBAL_SHIFTS_ON"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("global_shifts_on",atoi(in_string.c_str()));
				}

				if (in_string == "GLOBAL_SHIFTS"){
					std::vector< std::vector<double> > temp_shifts;
          temp_shifts.resize(opts.getInt("num_sheets"));
					for (int i = 0; i < opts.getInt("num_sheets"); ++i){
						getline(in_line,in_string,' '); // skip "=" or "|" spacing
						temp_shifts[i].resize(3);
						for (int d = 0; d < 3; ++d){
							getline(in_line,in_string,' ');
							double var = atof(in_string.c_str());
							temp_shifts[i][d] = var;
						}
						printf("global_shifts[%d] = [%lf, %lf, %lf] \n",i,temp_shifts[i][0], temp_shifts[i][1], temp_shifts[i][2] );

					}

					opts.setParam("global_shifts",temp_shifts);
				}

			}
		}
		in_file.close();

		if (opts.getInt("poly_order")%4 != 0){
			printf("!!WARNING!!: poly_order = %d is NOT divisible 4 (needed for KPM iterative method) \nLEO2D Quiting... \n",opts.getInt("poly_order"));
			return -1;
		}

		if (opts.getInt("solver_type") == 3 && opts.getInt("nShifts")%2 == 0){
			opts.setParam("nShifts",opts.getInt("nShifts") + 1);
			printf("!!WARNING!!: Setting nShifts to an odd number for the vacancy sweep method! nShifts = %d \n",opts.getInt("nShifts"));
		}

		/*
		if (num_sheets > 1 && boundary_condition == 1){
			printf("!!WARNING: multiple sheet periodic run detected! Periodic boundary conditions not yet implemented for interlayer coupling! \n");
			return -1;
		}
		*/

		// update solver_space (i.e. sheets need to know solver_space, but no gaurentee it was set before sdata were input)
		for (int i = 0; i < num_sheets; ++i){
			s_data[i].solver_space = opts.getInt("solver_space");
			s_data[i].strain_type = opts.getInt("strain_type");
			s_data[i].sheet_index = i;
			s_data[i].opts = opts;
		}

		// Create the locality object with the sheet input data
		Locality loc(s_data,heights,angles);

		// Start MPI within Locality object on each processor
		int multi_rank_job = loc.initMPI(argc, argv);
		if (multi_rank_job == -1){
			printf("Error: Only 1 MPI rank detected (need to run with n > 1).\n");
			loc.finMPI();
			return -1;
		}

		// Simulation's solver is set with setup call to Locality object
		loc.setup(opts);

		// Builds the geometry of the problem on the root node and then sends them to workers via MPI
		loc.constructGeom();

		// Post processing operations. Save prints timing information from each node. File saves happen on the MPI loop from the root node!
		loc.save();

		// End MPI processes and finish
		loc.finMPI();


	} else {
		printf("Input file '%s' not found, stopping, \n",argv[1]);

	}

	// No errors hopefully!
	return 0;

}
