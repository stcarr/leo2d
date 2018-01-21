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
					}

				if (in_string == "MINSIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					min_size = atoi(in_string.c_str());

					min.push_back(min_size);
					min.push_back(min_size);
					min.push_back(min_size);
					}

				if (in_string == "BOUNDARY_CONDITION"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					boundary_condition = atoi(in_string.c_str());
					opts.setParam("boundary_condition",boundary_condition);
				}

				if (in_string == "SUPERCELL_M_N"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					int m = atoi(in_string.c_str());
					getline(in_line,in_string,' ');
					int n = atoi(in_string.c_str());
					opts.setParam("m_supercell",m);
					opts.setParam("n_supercell",n);
					int o = 1;
					opts.setParam("supercell_type",o);
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

				if (in_string == "K_SAMPLING") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("k_sampling",atoi(in_string.c_str()));
				}

				if (in_string == "K_TYPE") {
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

				if (in_string == "NUM_SHEETS") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					num_sheets = atoi(in_string.c_str());
					s_data.resize(num_sheets);
					heights.resize(num_sheets);
					angles.resize(num_sheets);
					opts.setParam("num_sheets",num_sheets);
					}

				if (in_string == "START_SHEET") {
					getline(in_line,in_string,' ');
					current_sheet = atoi(in_string.c_str()) - 1;
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
					Materials::Mat mat = Materials::string_to_mat(in_string);
					s_data[current_sheet] = Sdata(mat,min,max,boundary_condition,0,0,strain_file);
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
					}
				}

				if (in_string == "STRAIN_LAMBDA"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("strain_lambda",atof(in_string.c_str()));
				}

				if (in_string == "GSFE_Z_ON"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("gsfe_z_on",atoi(in_string.c_str()));
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

				if (in_string == "D_KPM_DOS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("d_kpm_dos",atoi(in_string.c_str()));
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
		}

		// update supercell if BC
		if (boundary_condition == 1){

			int type = opts.getInt("supercell_type");

			// Standard fixed supercell definition
			if (type == 0) {
				std::vector< std::vector<double> > sc = opts.getDoubleMat("supercell");
				for (int i = 0; i < s_data.size(); ++i){
					double theta = angles[i];
					std::vector< std::vector<double> > sc_here;
					sc_here.resize(2);
					sc_here[0].resize(2);
					sc_here[1].resize(2);

					sc_here[0][0] =  sc[0][0]*cos(theta) + sc[0][1]*sin(theta);
					sc_here[0][1] = -sc[0][0]*sin(theta) + sc[0][1]*cos(theta);
					sc_here[1][0] =  sc[1][0]*cos(theta) + sc[1][1]*sin(theta);
					sc_here[1][1] = -sc[1][0]*sin(theta) + sc[1][1]*cos(theta);
					s_data[i].supercell = sc_here;

				}
			} else if (type == 1) { // (M,N) Supercell type
				if ((int)s_data.size() > 2){
					printf("!!WARNING!!: (M,N) Supercell method not defined for more than 2 sheets \nLEO2D Quiting... \n");
					return -1;
				}
				int M = opts.getInt("m_supercell");
				int N = opts.getInt("n_supercell");

				double theta = acos((N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
				if (M < N){
					theta = -theta;
				}
				printf("supercell theta = %lf degrees (acos(%lf) )\n",360.0*theta/(2.0*M_PI), (N*N + 4*N*M + M*M)/(2.0*(N*N + N*M + M*M)));
				// we assume unitCell is same for both sheets...
				for (int i = 0; i < s_data.size(); ++i){

					std::vector< std::vector<double> > base_unitCell;
					base_unitCell.resize(2);
					for(int n = 0; n < 2; ++n){
						base_unitCell[n].resize(2);
						for(int m = 0; m < 2; ++m){
							base_unitCell[n][m] = Materials::lattice(s_data[0].mat)[n][m];
						}
					}

					int A1_num_a1;
					int A1_num_a2;
					int A2_num_a1;
					int A2_num_a2;
					if (i == 0){
						printf("i = 0\n");
						angles[0] = 0;
						A1_num_a1 = N;
						A1_num_a2 = M;
						A2_num_a1 = -M;
						A2_num_a2 = (M+N);
					} else if (i == 1){
						printf("i = 1\n");
						angles[1] = theta;
						A1_num_a1 = M;
						A1_num_a2 = N;
						A2_num_a1 = -N;
						A2_num_a2 = (M+N);
					}

					std::vector< std::vector<double> > unitCell;
					unitCell.resize(2);
					unitCell[0].resize(2);
					unitCell[1].resize(2);

					// Not updating the twist-angle is correct, as this will make sure the supercell is always
					// saved in coordinates relative to the individual sheet, not the overall heterostructure.
					double theta_here = 0.0;

					unitCell[0][0] = cos(theta_here)*base_unitCell[0][0] - sin(theta_here)*base_unitCell[0][1];
					unitCell[0][1] = sin(theta_here)*base_unitCell[0][0] + cos(theta_here)*base_unitCell[0][1];
					unitCell[1][0] = cos(theta_here)*base_unitCell[1][0] - sin(theta_here)*base_unitCell[1][1];
					unitCell[1][1] = sin(theta_here)*base_unitCell[1][0] + cos(theta_here)*base_unitCell[1][1];

					std::vector< std::vector<double> > sc_here;
					sc_here.resize(2);
					sc_here[0].resize(2);
					sc_here[1].resize(2);

					sc_here[0][0] = A1_num_a1*unitCell[0][0] + A1_num_a2*unitCell[1][0];
					sc_here[0][1] = A1_num_a1*unitCell[0][1] + A1_num_a2*unitCell[1][1];
					sc_here[1][0] = A2_num_a1*unitCell[0][0] + A2_num_a2*unitCell[1][0];
					sc_here[1][1] = A2_num_a1*unitCell[0][1] + A2_num_a2*unitCell[1][1];

					std::vector< std::vector<int> > sc_stride_here;
					sc_stride_here.resize(2);
					sc_stride_here[0].resize(2);
					sc_stride_here[1].resize(2);

					sc_stride_here[0][0] = A1_num_a1;
					sc_stride_here[0][1] = A1_num_a2;
					sc_stride_here[1][0] = A2_num_a1;
					sc_stride_here[1][1] = A2_num_a2;

					/*
					std::vector< std::vector<double> > local_sc_here;
					local_sc_here.resize(2);
					local_sc_here[0].resize(2);
					local_sc_here[1].resize(2);

					double theta_here = -angles[i];

					local_sc_here[0][0] = cos(theta_here)*sc_here[0][0] - sin(theta_here)*sc_here[0][1];
					local_sc_here[0][1] = sin(theta_here)*sc_here[0][0] + cos(theta_here)*sc_here[0][1];
					local_sc_here[1][0] = cos(theta_here)*sc_here[1][0] - sin(theta_here)*sc_here[1][1];
					local_sc_here[1][1] = sin(theta_here)*sc_here[1][0] + cos(theta_here)*sc_here[1][1];
					*/
					printf("unitCell  = [%lf %lf; %lf %lf]\n",unitCell[0][0],unitCell[0][1],unitCell[1][0],unitCell[1][1]);
					printf("supercell = [%lf %lf; %lf %lf]\n", sc_here[0][0], sc_here[0][1], sc_here[1][0], sc_here[1][1]);
					//printf("local_supercell = [%lf %lf; %lf %lf]\n", local_sc_here[0][0], local_sc_here[0][1], local_sc_here[1][0], local_sc_here[1][1]);
					printf("sc_stride = [%d %d; %d %d]\n", sc_stride_here[0][0], sc_stride_here[0][1], sc_stride_here[1][0], sc_stride_here[1][1]);


					s_data[i].supercell = sc_here;
					s_data[i].supercell_stride = sc_stride_here;

					// Since sheet[0] has 0 twist, we can use it's supercell as the supercell for hstruct!!
					if (i == 0){
						opts.setParam("supercell",sc_here);

						int matrix_pos_save = opts.getInt("matrix_pos_save");
						if (matrix_pos_save == 1){

							string job_name = opts.getString("job_name");
							const char* extension = "_supercell.dat";

							ofstream outFile;
							outFile.open ((job_name + extension).c_str());
							outFile << 	"0, 0, 0, " << sc_here[0][0] << ", " << sc_here[0][1] << ", 0\n" <<
													"0, 0, 0, " << sc_here[1][0] << ", " << sc_here[1][1] << ", 0\n" <<
													sc_here[0][0] << ", " << sc_here[0][1] << ", 0, " <<
													sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n" <<
													sc_here[1][0] << ", " << sc_here[1][1] << ", 0, " <<
													sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n";
							outFile.close();
						}
					}

				}


			} else if (type == 2) { // (X,Y) Grid type

				int X = opts.getInt("x_supercell");
				int Y = opts.getInt("y_supercell");

				std::vector< std::vector<double> > unitCell;
				unitCell.resize(2);
				for(int n = 0; n < 2; ++n){
					unitCell[n].resize(2);
					for(int m = 0; m < 2; ++m){
						unitCell[n][m] = Materials::lattice(s_data[0].mat)[n][m];
					}
				}

				std::vector< std::vector<double> > sc_here;
				sc_here.resize(2);
				sc_here[0].resize(2);
				sc_here[1].resize(2);

				sc_here[0][0] = X*unitCell[0][0];
				sc_here[0][1] = X*unitCell[0][1];
				sc_here[1][0] = Y*unitCell[1][0];
				sc_here[1][1] = Y*unitCell[1][1];

				std::vector< std::vector<int> > sc_stride_here;
				sc_stride_here.resize(2);
				sc_stride_here[0].resize(2);
				sc_stride_here[1].resize(2);

				sc_stride_here[0][0] = X;
				sc_stride_here[0][1] = 0;
				sc_stride_here[1][0] = 0;
				sc_stride_here[1][1] = Y;

				printf("unitCell  = [%lf %lf; %lf %lf]\n",unitCell[0][0],unitCell[0][1],unitCell[1][0],unitCell[1][1]);
				printf("supercell = [%lf %lf; %lf %lf]\n", sc_here[0][0], sc_here[0][1], sc_here[1][0], sc_here[1][1]);
				printf("sc_stride = [%d %d; %d %d]\n", sc_stride_here[0][0], sc_stride_here[0][1], sc_stride_here[1][0], sc_stride_here[1][1]);

				for (int i = 0; i < (int)s_data.size(); ++i){
					s_data[i].supercell = sc_here;
					s_data[i].supercell_stride = sc_stride_here;
				}

				// Since sheet[0] has 0 twist, we can use it's supercell as the supercell for hstruct!!
				opts.setParam("supercell",sc_here);

				int matrix_pos_save = opts.getInt("matrix_pos_save");
				if (matrix_pos_save == 1){

					string job_name = opts.getString("job_name");
					const char* extension = "_supercell.dat";

					ofstream outFile;
					outFile.open ((job_name + extension).c_str());
					outFile << 	"0, 0, 0, " << sc_here[0][0] << ", " << sc_here[0][1] << ", 0\n" <<
											"0, 0, 0, " << sc_here[1][0] << ", " << sc_here[1][1] << ", 0\n" <<
											sc_here[0][0] << ", " << sc_here[0][1] << ", 0, " <<
											sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n" <<
											sc_here[1][0] << ", " << sc_here[1][1] << ", 0, " <<
											sc_here[0][0]+sc_here[1][0] << ", " << sc_here[0][1]+sc_here[1][1] << ", 0\n";
					outFile.close();
				}

			}
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
