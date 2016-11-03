/* 
 * File:   main.cpp
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:45 PM
 */

#include "locality.h"

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
	int mat = 0;

	// ---------------------------------------------------------
	// Next three Categories define a single sheet's information
	// ---------------------------------------------------------
	
	// Unit cell information, gets put into an sdata object
	double a = 0;
	std::vector<std::vector<double> > unitCell;
	unitCell.resize(3);
	for (int i = 0; i < 3; ++i)
		unitCell[i].resize(3);

	// number of orbitals per unit cell
	int num_orbitals = 0;
	std::vector<int> types;
	std::vector<std::vector<double> > pos;
	
	// Height (in angstroms) and twist angle (CCW, in radians)
	double height = 0;
	double angle = 0;
	
	// ------------------------------------------------------
	// Solver methods and b-shift options in Loc_params class
	// ------------------------------------------------------
	
	Loc_params opts;
	
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
				
				}
					
				if (in_string == "NUM_SHEETS") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					num_sheets = atoi(in_string.c_str());
					s_data.resize(num_sheets);
					heights.resize(num_sheets);
					angles.resize(num_sheets);
					}
				
				if (in_string == "START_SHEET") {
					getline(in_line,in_string,' ');
					current_sheet = atoi(in_string.c_str()) - 1;
				}
				
				if (in_string == "END_SHEET") {
					getline(in_line,in_string,' ');
					if (current_sheet == atoi(in_string.c_str()) - 1) {
						s_data[current_sheet] = Sdata(unitCell,types,pos,min,max,mat,boundary_condition,0);
						heights[current_sheet] = height;
						angles[current_sheet] = angle;
					}
				}
				
				if (in_string == "MATERIAL") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					mat = atoi(in_string.c_str());
				}
				
				if (in_string == "ALPHA") {
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					a = atof(in_string.c_str());
				}
				
				if (in_string == "UNITCELL1"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < 3; ++i) {
						getline(in_line,in_string,' ');
						unitCell[i][0] = a*atof(in_string.c_str());
					}
				}
				
				if (in_string == "UNITCELL2"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < 3; ++i) {
						getline(in_line,in_string,' ');
						unitCell[i][1] = a*atof(in_string.c_str());
					}
				}
				
				if (in_string == "UNITCELL3"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < 3; ++i) {
						getline(in_line,in_string,' ');
						unitCell[i][2] = a*atof(in_string.c_str());
					}
				}
				
				if (in_string == "NUM_ORBITALS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					num_orbitals = atoi(in_string.c_str());
					types.resize(num_orbitals);
					pos.resize(num_orbitals);
				}
				
				if (in_string == "TYPES"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < num_orbitals; ++i) {
						getline(in_line,in_string,' ');
						types[i] = atoi(in_string.c_str());
					}
				}
					
				if (in_string == "POS"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < num_orbitals; ++i) {
						pos[i].resize(3);
						for (int j = 0; j < 3; ++j) {
							getline(in_line,in_string,' ');
							pos[i][j] = a*atof(in_string.c_str());
						}
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
					angle = atof(in_string.c_str());
				}
				
				
				
				if (in_string == "INTRA_SEARCHSIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("intra_searchsize",atoi(in_string.c_str()));
				}
				
				if (in_string == "INTER_SEARCHSIZE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					opts.setParam("inter_searchsize",atoi(in_string.c_str()));
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
					} else if (in_string == "VD_CENTER"){
						opts.setParam("solver_type",3);
					} else if (in_string == "VD_FILE"){
						opts.setParam("solver_type",4);
					}
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
				
				if (in_string == "SOLVER_SPACE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					if (in_string[0] == 'R'){
						opts.setParam("solver_space",0);
					} else if (in_string[0] == 'M'){
						opts.setParam("solver_space",1);
					}
				}
				
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
					int* temp_sheets = new int[opts.getInt("num_shift_sheets")];
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
	}
	
	if (opts.getInt("poly_order")%4 != 0){
		printf("!!WARNING!!: poly_order = %d is NOT divisible 4 (needed for KPM iterative method) \n Quiting... \n",opts.getInt("poly_order"));
		return -1;
	}
	
	if (opts.getInt("solver_type") == 3 && opts.getInt("nShifts")%2 == 0){
		opts.setParam("nShifts",opts.getInt("nShifts") + 1);
		printf("!!WARNING!!: Setting nShifts to an odd number for the vacancy sweep method! nShifts = %d \n",opts.getInt("nShifts"));
	}
	
	if (num_sheets > 1 && boundary_condition == 1){
		printf("!!WARNING: multiple sheet periodic run detected! Periodic boundary conditions not yet implemented for interlayer coupling! \n");
		return -1;
	}
	
	// update solver_space (i.e. sheets need to know solver_space, but no gaurentee it was set before sdata were input)
	for (int i = 0; i < num_sheets; ++i){
		s_data[i].solver_space = opts.getInt("solver_space");
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
	
	// Builds the geometery of the problem on the root node and then sends them to workers via MPI
	loc.constructGeom();
	
	// Post processing operations. Save prints timing information from each node. File saves happen on the MPI loop from the root node!
	loc.save();
	
	// End MPI processes and finish
	loc.finMPI();
	
	// No errors hopefully!
	return 0;

}
