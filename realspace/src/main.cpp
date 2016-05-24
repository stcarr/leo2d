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
	
	// Gets pre-pended to all output files
	std::string job_name = "HSTRUCT_JOB";
	
	// Determines the grid size (from min_size to max_size) which the simulation attempts to populate using a geometric condition (currently checks r < max_size)
	int min_size = -50;
	int max_size = 50;
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
	
	// ----------------------------------
	// Solver methods and b-shift options
	// ----------------------------------
	
	// Number of b-shifts to perform in one direction (uniform grid sample over the first sheets unit cell is performed via MPI)
	int nShifts = 1;
	
	// Solver information (type = 0 is FILTLAN local eigensolve, = 2 is Chebyshev spectrum sample.
	int solver_type = 0;
	double interval_start = -1;
	double interval_end = 1;
	
	// FILTLAN settings
	int num_eigs = 1;
	
	// Chebyshev settings
	int num_samples = 1;
	double energy_rescale = 15;
	double energy_shift = 0;
	double cheb_width = 0.2;
	int poly_order = 3000;
	
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
					job_name = in_string;
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
						s_data[current_sheet] = Sdata(unitCell,types,pos,min,max,mat);
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
				
				if (in_string == "NSHIFTS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					nShifts = atoi(in_string.c_str());
				}
				
				if (in_string == "NUM_EIGS"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					num_eigs = atoi(in_string.c_str());
				}
				if (in_string == "NUM_SAMPLES"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					num_samples = atoi(in_string.c_str());
				}

				if (in_string == "INTERVAL_START"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					interval_start = atof(in_string.c_str());
				}

				if (in_string == "INTERVAL_END"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					interval_end = atof(in_string.c_str());
				}
				
				if (in_string == "ENERGY_RESCALE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					energy_rescale = atof(in_string.c_str());
				}
				
				if (in_string == "ENERGY_SHIFT"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					energy_shift = atof(in_string.c_str());
				}
				
				if (in_string == "CHEB_WIDTH"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					cheb_width = atof(in_string.c_str());
				}
				
				if (in_string == "POLY_ORDER"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					poly_order = atoi(in_string.c_str());
				}
				
				if (in_string == "SOLVER_TYPE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					solver_type = atoi(in_string.c_str());
				}				
				
				
			}
		}
		in_file.close();
	}
	
	// Create the locality object with the sheet input data
	Locality loc(s_data,heights,angles);
	
	// Simulation's solver is set with setup call to Locality object
	loc.setup(job_name,nShifts, num_eigs, num_samples, interval_start, interval_end, energy_rescale, energy_shift, cheb_width, poly_order, solver_type);
	
	// Start MPI within Locality object on each processor
	loc.initMPI(argc, argv);
	
	// Builds the geometery of the problem on the root node and then sends them to workers via MPI
	loc.constructGeom();
	
	// Oprn post processing operations. Plot does nothing currently, save prints timing information from each node. File saves happen on the MPI loop from the root node!
	loc.plot();
	loc.save();
	
	// End MPI processes and finish
	loc.finMPI();
	
	// No errors hopefully!
	return 0;

}
