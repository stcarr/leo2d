/* 
 * File:   main.cpp
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:45 PM
 */

#include "locality.h"

//#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <sstream>
//#include <mpi.h>

//#include <petscksp.h>
//#include <slepceps.h>



using namespace std;
 
int main(int argc, char** argv) {

	// char help[] = "what this program does in brief can go here.";

	// PetscInitialize(&argc,&argv,(char*)0,help);
	//SlepcInitialize(&argc,&argv,(char*)0,help);

	// Generate input for simulation.
	//
	// Eventually will want this read from an input file
	// so that we don't need to recompile it each time

	
	// Shape of the grid to sample for atom population
	int min_size = -300;
	int max_size = 300;
	std::vector<int> min;
	std::vector<int> max;
	
	
	int num_sheets = 0;
	int current_sheet = -1;
	vector<Sdata> s_data;
	vector<double> heights, angles;

	double a = 0;
	std::vector<std::vector<double> > unitCell;
	unitCell.resize(3);
	for (int i = 0; i < 3; ++i)
		unitCell[i].resize(3);

	int num_orbitals = 0;
	std::vector<int> types;
	std::vector<std::vector<double> > pos;
	
	double height = 0;
	double angle = 0;
	
	int nShifts = 1;
	int num_eigs = 1;
	int num_samples = 1;
	double interval_start = -1;
	double interval_end = 1;
	int solver_type = 0;
	
	string line;
	ifstream in_file;
	in_file.open(argv[1]);
	if (in_file.is_open())
	{
		while ( getline(in_file,line) )
		{
			//cout << line << "\n";
			
			istringstream in_line(line);
			string in_string;
			
			while ( getline(in_line, in_string, ' ') )	{
			
				//cout << in_string << "\n";
				
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
						s_data[current_sheet] = Sdata(unitCell,types,pos,min,max);
						heights[current_sheet] = height;
						angles[current_sheet] = angle;
					}
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
						unitCell[0][i] = a*atof(in_string.c_str());
					}
				}
				
				if (in_string == "UNITCELL2"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < 3; ++i) {
						getline(in_line,in_string,' ');
						unitCell[1][i] = a*atof(in_string.c_str());
					}
				}
				
				if (in_string == "UNITCELL3"){
					getline(in_line,in_string,' ');
					for (int i = 0; i < 3; ++i) {
						getline(in_line,in_string,' ');
						unitCell[2][i] = a*atof(in_string.c_str());
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
							pos[i][j] = atof(in_string.c_str());
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
				
				if (in_string == "SOLVER_TYPE"){
					getline(in_line,in_string,' ');
					getline(in_line,in_string,' ');
					solver_type = atoi(in_string.c_str());
				}				
				
				
			}
		}
		in_file.close();
	}
	
	/*
	// Unit Cell: is just test information, basic square lattice
	std::vector<std::vector<double> > unitCell;
	std::vector<double> a1,a2,a3;

	double a = 2.46;
	
	double unitCell_in[3][3] = 
	{
		{a,0,0},
		{a/2.0,a*sqrt(3)/2.0,0},
		{0,0,1}
	};
	
	for (int i = 0; i < 3; i++){
		a1.push_back(unitCell_in[0][i]);
		a2.push_back(unitCell_in[1][i]);
		a3.push_back(unitCell_in[2][i]);
	}
	
	unitCell.push_back(a1);
	unitCell.push_back(a2);
	vector<int> types;
	unitCell.push_back(a3);
	
	// Number and atomic weight of atoms
	int num_atoms = 2;
	types.push_back(6); // 6 is carbon (atomic #)
	types.push_back(6);
	
	// Position of atoms in the Unit Cell
	vector<vector<double> > pos;
	pos.resize(num_atoms);
	for (int i = 0; i < num_atoms; ++i)
		pos[i].resize(3);

	pos[0][0] = 0.0;
	pos[0][1] = 0.0;
	pos[0][2] = 0.0;
	
	pos[1][0] = a/2.0;
	pos[1][1] = a/(2.0*sqrt(3));
	pos[1][2] = 0;
	
	std::vector<int> min;
	min.push_back(min_size);
	min.push_back(min_size);
	min.push_back(min_size);
	
	std::vector<int> max;
	max.push_back(max_size);
	max.push_back(max_size);
	max.push_back(max_size);
	
	// Save all this info into a wrapper data type
	Sdata s_in(unitCell,types,pos,min,max);
	
	// Heterostructure input information
	vector<double> heights, angles;
	vector<Sdata> s_data;
	
	s_data.push_back(s_in);
	s_data.push_back(s_in);
	
	heights.push_back(0);
	heights.push_back(6);
	
	angles.push_back(0);
	angles.push_back(0.1); // roughly 6 degree angle
	*/
	
	
	// Create the locality object with the input data
	Locality loc(s_data,heights,angles);
	
	// Simulation is controlled with following calls to Locality object
	loc.setup(nShifts, num_eigs, num_samples, interval_start, interval_end, solver_type);
	
	// Start MPI within Locality object on each processor
	loc.initMPI(argc, argv);
	
	// Builds q-points and lists the non-zero elements of the Hamiltonian
	// on the root node and then sends them to workers via MPI
	//
	// !! Need to also make it send b-shift information to each worker !!
	loc.constructGeom();
	
	// Post processing operations
	//
	// !! NOT STARTED !!
	loc.plot();
	loc.save();
	
	// End MPI processes and finish
	loc.finMPI();
	return 0;

}
