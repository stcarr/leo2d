/* 
 * File:   main.cpp
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:45 PM
 */

#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include "locality.h"

#include <mpi.h>

#include <petscksp.h>
#include <slepceps.h>



using namespace std;
 
int main(int argc, char** argv) {

	// Generate input for simulation.
	//
	// Eventually will want this read from an input file
	// so that we don't need to recompile it each time

	
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
	unitCell.push_back(a3);
	
	// Number and atomic weight of atoms
	int num_atoms = 2;
	vector<int> types;
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
	

	// Shape of the grid to sample for atom population
	int min_size = -2;
	int max_size = 2;
	
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
	vector<Sdata> sdata;
	
	sdata.push_back(s_in);
	sdata.push_back(s_in);
	
	heights.push_back(0);
	heights.push_back(3);
	
	angles.push_back(0);
	angles.push_back(1);
	
	
	// Create the locality object with the input data
	Locality loc(sdata,heights,angles);
	
	// Simulation is controlled with following calls to Locality object
	
	// Start MPI within Locality object on each processor
	loc.initMPI(argc, argv);
	
	// Currently empty function, when written will allow us to do 
	// multiple runs without resetting MPI 
	loc.setup();
	
	
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
