/* 
 * File:   hstruct.cpp
 * Author: Stephen
 * 
 * Created on January 13, 2016, 3:16 PM
 */

#include "hstruct.h"

#include <fftw3.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <math.h>

// ---------------------
// Use this constructor!
// ---------------------
Hstruct::Hstruct(std::vector<Sheet> sheets_in,std::vector<double> angles_in,std::vector<double> heights_in, int solver_space_in) {
    std::vector<double> blank_shift;
	// Create a shift corresponding to 0 so we can initialize
    for(int j = 0; j < 3; ++j)
        blank_shift.push_back(0);
	
	// Save each sheet's information
    for(int i = 0; i < sheets_in.size(); ++i){
        sheets.push_back(sheets_in[i]);
        heights.push_back(heights_in[i]);
        angles.push_back(angles_in[i]);
        shifts.push_back(blank_shift);
    }
    max_sheets = sheets.size();
	solver_space = solver_space_in;
	
	// For Momentum space:
	// Need to have sheet 0 get sheet 1's reciprocal lattice before sheet 0 constructs its geometry (and vice-versa) 
	// Currently we assume that sheet 0 and sheet 1 are identical (up to a twist), so we do not "swap" their reciprocal lattices.
	
	if (solver_space == 1) {
		printf("!!NOTE!!: Momentum-space only implemented for bilayers (check constructor in hstruct.cpp for more info) \n");
	}
    
	// Call to generate index information
	setIndex();
        
    
}

Hstruct::Hstruct(const Hstruct& orig) {
}

Hstruct::~Hstruct() {
}

// ---------------------------------------------------------------
// Figures out global index from each sheets local indexing scheme
// ---------------------------------------------------------------
void Hstruct::setIndex(){
    
	// First figure out global max index
	int k_total = 0;
    for(int s = 0; s < max_sheets; ++s){
        k_total += sheets[s].getMaxIndex();
    }
    
    max_index = k_total;
	
	// Figure out grid position at each orbital index
	index_array.resize(max_index);
	for(int k = 0; k < max_index; ++k){
	
		int i_here, j_here, l_here, s_here;	
		int index_counter = 0;
		
		// Find the sheet that index k is on
		s_here = 0;	
		for(int s = 0; s < sheets.size(); ++s) {
			if (k < index_counter + sheets[s].getMaxIndex()){
				s_here = s;
				break;
			}

			index_counter += sheets[s].getMaxIndex();
		}
		
		// get local grid information from the sheet object
		i_here = sheets[s_here].indexToGrid(k - index_counter, 0);
		j_here = sheets[s_here].indexToGrid(k - index_counter, 1);
		l_here = sheets[s_here].indexToGrid(k - index_counter, 2);
		
		// index array is a 1-d 4*max_index size (use 4*k + 0,1,2,3 for each column)
		index_array[k].push_back(i_here);
		index_array[k].push_back(j_here);
		index_array[k].push_back(l_here);
		index_array[k].push_back(s_here);
		
	}
	
}

int Hstruct::getMaxIndex(){
    return max_index;
}

// ----------------------------------------------------------------------------
// Set shift of a sheet, ends up not being used in this locality implementation
// ----------------------------------------------------------------------------
void Hstruct::setShift(int sheet, std::vector<double> b){
    shifts[sheet] = b;
}

// --------------------------------
// Find the dim position of index k
//  (for dim: x = 0, y = 1, z = 2)
// --------------------------------
double Hstruct::posAtomIndex(int k, int dim){
    if (k > max_index){
        return 0;
    }
    
    bool findSheet = true;
    int current_sheet = 0;
    int current_index = 0;
    int s;

    while(findSheet){
        if (k < (current_index + sheets[current_sheet].getMaxIndex())){
            s = current_sheet;
            findSheet = false;
        }else {
            current_index += sheets[current_sheet].getMaxIndex();
            current_sheet += 1;
        }
    }
    	
    double local_x = sheets[s].posAtomIndex(k - current_index,0);
    double local_y = sheets[s].posAtomIndex(k - current_index,1);
    double local_z = sheets[s].posAtomIndex(k - current_index,2);
	
    double theta;
	
	// If we are in Momentum-space, we need to have sheet 0 have the geometry of sheet 1, so we swap the angles
	
	if (solver_space == 0) {
		theta = angles[s];
	} else if (solver_space == 1) {
		if (s == 0) {
			theta = angles[1];
		} else if (s == 1) {
			theta = angles[0];
		}
	}

	// NEED TO FIX SHIFTS !! (but they are not used in current implementation)
	// should rotate AFTER shifts, not before
	
    if (dim == 0){
        double x = local_x*cos(theta) - local_y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0) + shifts[s][2]*sheets[s].getUnit(2,0);
        return x;
    }
    if (dim == 1){
        double y = local_y*cos(theta) + local_x*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1) + shifts[s][2]*sheets[s].getUnit(2,1);
        return y;
    }
    if (dim == 2){
        double z = local_z + heights[s] + shifts[s][0]*sheets[s].getUnit(0,2) + shifts[s][1]*sheets[s].getUnit(1,2) + shifts[s][2]*sheets[s].getUnit(2,2);
		return z;
    }
}

// --------------------------------------------
// Find the unit cell on sheet s closest to pos
//           for dim: (i = 0, j = 1)
// --------------------------------------------
int Hstruct::findNearest(double (&pos)[3],int s, int dim){

	double theta = angles[s];
	double x = pos[0];
	double y = pos[1];
	double z = pos[2];
	
	double a_inv[2][2];
	
	for (int i = 0; i < 2; ++i){
		for (int j = 0; j < 2; ++j){
			a_inv[i][j] = sheets[s].getInverse(i,j);
		}
	}
	
	// NEED TO FIX SHIFTS (but they are not used in current implementation)!!
	// should rotate AFTER shifts, not before

	if (dim == 0){
		double x_new = x*cos(theta) + y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0) + shifts[s][2]*sheets[s].getUnit(2,0);
		double y_new = y*cos(theta) - x*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1) + shifts[s][2]*sheets[s].getUnit(2,1);
		double i_new = a_inv[0][0]*x_new + a_inv[0][1]*y_new;
		int i = std::min(std::max(int(floor(i_new)),sheets[s].getShape(0,0)),sheets[s].getShape(1,0)) - sheets[s].getShape(0,0);
		return i;
	}
	
	if (dim == 1){
		double x_new = x*cos(theta) + y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0) + shifts[s][2]*sheets[s].getUnit(2,0);
		double y_new = y*cos(theta) - x*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1) + shifts[s][2]*sheets[s].getUnit(2,1);
		double j_new = a_inv[1][0]*x_new + a_inv[1][1]*y_new;
		int j = std::min(std::max(int(floor(j_new)),sheets[s].getShape(0,1)),sheets[s].getShape(1,1)) - sheets[s].getShape(0,1);
		return j;
	}
	
}

std::vector<std::vector<int> > Hstruct::getIndexArray(){
	return index_array;
}

// -----------------------------------------
// Allocates interlayer candidate pairs into 
//     the given 2-dimensional pair_array
// -----------------------------------------
void Hstruct::getInterPairs(std::vector<std::vector<int> > &pair_array, int searchsize){

	// We search over a searchsize x searchsize sized grid of unitcells
	//int searchsize = 6;
	
	// We do not save pairs that are farther apart than this (in Angstroms)
	double inter_cutoff = 12.5;
	
	// loop over all orbitals (kh = "k here")
	for (int kh = 0; kh < max_index; ++kh){
		
		// Get the current grid information (ih = i "here")
		int ih = index_array[kh][0];
		int jh = index_array[kh][1];
		int lh = index_array[kh][2];
		int sh = index_array[kh][3];
		
		// And the position information
		double pos_here[3];
		pos_here[0] = posAtomIndex(kh,0);
		pos_here[1] = posAtomIndex(kh,1);
		pos_here[2] = posAtomIndex(kh,2);
	
		// if we are not on the "lowest" sheet, we look for pairs from the sheet above
		if (sh > 0) {
		
			// "0" determines the center of our search range on sheet s0
			int i0 = findNearest(pos_here, sh - 1, 0);
			int j0 = findNearest(pos_here, sh - 1, 1);
			int s0 = sh - 1;

			// find the base_index we need to add to the local sheet index on sheet s0
			int base_index = 0;
			for(int s = 0; s < s0; ++s)
				base_index += sheets[s].getMaxIndex();

			// Now loop over a search-size sized grid on sheet s0
			for (int i = std::max(0, i0 - searchsize); i < std::min(sheets[s0].getShape(1,0)  - sheets[s0].getShape(0,0), i0 + searchsize); ++i) {
				for (int j = std::max(0,j0 - searchsize); j < std::min(sheets[s0].getShape(1,1) - sheets[s0].getShape(0,1), j0 + searchsize); ++j) {
					for (int l = 0; l < sheets[s0].getNumAtoms(); ++l) {
						// Current grid
						int grid_2[3] = {i,j,l};
						
						// Current index
						int k2 = sheets[s0].gridToIndex(grid_2);
						
						// Check if an orbital exists here (i.e. k = -1 if no orbital exists at that grid position)
						if (k2 != -1){
						
							// Get current position
							double x2 = posAtomIndex(k2+base_index,0);
							double y2 = posAtomIndex(k2+base_index,1);
							
							// If the positions are within the cutoff range we save their indices as a pair for our tight-binding model
							if ((x2 - pos_here[0])*(x2 - pos_here[0]) + (y2 - pos_here[1])*(y2 - pos_here[1]) < inter_cutoff*inter_cutoff) {
								std::vector<int> pair_here;
								pair_here.push_back(kh);
								pair_here.push_back(k2 + base_index);
								pair_array.push_back(pair_here);
							}
						}
					}
				}
			}
		}
		
		// if we are not on the "highest" sheet, we look for pairs from the sheet below
		if (sh < max_sheets - 1) {
		
			// "0" determines the center of our search range on sheet s0
			int i0 = findNearest(pos_here, sh + 1, 0);
			int j0 = findNearest(pos_here, sh + 1, 1);
			int s0 = sh + 1;
			
			// prints how index 0 is paired to the grid above it (used as a troubleshooting print)
			//if (kh == 0)
			//	printf("index %d: [%d,%d,%d] to [%d,%d,%d] \n",kh,ih,jh,sh,i0,j0,s0);
		
			// find the base_index we need to add to the local sheet index on sheet s0
			int base_index = 0;
			for(int s = 0; s < s0; ++s)
				base_index += sheets[s].getMaxIndex();
			
			// Now loop over a search-size sized grid on sheet s0
			for (int i = std::max(0, i0 - searchsize); i < std::min(sheets[s0].getShape(1,0)  - sheets[s0].getShape(0,0), i0 + searchsize); ++i) {
				for (int j = std::max(0,j0 - searchsize); j < std::min(sheets[s0].getShape(1,1) - sheets[s0].getShape(0,1), j0 + searchsize); ++j) {
					for (int l = 0; l < sheets[s0].getNumAtoms(); ++l) {
						// Current grid
						int grid_2[3] = {i,j,l};
												
						// Current index
						int k2 = sheets[s0].gridToIndex(grid_2);
						
						// Check if an orbital exists here (i.e. k = -1 if no orbital exists at that grid position)						
						if (k2 != -1){

							// Get current position
							double x2 = posAtomIndex(k2+base_index,0);
							double y2 = posAtomIndex(k2+base_index,1);
							
							// If the positions are within the cutoff range we save their indices as a pair for our tight-binding model							
							if ((x2 - pos_here[0])*(x2 - pos_here[0]) + (y2 - pos_here[1])*(y2 - pos_here[1]) < inter_cutoff*inter_cutoff) {
								std::vector<int> pair_here;
								pair_here.push_back(kh);
								pair_here.push_back(k2 + base_index);
								pair_array.push_back(pair_here);
							}
						}
					}
				}
			}
		}
	}
 
}

// ----------------------------------------------------------
//   Gets the global index given the grid information input
// [i,j,o,s] (i,j is 2D grid, o is orbital, s is sheet index)
// ----------------------------------------------------------
int Hstruct::gridToIndex(int (&grid_index)[4]) {
	
	int i = grid_index[0];
	int j = grid_index[1];
	int o = grid_index[2];
	int s = grid_index[3];
	
	int sheet_grid[3] = {i,j,o};
	
	// Just make sure we are not out-of-range from sheet index
	if (s < 0 || s >= max_sheets)
		return -1;
	
	// get local index on from [i,j,o] on sheet s
	int temp_index = sheets[s].gridToIndex(sheet_grid);
	
	// if it doesnt exist on the sheet, just return -1
	if (temp_index == -1)
		return -1;
	
	// otherwise find the global index by adding each previous sheet's max index
	for (int x = 0; x < s; ++x) {
		temp_index += sheets[x].getMaxIndex();
	}
	
	return temp_index;

}

void Hstruct::getIntraPairs(std::vector<int> &array_i, std::vector<int> &array_j, std::vector<double> &array_t, int searchsize) {

	int current_index = 0;
	for (int s = 0; s < max_sheets; ++s){
		sheets[s].getIntraPairs(array_i,array_j,array_t,current_index,searchsize);
		current_index += sheets[s].getMaxIndex();
	}

	/*
	int searchsize = 10;
	
	
	for (int kh = 0; kh < max_index; ++kh){
		
		// Manually add the diagonal elements, they are 0 in every model (energy offset is added in locality.cpp, not here)
		array_i.push_back(kh);
		array_j.push_back(kh);
		array_t.push_back(0);
		
		int i0 = index_array[kh][0];
		int j0 = index_array[kh][1];
		int l0 = index_array[kh][2];
		int s0 = index_array[kh][3];
		
		int mat = sheets[s0].getMat();
		
		double pos_here[3];
		pos_here[0] = posAtomIndex(kh,0);
		pos_here[1] = posAtomIndex(kh,1);
		pos_here[2] = posAtomIndex(kh,2);
		
		for (int i = std::max(0, i0 - searchsize); i < std::min(sheets[s0].getShape(1,0)  - sheets[s0].getShape(0,0), i0 + searchsize); ++i) {
			for (int j = std::max(0,j0 - searchsize); j < std::min(sheets[s0].getShape(1,1) - sheets[s0].getShape(0,1), j0 + searchsize); ++j) {
				for (int l = 0; l < sheets[s0].getNumAtoms(); ++l) {
				
					int grid_2[4] = {i,j,l,s0};
					int k2 = gridToIndex(grid_2);

					if (k2 != -1){
						double x2 = posAtomIndex(k2,0);
						double y2 = posAtomIndex(k2,1);
						double z2 = posAtomIndex(k2,2);
						
						int mat2 = sheets[index_array[k2][3]].getMat(); 
						
						double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, mat);
						if (t != 0){
							array_i.push_back(kh);
							array_j.push_back(k2);
							array_t.push_back(t);
						}
					}
				}
			}
		}
	}
	*/
}

// ---------------------------------------------------------------
//      Gives the positions of each orbital in dimension "dim"
// Used for saving the positions of all atoms to file for plotting
// ---------------------------------------------------------------
void Hstruct::getIndexToPos(double* array_in,int dim){
	
	for (int k = 0; k < max_index; ++k){
		array_in[k] = posAtomIndex(k,dim);
	}	

}

std::vector<std::vector<int> > Hstruct::getVacancyList(int center_index, int nShifts){

	std::vector<std::vector<int> > v_list;
	int center_grid[4];
	for (int i = 0; i < 4; ++i){
		center_grid[i] = index_array[center_index][i];
	}
	
	int shift_avg = (nShifts - 1)/2;
	
	for (int i = 0; i < nShifts; ++i){
		for (int j = 0; j < nShifts; ++j){
			
			std::vector<int> temp_v;
			
			int i_shift = i - shift_avg;
			int j_shift = j - shift_avg;
			
			if (i_shift != 0 || j_shift != 0){
				for (int k = 5; k < 11; ++k){
					int temp_grid[4];
					temp_grid[0] = center_grid[0] + i_shift;
					temp_grid[1] = center_grid[1] + j_shift;
					temp_grid[2] = k;
					temp_grid[3] = center_grid[3];
					temp_v.push_back(gridToIndex(temp_grid));
				}
			}				
			v_list.push_back(temp_v);
		}
	}
	
	// Old vacancy creation method
	/*
	for (int j = 0; j < nShifts*nShifts; ++j){
		std::vector<int> temp_v;
		for (int k = 5; k < 8; ++k){
			int temp_grid[4];
			temp_grid[0] = center_grid[0];
			temp_grid[1] = center_grid[1] + 1 + j;
			temp_grid[2] = k;
			temp_grid[3] = center_grid[3];
			temp_v.push_back(gridToIndex(temp_grid));
		}
		temp_v_list.push_back(temp_v);
	}
	*/
	
	return v_list;
	
	
}

std::vector< std::vector<int> > Hstruct::getTargetList(Loc_params opts){

	std::vector<std::vector<int> > t_list;
	
	int solver_type = opts.getInt("solver_type");
	int strain_type = opts.getInt("strain_type");
	
	if (solver_type == 1 || solver_type == 2) {
	
		std::vector<int> temp_list;
	
		int num_target_sheets = opts.getInt("num_target_sheets");
		std::vector<int> target_sheets = opts.getVecInt("target_sheets");
		
		for (int s_index = 0; s_index < num_target_sheets; ++s_index){
			
			int target_sheet = target_sheets[s_index];
			
			int target_x_offset = ( sheets[target_sheet].getShape(1,0) - sheets[target_sheet].getShape(0,0) ) / 2;
			int target_y_offset = ( sheets[target_sheet].getShape(1,1) - sheets[target_sheet].getShape(0,1) ) / 2;
			int center_grid[4] = {target_x_offset,target_y_offset,0,target_sheet};
			int center_index = gridToIndex(center_grid);
			
			int num_orbs = sheets[target_sheet].getNumAtoms();
			
			for (int orb = 0; orb < num_orbs; ++orb){
					temp_list.push_back(center_index + orb);
			}
		}
		
		t_list.push_back(temp_list);
		
	} else if (solver_type == 3) {

		// Selects a 3x3 grid of unit-cells as targets, centered at the origin of the disk for each sheet

		std::vector<int> temp_list;

		int num_target_sheets = opts.getInt("num_target_sheets");
		std::vector<int> target_sheets = opts.getVecInt("target_sheets");

		for (int s_index = 0; s_index < num_target_sheets; ++s_index){
			int target_sheet = target_sheets[s_index];


			int target_x_offset = ( sheets[target_sheet].getShape(1,0) - sheets[target_sheet].getShape(0,0) ) / 2;
			int target_y_offset = ( sheets[target_sheet].getShape(1,1) - sheets[target_sheet].getShape(0,1) ) / 2;

			for (int i = -1; i < 2; ++i){
					for (int j = -1; j < 2; ++j){
							int temp_grid[4] = {target_x_offset + i,target_y_offset + j,0,target_sheet};
							int temp_index = gridToIndex(temp_grid);

							int num_orbs = sheets[target_sheet].getNumAtoms();

							for (int orb = 0; orb < num_orbs; ++orb){
									temp_list.push_back(temp_index + orb);
							}
					}
			}
				
		}
	}
	// for strain jobs we do a grid of targets around the center orbital, controlled by two free parameters given below.
	else if (solver_type == 5) {
		
		// target sampling grid size, makes(2n+1)^2 samples
		
		// int tsg = 2;
		int tsg = 4;
		// target sampling spacing, # of unit cells between each sample
		int tss = 10;
		
		int num_target_sheets = opts.getInt("num_target_sheets");
		std::vector<int> target_sheets = opts.getVecInt("target_sheets");
		
		for (int i = -tsg; i < tsg+1; ++i){
			for (int j = -tsg; j < tsg+1; ++j){
			
				std::vector<int> temp_list;
				
				for (int s_index = 0; s_index < num_target_sheets; ++s_index){
			
					int target_sheet = target_sheets[s_index];
					int num_orbs = sheets[target_sheet].getNumAtoms();
					
					int target_x_offset = ( sheets[target_sheet].getShape(1,0) - sheets[target_sheet].getShape(0,0) ) / 2;
					int target_y_offset = ( sheets[target_sheet].getShape(1,1) - sheets[target_sheet].getShape(0,1) ) / 2;
					int temp_grid[4] = {target_x_offset + i*tss, target_y_offset + j*tss,0,target_sheet};
					int temp_index = gridToIndex(temp_grid);
				
					if (temp_index != -1){
						for (int orb = 0; orb < num_orbs; ++orb){
							temp_list.push_back(temp_index + orb);
						}
					}
					
				}
				
					
				t_list.push_back(temp_list);
			}
		}
	
	}
	
	return t_list;

}

void Hstruct::makeInterFFTFile(int n_x, int n_y, int L_x, int L_y, int length_x,int length_y, std::string fft_file){

	// !!!!!!!! WARNING !!!!!!!!! 
	//This is hard-coded for a "perfect" (no E,B,vacancies, or strain) twisted bilayer, momentum-space is not expected to work well for more general systems, use real-space instead!
	//
	
	// n_i is the real-space discretization in that direction in the first 1x1 grid
	// L_i is the number of 1x1 grids to append to the input before FFT in each direction. i.e. in R: <--- L_x -- | Origin 1x1 grid | -- L_x ---> 
	// length_i is the total length in the i direction in momentum-space to save, i.e. in K: <----- L_i ------>
	// o_1,o_2 are the orbitals of interest for this FFT
	// theta is the relative angle between 1 and 2*L_x
	// fft_file is the desired file name for the output
	
	int num_orb_1 = sheets[0].getNumAtoms();
	int num_orb_2 = sheets[1].getNumAtoms();
	
	double theta = angles[1] - angles[0];
	
	
	std::ofstream fout(fft_file.c_str());
	
	fout << num_orb_1 << " " << num_orb_2 << std::endl;
	fout << "\n";

	for (int o1 = 0; o1 < num_orb_1; ++o1){
		for (int o2 = 0; o2 < num_orb_2; ++o2){
			
			int mat1 = sheets[0].getMat();
			int mat2 = sheets[1].getMat();
			int z1 = heights[0];
			int z2 = heights[1];
			
			double o1_shift_x = sheets[0].getOrbPos(o1,0);
			double o1_shift_y = sheets[0].getOrbPos(o1,1);
			double o2_shift_x = sheets[1].getOrbPos(o2,0);
			double o2_shift_y = sheets[1].getOrbPos(o2,1); 
			
			double x_pos,y_pos;
			int x_size = n_x*(2*L_x)+2;
			int y_size = n_y*(2*L_y)+2;
			int y_size2 = y_size/2+1; // r2c y-direction size
			double* in = new double[x_size*y_size];

			double dx = 1.0/n_x;
			double dy = 1.0/n_y;
			fftw_complex *out;
			out = (fftw_complex*) fftw_malloc(x_size*y_size2*sizeof(fftw_complex));

			fftw_plan p;
			p = fftw_plan_dft_r2c_2d(x_size,y_size,in,out,FFTW_MEASURE);

			for (int i = 0; i < x_size; i++){
				for (int j = 0; j < y_size; j++){
				
					if (i < x_size/2)
						x_pos = -dx*i;
					else
						x_pos = -dx*i + dx*(x_size-1);
						
					if (j < y_size/2)
						y_pos = -dy*j;
					else
						y_pos = -dy*j + dy*(y_size-1);
						
					in[j + i*y_size] = interlayer_term(0 + o1_shift_x, 0 + o1_shift_y, z1, x_pos + o2_shift_x, y_pos + o2_shift_y, z2, o1, o2, 0, theta, mat1, mat2);
					
				}
			}
			
			fftw_execute(p);

			double rescale = 4*M_PI*M_PI*n_x*n_y;
			
			for (int i = 0; i < x_size*y_size2; i++){
				out[i][0] = out[i][0]/rescale;
				out[i][1] = out[i][1]/rescale;
			}

		   //save_fftw_complex(out,x_size,y_size2,L_x*length_x,L_y*length_y,fft_file);

			// Now we save the FFT to file
			// Need defined:

			//		fftw_complex* out 	: output from fftw3 = out
			//		int x_size, y_size 	: 					= x_size, y_size2
			// 		int x_L, y_L		: 					= L_x*length_x, L_y*length_y
			//		string file			: filename of output= fft_file
			
			int x_L = L_x*length_x;
			int y_L = L_y*length_y;
			y_size = y_size2;
			
			// So redefine: y_size = y_size2, file = fft_file
			// and define: x_L = L_x*length, y_L = L_y*length
			
			fout << o1 << " " << o2 << " " << x_L << " " << y_L << " 0" << std::endl;
			
			// save real data to file in plain-text matrix format
			for (int i = 0; i < 2*x_L; i++) {
				for (int j = 0; j < y_L; j++) {
					if (i < x_L)
						fout << out[j+(x_size - x_L+i)*y_size][0] << " ";
					else
						fout << out[j+(i-x_L)*y_size][0] << " ";
					}
					
				fout << std::endl;
			}
			
			fout << "\n";
			
			fout << o1 << " " << o2 << " " << x_L << " " << y_L << " 1" << std::endl;

			// save complex data to file in plain-text matrix format
			for (int i = 0; i < 2*x_L; i++) {
				for (int j = 0; j < y_L; j++) {
					if (i < x_L)
						fout << out[j+(x_size-x_L+i)*y_size][1] << " ";
					else
						fout << out[j+(i-x_L)*y_size][1] << " ";
					}
					
				fout << std::endl;
			}
			
			fout << "\n";

			fftw_destroy_plan(p);
			fftw_free(out);
			
		}
	}

	fout.close();
}







