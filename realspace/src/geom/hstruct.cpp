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

int Hstruct::indexToSheet(int k){

  if (k < 0 || k > max_index){
    return -1;
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

}

// --------------------------------
// Find the dim position of index k
//  (for dim: x = 0, y = 1, z = 2)
// --------------------------------
double Hstruct::posAtomIndex(int k, int dim){
    if (k < 0 || k > max_index){
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

	// If we are in Momentum-space, we need sheet 0 to have the geometry of sheet 1, so we swap the angles
	if (solver_space == 0) {
		theta = angles[s];
	} else if (solver_space == 1) {
		if (s == 0) {
			theta = angles[1];
		} else if (s == 1) {
			theta = angles[0];
		}
	}

    if (dim == 0){
        double x = local_x*cos(theta) - local_y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0);
        return x;
    }
    if (dim == 1){
        double y = local_x*sin(theta) + local_y*cos(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1);
        return y;
    }
    if (dim == 2){
        double z = local_z + heights[s];
		return z;
    }
}

// --------------------------------------------
// Find the unit cell on sheet s closest to pos
//           for dim: (i = 0, j = 1)
// --------------------------------------------
int Hstruct::findNearest(double (&pos)[3],int s, int dim){

	double theta = 0;

	if (solver_space == 0) {
		theta = -angles[s];
	} else if (solver_space == 1) {
		if (s == 0) {
			theta = -angles[1];
		} else if (s == 1) {
			theta = -angles[0];
		}
	}

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
		double x_new = x*cos(theta) - y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0);
		double y_new = y*cos(theta) + x*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1);
		double i_new = a_inv[0][0]*x_new + a_inv[0][1]*y_new;
		int i = std::min(std::max(int(floor(i_new)),sheets[s].getShape(0,0)),sheets[s].getShape(1,0)) - sheets[s].getShape(0,0);
		return i;
	}

	if (dim == 1){
		double x_new = x*cos(theta) - y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0);
		double y_new = y*cos(theta) + x*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1);
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
void Hstruct::getInterPairs(std::vector<std::vector<int> > &pair_array, std::vector<std::vector<double> > &supercell_vecs, Job_params opts){

  // We search over a searchsize x searchsize sized grid of unitcells
  int searchsize = Materials::inter_search_radius(sheets[0].getMat());
  // We do not save pairs that are farther apart than this (in Angstroms)
  double inter_cutoff = 12.5;

  int boundary_condition = opts.getInt("boundary_condition");
  std::vector< std::vector<double> > supercell;
  if (boundary_condition == 1){
    supercell = opts.getDoubleMat("supercell");
  }


	// loop over all orbitals (kh = "k here")
	for (int kh = 0; kh < max_index; ++kh){

    std::vector<std::vector<int> > kh_pair_array;
    std::vector<std::vector<double> > kh_supercell_vecs;

		//printf("kh = %d \n",kh);

		// Get the current grid information (ih = i "here")
		int ih = index_array[kh][0];
		int jh = index_array[kh][1];
		int lh = index_array[kh][2];
		int sh = index_array[kh][3];

		// And the position information
		double kh_pos_here[3];
		kh_pos_here[0] = posAtomIndex(kh,0);
		kh_pos_here[1] = posAtomIndex(kh,1);
		kh_pos_here[2] = posAtomIndex(kh,2);
		// If in momentum space, take reciprocal, i.e. K couples to -K (NOT K couples to K)
		if (solver_space == 1){
			kh_pos_here[0] = -kh_pos_here[0];
			kh_pos_here[1] = -kh_pos_here[1];
			kh_pos_here[2] = -kh_pos_here[2];
		}


    int sc_x_stride = 0;
    int sc_y_stride = 0;

    if (boundary_condition == 1){
      sc_x_stride = 1;
      sc_y_stride = 1;
    }

    // We check all nearby supercells
    for (int dx = -sc_x_stride; dx < sc_x_stride+1; ++dx){
      for (int dy = -sc_y_stride; dy < sc_y_stride+1; ++dy){

        double pos_here[3];

        // compute the supercell vector (for k sampling usually)
        std::vector<double> sc_vec;
        sc_vec.resize(2);
        sc_vec[0] = 0;
        sc_vec[1] = 0;

        if (boundary_condition == 1){
          sc_vec[0] = dx*supercell[0][0] + dy*supercell[1][0];
          sc_vec[1] = dx*supercell[0][1] + dy*supercell[1][1];
        }

        pos_here[0] = kh_pos_here[0] + sc_vec[0];
        pos_here[1] = kh_pos_here[1] + sc_vec[1];
        pos_here[2] = kh_pos_here[2];

    		// if we are not on the "lowest" sheet, we look for pairs from the sheet above
    		if (sh > 0) {

    			// "0" determines the center of our search range on sheet s0
    			int i0 = findNearest(pos_here, sh - 1, 0);
    			int j0 = findNearest(pos_here, sh - 1, 1);
    			int s0 = sh - 1;


    			/*
    			// For debugging/checking findNearest()
    			int t_grid[4];
    			t_grid[0] = i0;
    			t_grid[1] = j0;
    			t_grid[2] = 0;
    			t_grid[3] = s0;
    			int c_k = gridToIndex(t_grid);

    			printf("[%lf, %lf] -> [%lf, %lf] \n",pos_here[0],pos_here[1],posAtomIndex(c_k,0),posAtomIndex(c_k,1));
    			*/

    			//printf("nearest = [%d, %d, %d] \n", i0, j0, s0);


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
    						/*
    						if (i == i0 && j == j0){
    							printf("center k2 = %d \n", k2+base_index);
    						}
    						if (kh == 90)
    							printf("k2 = %d \n",k2+base_index);
    						*/


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
    								kh_pair_array.push_back(pair_here);
                    kh_supercell_vecs.push_back(sc_vec);
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

    			//printf("nearest = [%d, %d, %d] \n", i0, j0, s0);

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
    						/*
    						if (i == i0 && j == j0){
    							printf("center k2 = %d \n", k2+base_index);
    						}
    						*/
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
    								kh_pair_array.push_back(pair_here);
                    kh_supercell_vecs.push_back(sc_vec);

    							}
    						}
    					}
    				}
    			}
    		}


      } // end of dy loop
    } // end of dx loop

    // Now we reorder kh_pair_array and kh_supercell_vecs
    if (boundary_condition == 1){
      orderPairs(kh_pair_array,kh_supercell_vecs);
    }

    for (int i = 0; i < kh_pair_array.size(); ++i){
      pair_array.push_back(kh_pair_array[i]);
      supercell_vecs.push_back(kh_supercell_vecs[i]);
    }

	} // end of kh loop

}

void Hstruct::orderPairs(std::vector< std::vector<int> >& pairs, std::vector< std::vector<double> >& sc_vecs) {

  int num_pairs = (int) pairs.size();

  // a bubblesort should be good enough for this (can test timing later)

  for (int i = num_pairs-1; i >= 0; --i){
    for (int j = 1; j < i+1; ++j){
      // check if the (n-1)'th pair index is larger than the n'th index
      if (pairs[j-1][1] > pairs[j][1]){

        // if so, swap them

        std::vector<int> temp_pair;
        temp_pair.resize(2);
        std::vector<double> temp_sc;
        temp_sc.resize(2);

        temp_pair = pairs[j-1];
        temp_sc = sc_vecs[j-1];

        pairs[j-1] = pairs[j];
        sc_vecs[j-1] = sc_vecs[j];

        pairs[j] = temp_pair;
        sc_vecs[j] = temp_sc;

      }

    }
  }

}

// -----------------------------------------
//     Allocates the Shift configuration
//        for all atoms in a bilayer
// -----------------------------------------
void Hstruct::getShiftConfigs(std::vector<std::vector<double> > &config_array, Job_params opts){

  if (max_sheets == 2){

    std::vector< std::vector< std::vector<double> > > unit_cells;
    std::vector< std::vector< std::vector<double> > > inv_unit_cells;

    unit_cells.resize(2);
    inv_unit_cells.resize(2);
    for (int s = 0; s < 2; ++s){

      double theta = angles[s];

      unit_cells[s].resize(2);
      inv_unit_cells[s].resize(2);
      for (int x = 0; x < 2; ++x){
        unit_cells[s][x].resize(2);
        inv_unit_cells[s][x].resize(2);
      }

      unit_cells[s][0][0] = cos(theta)*sheets[s].getUnit(0,0) - sin(theta)*sheets[s].getUnit(0,1);
      unit_cells[s][0][1] = sin(theta)*sheets[s].getUnit(0,0) + cos(theta)*sheets[s].getUnit(0,1);
      unit_cells[s][1][0] = cos(theta)*sheets[s].getUnit(1,0) - sin(theta)*sheets[s].getUnit(1,1);
      unit_cells[s][1][1] = sin(theta)*sheets[s].getUnit(1,0) + cos(theta)*sheets[s].getUnit(1,1);

      double det = unit_cells[s][0][0]*unit_cells[s][1][1] - unit_cells[s][0][1]*unit_cells[s][1][0];

    	inv_unit_cells[s][0][0] = unit_cells[s][1][1]/det;
    	inv_unit_cells[s][0][1] = -unit_cells[s][1][0]/det;
    	inv_unit_cells[s][1][0] = -unit_cells[s][0][1]/det;
    	inv_unit_cells[s][1][1] = unit_cells[s][0][0]/det;

    }

    config_array.resize(max_index);
    for (int k = 0; k < max_index; ++k){
      config_array[k].resize(2);
      double pos_here[3];
      for (int d = 0; d < 3; ++d){
        pos_here[0] = posAtomIndex(k,d);
      }
      int s_here = indexToSheet(k);
      int s = -1;
      if (s_here == 0){
        s = 1;
      } else if (s_here == 1){
        s = 0;
      }

      int new_i = findNearest(pos_here,s,0);
      int new_j = findNearest(pos_here,s,1);

      int found = 0;
      // now we try to find where atom k is in relation to the other sheet
      for (int i = new_i-2; i < new_i + 3; ++i){
        for (int j = new_j-2; j < new_j + 3; ++j){

          if (found == 0){

            int new_grid[3];
            new_grid[0] = i;
            new_grid[1] = j;
            new_grid[2] = 0;

            double new_pos[2];
            new_pos[0] = sheets[s].posAtomGrid(new_grid,0);
            new_pos[1] = sheets[s].posAtomGrid(new_grid,1);

            double config_real[2];
            config_real[0] = pos_here[0] - new_pos[0];
            config_real[1] = pos_here[1] - new_pos[1];

            double config[2];
            config[0] = config_real[0]*inv_unit_cells[s][0][0] + config_real[1]*inv_unit_cells[s][0][1];
            config[1] = config_real[0]*inv_unit_cells[s][1][0] + config_real[1]*inv_unit_cells[s][1][1];

            if (config[0] >= 0 && config[0] < 1 && config[1] >= 0 && config[1] < 1){
              config_array[k][0] = config[0];
              config_array[k][1] = config[1];

              found = 1;
            }

          }

        }
      }

      if(found == 0){
        throw std::runtime_error("Hstruct::getShiftConfigs failed to find an atoms shift configuration!!");
      }
    }

  } else {
    throw std::runtime_error("Hstruct::getShiftConfigs only implemented for bilayer!!\n");
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

void Hstruct::getIntraPairs(std::vector<int> &array_i, std::vector<int> &array_j, std::vector<double> &array_t, std::vector< std::vector<double> > &sc_vecs, Job_params opts) {

	int current_index = 0;
	for (int s = 0; s < max_sheets; ++s){
		sheets[s].getIntraPairs(array_i,array_j,array_t,sc_vecs,opts,current_index);
		current_index += sheets[s].getMaxIndex();
	}

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

std::vector<std::vector<double> > Hstruct::getSupercellVecs(){



}

double Hstruct::getUnitArea(int s){

	double a11 = sheets[s].getUnit(0,0);
	double a12 = sheets[s].getUnit(1,0);
	double a21 = sheets[s].getUnit(0,1);
	double a22 = sheets[s].getUnit(1,1);

	double area = (a11*a22 - a12*a21);
	if (area > 0){
		return area;
	} else {
		return -area;
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

std::vector< std::vector<int> > Hstruct::getTargetList(Job_params opts){

	std::vector<std::vector<int> > t_list;

	int solver_type = opts.getInt("solver_type");
	int strain_type = opts.getInt("strain_type");

	if (solver_type == 1 || solver_type == 2) {

		std::vector<int> temp_list;

		int num_target_sheets = opts.getInt("num_target_sheets");
		std::vector<int> target_sheets = opts.getIntVec("target_sheets");

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
		std::vector<int> target_sheets = opts.getIntVec("target_sheets");

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
		std::vector<int> target_sheets = opts.getIntVec("target_sheets");

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

				if (!temp_list.empty())
					t_list.push_back(temp_list);
			}
		}

	}

	return t_list;

}

void Hstruct::makeInterFFTFile(int n_x, int n_y, int L_x, int L_y, int length_x, int length_y, double A1, double A2, std::string fft_file){

	// !!!!!!!! WARNING !!!!!!!!!
	//This is hard-coded for a "perfect" (no E,B,vacancies, or strain) twisted bilayer, momentum-space is not expected to work well for more general systems, use real-space instead!
	//

	// n_i is the real-space discretization in that direction in the first 1x1 grid
	// L_i is the number of 1x1 grids to append to the input before FFT in each direction. i.e. in R: <--- L_x -- | Origin 1x1 grid | -- L_x --->
	// length_i is the total length in the i direction in momentum-space to save, i.e. in K: <----- L_i ------>
	// o_1,o_2 are the orbitals of interest for this FFT
	// theta is the relative angle between 1 and 2*L_x
	// A1, A2 are real-space areas of each layer's unit cell
	// fft_file is the desired file name for the output

	int num_orb_1 = sheets[0].getNumAtoms();
	int num_orb_2 = sheets[1].getNumAtoms();

  Materials::Mat mat1 = sheets[0].getMat();
  Materials::Mat mat2 = sheets[1].getMat();

	double angle1 = angles[0];
  double angle2 = angles[1];


	std::ofstream fout(fft_file.c_str());

	fout << num_orb_1 << " " << num_orb_2 << std::endl;
	fout << "\n";

	for (int o1 = 0; o1 < num_orb_1; ++o1){
		for (int o2 = 0; o2 < num_orb_2; ++o2){

			double z1 = heights[0];
			double z2 = heights[1];

			double o1_shift_x = sheets[0].getOrbPos(o1,0);
			double o1_shift_y = sheets[0].getOrbPos(o1,1);
			double o1_shift_z = sheets[0].getOrbPos(o1,2);
			double o2_shift_x = sheets[1].getOrbPos(o2,0);
			double o2_shift_y = sheets[1].getOrbPos(o2,1);
			double o2_shift_z = sheets[1].getOrbPos(o2,2);

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

					//printf("[%lf, %lf, %lf, %lf, %lf, %lf, %d, %d, 0, %lf, %d, %d] \n", o1_shift_x, o1_shift_y, z1, x_pos+o2_shift_x, y_pos+o2_shift_y, z2, o1, o2, theta, mat1, mat2);
          std::array<double, 3> disp = {{ x_pos, y_pos, z2-z1 }};
					in[j + i*y_size] = Materials::interlayer_term(o1, o2, disp, angle1, angle2, mat1, mat2)/(sqrt(A1*A2));

				}
			}

			// !!! START DEBUG !!!
			// Debug on the "in" matrix
			/*

			// consider [0,1] off-diagonal term (was causing issues).
			if (o1 == 0 && o2 == 1){
				std::ofstream fout_debug("interlayer_input_1_to_2.dat");
				for (int i = 0; i < x_size; i++){
					for (int j = 0; j < y_size; j++){
						fout_debug << in[j + i*y_size] << " ";
					}
					fout_debug << std::endl;
				}
				fout_debug.close();
			}

			*/
			// !!! END DEBUG !!!

			fftw_execute(p);


			double rescale = n_x*n_y;

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

			delete in;

		}
	}

	fout.close();
}