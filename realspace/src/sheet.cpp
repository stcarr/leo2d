/* 
 * File:   sheet.cpp
 * Author: Stephen
 * 
 * Created on January 4, 2016, 4:32 PM
 */

#include "sheet.h"
#include <math.h>
#include <stdio.h>

// --------------------------------------
// Old constructor, no use of class Sdata
// --------------------------------------
Sheet::Sheet(std::vector<std::vector<double> > _a, std::vector<int> _types, std::vector<std::vector<double> > _pos, std::vector<int> _min, std::vector<int> _max) {
    a = _a;
    min_shape = _min;
    max_shape = _max;
	
	// "atoms" <--> "orbitals"
	// (i.e. can have more than 1 "atom" at the same position in real-space, but they get a different "atom_type" flag)
    atom_types = _types;
    atom_pos = _pos;
	
	// Set indexing
    setIndex();
	
	// Determine inverse of the grid -> position matrix
	// (i.e. get the position -> grid matrix)
	setInverse();
           
}

// ---------------------
// Use this constructor!
// ---------------------
Sheet::Sheet(Sdata input){
	a = input.a;
	min_shape = input.min_shape;
	max_shape = input.max_shape;
	
	// "atoms" <--> "orbitals"
	// (i.e. can have more than 1 "atom" at the same position in real-space, but they get a different "atom_type" flag)
	atom_types = input.atom_types;
	atom_pos = input.atom_pos;
	
	mat = input.mat;
	
	// Set indexing
	setIndex();
	
	// Determine inverse of the grid -> position matrix
	// (i.e. get the position -> grid matrix)
	setInverse();
	
}

// ----------------
// Copy Constructor
// ----------------
Sheet::Sheet(const Sheet& orig) {
    a = orig.a;
    max_shape = orig.max_shape;
    min_shape = orig.min_shape;
    atom_types = orig.atom_types;
    atom_pos = orig.atom_pos;

	max_index = orig.max_index;
    grid_array = orig.grid_array;
    index_array = orig.index_array;
	
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			a_inverse[i][j] = orig.a_inverse[i][j];
}

Sheet::~Sheet() {
}

// -----------------------------------------------------------------
// Generates max_index, "index -> grid" and "grid -> index" matrices
// -----------------------------------------------------------------
void Sheet::setIndex(){
    
	// Determine the size of the grid which will be populated with orbitals
    int height = max_shape[0] - min_shape[0];
    int width = max_shape[1] - min_shape[1];
    int depth = atom_types.size();

	// k represents the index of the orbital
    int k = 0;
 
	// We generate a 3D array which gives the index at grid position (i,j,l)
	// loop over the grid size, (i,j,l), where l is over the orbitals in the unit cell
    grid_array.resize(height);
    for (int i = 0; i < height; ++i) {
        grid_array[i].resize(width);

        for (int j = 0; j < width; ++j){
            grid_array[i][j].resize(depth);
      
            for (int l = 0; l < depth; ++l){
			
				// Find the position of this orbital
                double pos[3];
                int index_here[3] = {i,j,l};
                for (int x =0; x < 3; ++x){
                    pos[x] = posAtomGrid(index_here,x);
                }
				
				// Add it if it is valid to the shape of the sheet
                if(checkShape(pos)) {
                    grid_array[i][j][l] = k;
					++k;
                } 
				// Otherwise it gets index -1 to signify no orbital is there
                else {
                    grid_array[i][j][l] = -1;
                }
            }
        }
    }
    
	// save the max_index
    max_index = k;
	
	// Now we create the inverse matrix, i.e. that returns a grid position for a specific k
    index_array.resize(k);
	
	// it is a 2D array of shape max_index x 3
    for (int y = 0; y < k; ++y){
        index_array[y].resize(3);
    }
    
	// Loop over the grid_array and use it to populate our index_array
    int temp_grid[3] = {0,0,0};
    for (int i = 0; i < height; ++i) {
        for (int j = 0; j < width; ++j){
            for (int l = 0; l < depth; ++l){
                int k_here = grid_array[i][j][l];
                if (k_here != -1){
                    int index_here[3] = {i,j,l};
                    for (int x = 0; x < 3; ++x){
                        index_array[k_here][x] = index_here[x];
                    }
                }
            }
        }
  
    }
}

// ---------------------------------------
// Returns the position from orbital index
//       dim: (x = 0, y = 1, z = 2)
// 		  just calls posAtomGrid!!
// ---------------------------------------
double Sheet::posAtomIndex(int index, int dim){
    int grid_here[3];
    
    for (int x = 0; x < 3; ++x){
        grid_here[x] = index_array[index][x];
    }
    
    return posAtomGrid(grid_here, dim);
}

// ------------------------------------
// Returns the position from grid_index
//      dim: (x = 0, y = 1, z = 2)
// ------------------------------------
double Sheet::posAtomGrid(int (&grid_index)[3],int dim){
    
    
        int i = grid_index[0] + min_shape[0];
        int j = grid_index[1] + min_shape[1];
        int l = grid_index[2];
        
		// Simply multiply the grid by the unit cell vectors 
		// and add the orbital's position in the unit cell
        double x = i*a[0][0] + j*a[0][1] + atom_pos[l][0];
        double y = i*a[1][0] + j*a[1][1] + atom_pos[l][1];
		
		// We assume that a1 and a2 have no z component (always do 2D materials on x-y plane...)
        double z = atom_pos[l][2];
        
        if (dim == 0)
            return x;
        if (dim == 1)
            return y;
        if (dim == 2)
            return z;
}

// ---------------------------------------------------------
// Checks if position pos is a valid location for an orbital
//          Currently just checks if r < max_shape
// ---------------------------------------------------------
bool Sheet::checkShape(double (&pos)[3]){
        
    if (pow(pos[0],2) + pow(pos[1],2) < pow(max_shape[0],2)){
        return true;}
    else{
        return false;}
}

// ------------------------------------------
// Returns the given index's grid information
//   dim: (i = 0, j = 1, orbital = 2)
// ------------------------------------------
int Sheet::indexToGrid(int k,int dim){
    return index_array[k][dim];
}

// --------------------------------------------
// Returns the given grid_index's orbital index
// --------------------------------------------
int Sheet::gridToIndex(int (&grid_index)[3]){
    int i = grid_index[0];
    int j = grid_index[1];
    int l = grid_index[2];
	
	
    if (i >= 0 && i < max_shape[0] - min_shape[0])
		if (j >= 0 && j < max_shape[1] - min_shape[1])
			if (l >= 0 && l < atom_types.size())
				return grid_array[i][j][l];
    
	return -1;
}

int Sheet::getMaxIndex(){
    return max_index;
}

double Sheet::getUnit(int vec, int dim){
    return a[vec][dim];
}

int Sheet::getShape(int type, int dim){

	if (type == 0)
		return min_shape[dim];
		
	if (type == 1)
		return max_shape[dim];
		
	return 0;
}

// ----------------------------------
// Atoms actually refers to "orbitals"
// ----------------------------------
int Sheet::getNumAtoms(){
	return atom_types.size();
}


// ----------------------
// Returns material index
// ----------------------
int Sheet::getMat(){
	return mat;
}

// ------------------------------------------------------------
// Creates the matrix which converts position -> unit cell grid
// ------------------------------------------------------------
void Sheet::setInverse(){
	
	double a11 = a[0][0];
	double a12 = a[0][1];
	double a21 = a[1][0];
	double a22 = a[1][1];
	double det_a = a11*a22 - a12*a21;
	
	a_inverse[0][0] = a22/det_a;
	a_inverse[0][1] = -a12/det_a;
	a_inverse[1][0] = -a21/det_a;
	a_inverse[1][1] = a11/det_a;
	
	
}

double Sheet::getInverse(int i, int j){
	return a_inverse[i][j];
}