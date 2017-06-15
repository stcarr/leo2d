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
Sheet::Sheet(std::vector<std::vector<double> > _a, std::vector<int> _types, std::vector<std::vector<double> > _pos, std::vector<int> _min, std::vector<int> _max, int bc) {
    a = _a;
    min_shape = _min;
    max_shape = _max;
	
	boundary_condition = bc;
	
	// "atoms" <--> "orbitals"
	// (i.e. can have more than 1 "atom" at the same position in real-space, but they get a different "atom_type" flag)
    atom_types = _types;
    atom_pos = _pos;
	
	// default to real-space, not momentum
	solver_space = 0;
	
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
	boundary_condition = input.boundary_condition;
	
	// "atoms" <--> "orbitals"
	// (i.e. can have more than 1 "atom" at the same position in real-space, but they get a different "atom_type" flag)
	atom_types = input.atom_types;
	atom_pos = input.atom_pos;
	
	mat = input.mat;
	
	solver_space = input.solver_space;
	strain_type = input.strain_type;
	strain_file = input.strain_file;
	
	// If momentum space, compute reciprocal lattice
	if (solver_space == 1){
		setReciprocal();
	}
	
	// Set indexing
	
	// no strain
	if (strain_type == 0) {
		setIndex();
	}
	
	// strain from a realspace basis of form x,y,z,i,j,o
	if (strain_type == 1) {
		loadIndexRealspace();
	}
	
	// strain from a configuration space basis of form b_x,b_y,o
	if (strain_type == 2) {
		loadIndexConfiguration();
	}
	
	
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
	boundary_condition = orig.boundary_condition;
    atom_types = orig.atom_types;
    atom_pos = orig.atom_pos;
	mat = orig.mat;
	solver_space = orig.solver_space;
	strain_type = orig.strain_type;
	strain_file = orig.strain_file;

	max_index = orig.max_index;
    grid_array = orig.grid_array;
    index_array = orig.index_array;
	pos_array = orig.pos_array;
	
	for (int i = 0; i < 2; ++i)
		for (int j = 0; j < 2; ++j)
			a_inverse[i][j] = orig.a_inverse[i][j];
			
	if (solver_space == 1){
		b = orig.b;
	}
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
//    Sets Indexing via a .dat file of
//  real-space positions and grid coords.
// ---------------------------------------
void Sheet::loadIndexRealspace(){

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
			for(int d = 0; d < depth; ++ d){
				grid_array[i][j][d] = -1;
			}
		}
	}
	
		
	std::string line;
	std::ifstream in_file;
	in_file.open(strain_file.c_str());
	if (in_file.is_open())
	{
		while ( getline(in_file,line) ) {
			
			std::vector<std::string> val;
			std::stringstream ss(line); // Turn the string into a stream.
			std::string tok;

			while(getline(ss, tok, ' ')) {
				if (tok != ""){
					val.push_back(tok);
				}
			}
			
			//for (int i = 0; i < val.size(); i++)
				//printf("val[%d] = %s, length = %d \n", i,val[i].c_str(),val[i].size());
			
			double pos[3];
			
			// real-space position of this orbital
			pos[0] = atof(val[0].c_str());
			pos[1] = atof(val[1].c_str());
			pos[2] = atof(val[2].c_str());
			
			//printf("pos of atom %d is [%lf, %lf, %lf] \n",k,pos[0],pos[1],pos[2]);
			
			// grid position of this orbital
			int i = atoi(val[3].c_str()) - min_shape[0];
			int j = atoi(val[4].c_str()) - min_shape[1];
			int l = atoi(val[5].c_str()) - 1;
			
			// "grid" size should be roughly:
			// max in x dir: 200 + 400*sc_width Angstroms
			// max in y dir: 100 + 200*sc_height Angstroms
			// so we should always pick sc_height ~ 2*sc_width
			
			// these choices should work up to r ~ 1000 Angstroms.
			int sc_width = 2;
			int sc_height = 4;
			
			for (int x_sc = -sc_width; x_sc < sc_width+1; x_sc++) {
				for(int y_sc = -sc_height; y_sc < sc_height+1; y_sc++) {
			
					// Assuming data = [x,y,z,i,j,o]
					// We need to first add supercell offsets:
					// East  (+x): [402.1358, 0, 0, 164,  -1, 0]
					// North (+y): [0, 232.1734, 0, -54, 109, 0]
					
					int i_temp = i + 164*x_sc -  54*y_sc;
					int j_temp = j - 1*x_sc   + 109*y_sc;
					int l_temp = l;
					
					double pos_x = pos[0] + 402.1358*x_sc;
					double pos_y = pos[1] + 232.1734*y_sc;
					double pos_z = pos[2];
			
					// then we relabel the input so that it 
					// matches our unit-cell modelling for monolayer
					
					if (l_temp == 0){
						l_temp = 1;
						
					} else if (l == 1){
						l_temp = 0;
						i_temp = i_temp + 1;
						j_temp = j_temp + 1;
					}
					
					double new_pos[3];
					new_pos[0] = pos_x;
					new_pos[1] = pos_y;
					new_pos[2] = pos_z;
					
					if (i_temp >= 0 && j_temp >= 0  && i_temp < height && j_temp < width) {
					
						// Add it if it is valid to the shape of the sheet
						if(checkShape(new_pos)) {
							//printf("%d at [%lf, %lf, %lf] \n",k,pos[0],pos[1],pos[2]);
							grid_array[i_temp][j_temp][l_temp] = k;
							
							std::vector<double> temp_pos;
							temp_pos.push_back(pos_x);
							temp_pos.push_back(pos_y);
							temp_pos.push_back(pos_z);
							pos_array.push_back(temp_pos);
							
							//printf("%d, %lf, %lf, %lf, %d, %d, %d \n",k,pos_x,pos_y,pos_z,i_temp,j_temp,l_temp);
							
							++k;
						} 
						// Otherwise it gets index -1 to signify no orbital is there
						else {
							grid_array[i_temp][j_temp][l_temp] = -1;
						}
					}
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

void Sheet::loadIndexConfiguration(){

}

// ---------------------------------------
// Returns the position from orbital index
//       dim: (x = 0, y = 1, z = 2)
// 		  just calls posAtomGrid!!
// ---------------------------------------
double Sheet::posAtomIndex(int index, int dim){
    
	if (strain_type == 0) {
		int grid_here[3];
		
		for (int x = 0; x < 3; ++x){
			grid_here[x] = index_array[index][x];
		}
		
		return posAtomGrid(grid_here, dim);
	} else {
		
		return pos_array[index][dim];
	}
}

// ------------------------------------
// Returns the position from grid_index
//      dim: (x = 0, y = 1, z = 2)
// ------------------------------------
double Sheet::posAtomGrid(int (&grid_index)[3],int dim){
    
		if (strain_type != 0){
			return pos_array[ grid_array[ grid_index[0] ][ grid_index[1] ][ grid_index[2] ] ][ dim ];
		
		} else {
    
			int i = grid_index[0] + min_shape[0];
			int j = grid_index[1] + min_shape[1];
			int l = grid_index[2];
			
			double x;
			double y;
			double z;
			
			// Simply multiply the grid by the unit cell vectors 
			// and add the orbital's position in the unit cell
			//printf("a[1] = [%lf, %lf], a[2] = [%lf, %lf] \n",a[0][0],a[0][1], a[1][0], a[1][1]);
			if (solver_space == 0){
				x = i*a[0][0] + j*a[1][0] + atom_pos[l][0];
				y = i*a[0][1] + j*a[1][1] + atom_pos[l][1];
				z = atom_pos[l][2]; // We here assume that a1 and a2 have no z component (always do 2D materials on x-y plane...)

			} else if (solver_space == 1){
				x = i*b[0][0] + j*b[1][0];
				y = i*b[0][1] + j*b[1][1];		
				z = 0;
			}
			
			//printf("x,y,z = [%lf, %lf, %lf] \n",x,y,z);
			
			if (dim == 0)
				return x;
			if (dim == 1)
				return y;
			if (dim == 2)
				return z;
		}
}

// ---------------------------------------------------------
// Checks if position pos is a valid location for an orbital
//          Currently just checks if r < max_shape
// ---------------------------------------------------------
bool Sheet::checkShape(double (&pos)[3]){

	if (boundary_condition == 1){
		return true;
	}
	
    if (pow(pos[0],2) + pow(pos[1],2) < pow(max_shape[0],2)){
        return true;
	} else {
        return false;
	}
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

double Sheet::getOrbPos(int orb, int dim){
	return atom_pos[orb][dim];
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

void Sheet::getIntraPairs(std::vector<int> &array_i, std::vector<int> &array_j, std::vector<double> &array_t, int start_index, int searchsize){

	if (solver_space == 0) { 
		for (int kh = 0; kh < max_index; ++kh){
			
			int i0 = index_array[kh][0];
			int j0 = index_array[kh][1];
			int l0 = index_array[kh][2];
			
			double pos_here[3];
			pos_here[0] = posAtomIndex(kh,0);
			pos_here[1] = posAtomIndex(kh,1);
			pos_here[2] = posAtomIndex(kh,2);
			
			//if (kh == 0)
				//printf("index %d [%d,%d,%d]: pos_here = [%lf,%lf,%lf] \n",kh,i0,j0,l0,pos_here[0],pos_here[1],pos_here[2]);
			
			for (int i = std::max(0, i0 - searchsize); i < std::min(getShape(1,0)  - getShape(0,0), i0 + searchsize); ++i) {
				for (int j = std::max(0,j0 - searchsize); j < std::min(getShape(1,1) - getShape(0,1), j0 + searchsize); ++j) {
					for (int l = 0; l < getNumAtoms(); ++l) {
					
						int grid_2[3] = {i,j,l};
						int k2 = gridToIndex(grid_2);

						if (k2 != -1){
							double x2 = posAtomIndex(k2,0);
							double y2 = posAtomIndex(k2,1);
							double z2 = posAtomIndex(k2,2);
							
							double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
							if (t != 0 || kh == k2){
							
								//if (kh == 0 || k2 == 0)
									//printf("adding term at search [i,j,l] = [%d, %d, %d]: [%d,%d] = %lf \n",i,j,l,kh,k2,t);
									
								array_i.push_back(kh + start_index);
								array_j.push_back(k2 + start_index);
								array_t.push_back(t);
							}
						}
					}
				}
			}
			
			// periodic wrapping:
			// solves for the 8 unique cases of (i,j) being "outside" the direct grid of the simulation
			if (boundary_condition == 1){
			
				if (i0 - searchsize < 0){
				
					for(int i = i0 - searchsize + (getShape(1,0)  - getShape(0,0)); i < getShape(1,0)  - getShape(0,0); ++i){
					
						for (int j = std::max(0,j0 - searchsize); j < std::min(getShape(1,1) - getShape(0,1), j0 + searchsize); ++j) {
							for (int l = 0; l < getNumAtoms(); ++l) {
							
								int grid_2_for_index[3] = {i,j,l};
								int grid_2_for_pos[3]   = {i - (getShape(1,0) - getShape(0,0)), j, l};
								int k2 = gridToIndex(grid_2_for_index);

								if (k2 != -1){
									double x2 = posAtomGrid(grid_2_for_pos,0);
									double y2 = posAtomGrid(grid_2_for_pos,1);
									double z2 = posAtomGrid(grid_2_for_pos,2);
									
									double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
									if (t != 0 || kh == k2){
										array_i.push_back(kh + start_index);
										array_j.push_back(k2 + start_index);
										array_t.push_back(t);
									}
								}
							}
						}
						
						if (j0 - searchsize < 0){
						
							for (int j = j0 - searchsize + (getShape(1,1) - getShape(0,1)); j < getShape(1,1) - getShape(0,1); ++j) {
								for (int l = 0; l < getNumAtoms(); ++l) {
								
									int grid_2_for_index[3] = {i,j,l};
									int grid_2_for_pos[3]   = {i - (getShape(1,0) - getShape(0,0)), j - (getShape(1,1) - getShape(0,1)), l};
									int k2 = gridToIndex(grid_2_for_index);

									if (k2 != -1){
										double x2 = posAtomGrid(grid_2_for_pos,0);
										double y2 = posAtomGrid(grid_2_for_pos,1);
										double z2 = posAtomGrid(grid_2_for_pos,2);
										
										double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
										if (t != 0 || kh == k2){
											array_i.push_back(kh + start_index);
											array_j.push_back(k2 + start_index);
											array_t.push_back(t);
										}
									}
								}
							}
						}
						
						if (j0 + searchsize > getShape(1,1) - getShape(0,1)){
						
							for (int j = 0; j < j0 + searchsize - (getShape(1,1) - getShape(0,1)); ++j) {
								for (int l = 0; l < getNumAtoms(); ++l) {
								
									int grid_2_for_index[3] = {i,j,l};
									int grid_2_for_pos[3]   = {i - (getShape(1,0) - getShape(0,0)), j + (getShape(1,1) - getShape(0,1)), l};
									int k2 = gridToIndex(grid_2_for_index);

									if (k2 != -1){
										double x2 = posAtomGrid(grid_2_for_pos,0);
										double y2 = posAtomGrid(grid_2_for_pos,1);
										double z2 = posAtomGrid(grid_2_for_pos,2);
										
										double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
										if (t != 0 || kh == k2){
											array_i.push_back(kh + start_index);
											array_j.push_back(k2 + start_index);
											array_t.push_back(t);
										}
									}
								}
							}
						}
						
						
					}
				}
					
				
				if (i0 + searchsize > getShape(1,0)  - getShape(0,0)){
				
					for (int i = 0; i < i0 + searchsize - (getShape(1,0)  - getShape(0,0)); ++i){
					
									
						for (int j = std::max(0,j0 - searchsize); j < std::min(getShape(1,1) - getShape(0,1), j0 + searchsize); ++j) {
							for (int l = 0; l < getNumAtoms(); ++l) {
							
								int grid_2_for_index[3] = {i,j,l};
								int grid_2_for_pos[3]   = {i + (getShape(1,0) - getShape(0,0)), j, l};
								int k2 = gridToIndex(grid_2_for_index);

								if (k2 != -1){
									double x2 = posAtomGrid(grid_2_for_pos,0);
									double y2 = posAtomGrid(grid_2_for_pos,1);
									double z2 = posAtomGrid(grid_2_for_pos,2);
									
									double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
									if (t != 0 || kh == k2){
										array_i.push_back(kh + start_index);
										array_j.push_back(k2 + start_index);
										array_t.push_back(t);
									}
								}
							}
						}
						
						if (j0 - searchsize < 0){
						
							for (int j = j0 - searchsize + (getShape(1,1) - getShape(0,1)); j < getShape(1,1) - getShape(0,1); ++j) {
								for (int l = 0; l < getNumAtoms(); ++l) {
								
									int grid_2_for_index[3] = {i,j,l};
									int grid_2_for_pos[3]   = {i + (getShape(1,0) - getShape(0,0)), j - (getShape(1,1) - getShape(0,1)), l};
									int k2 = gridToIndex(grid_2_for_index);

									if (k2 != -1){
										double x2 = posAtomGrid(grid_2_for_pos,0);
										double y2 = posAtomGrid(grid_2_for_pos,1);
										double z2 = posAtomGrid(grid_2_for_pos,2);
										
										double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
										if (t != 0 || kh == k2){
											array_i.push_back(kh + start_index);
											array_j.push_back(k2 + start_index);
											array_t.push_back(t);
										}
									}
								}
							}
						}
						
						if (j0 + searchsize > getShape(1,1) - getShape(0,1)){
						
							for (int j = 0; j < j0 + searchsize - (getShape(1,1) - getShape(0,1)); ++j) {
								for (int l = 0; l < getNumAtoms(); ++l) {
								
									int grid_2_for_index[3] = {i,j,l};
									int grid_2_for_pos[3]   = {i + (getShape(1,0) - getShape(0,0)), j + (getShape(1,1) - getShape(0,1)), l};
									int k2 = gridToIndex(grid_2_for_index);

									if (k2 != -1){
										double x2 = posAtomGrid(grid_2_for_pos,0);
										double y2 = posAtomGrid(grid_2_for_pos,1);
										double z2 = posAtomGrid(grid_2_for_pos,2);
										
										double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
										if (t != 0 || kh == k2){
											array_i.push_back(kh + start_index);
											array_j.push_back(k2 + start_index);
											array_t.push_back(t);
										}
									}
								}
							}
						}

					}
				} // end of i0 + searchsize > size of grid
				
				if (j0 - searchsize < 0){
					for (int i = std::max(0, i0 - searchsize); i < std::min(getShape(1,0)  - getShape(0,0), i0 + searchsize); ++i) {
						for (int j = j0 - searchsize + (getShape(1,1) - getShape(0,1)); j < getShape(1,1) - getShape(0,1); ++j) {
							for (int l = 0; l < getNumAtoms(); ++l) {
							
								int grid_2_for_index[3] = {i,j,l};
								int grid_2_for_pos[3]   = {i, j - (getShape(1,1) - getShape(0,1)), l};
								int k2 = gridToIndex(grid_2_for_index);

								if (k2 != -1){
									double x2 = posAtomGrid(grid_2_for_pos,0);
									double y2 = posAtomGrid(grid_2_for_pos,1);
									double z2 = posAtomGrid(grid_2_for_pos,2);
									
									double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
									if (t != 0 || kh == k2){
										array_i.push_back(kh + start_index);
										array_j.push_back(k2 + start_index);
										array_t.push_back(t);
									}
								}
							}
						}
					}
				} // end of j0 - searchsize < 0
						
				if (j0 + searchsize > getShape(1,1) - getShape(0,1)){
					for (int i = std::max(0, i0 - searchsize); i < std::min(getShape(1,0)  - getShape(0,0), i0 + searchsize); ++i) {
						for (int j = 0; j < j0 + searchsize - (getShape(1,1) - getShape(0,1)); ++j) {
							for (int l = 0; l < getNumAtoms(); ++l) {
							
								int grid_2_for_index[3] = {i,j,l};
								int grid_2_for_pos[3]   = {i, j + (getShape(1,1) - getShape(0,1)), l};
								int k2 = gridToIndex(grid_2_for_index);

								if (k2 != -1){
									double x2 = posAtomGrid(grid_2_for_pos,0);
									double y2 = posAtomGrid(grid_2_for_pos,1);
									double z2 = posAtomGrid(grid_2_for_pos,2);
									
									double t = intralayer_term(pos_here[0], pos_here[1], pos_here[2], x2, y2, z2, l0, l, mat);
									if (t != 0 || kh == k2){
										array_i.push_back(kh + start_index);
										array_j.push_back(k2 + start_index);
										array_t.push_back(t);
									}
								}
							}
						}
					}
				} // end of j0 + searchsize > size of grid
			
			}
		}
	} else if (solver_space == 1){

		// Momentum space intralayer pairing is always block diagonal
		for (int i = 0; i < getShape(1,0)  - getShape(0,0); ++i) {
			for (int j = 0; j < getShape(1,1) - getShape(0,1); ++j) {
				
				
				// In each unitcell, find all valid orbitals and pair them to each other completely
				// For Momentum space, t is not used so just set it = 0;
				std::vector<int> orbs_here;
				
				for (int l = 0; l < getNumAtoms(); ++l) {
				
					int grid_2[3] = {i,j,l};
					int k2 = gridToIndex(grid_2);

					if (k2 != -1){
						orbs_here.push_back(k2);
					}
				}
				
				int n_orbs = (int) orbs_here.size();
				
				for (int o1 = 0; o1 < n_orbs; ++o1){
					for (int o2 = 0; o2 < n_orbs; ++o2){
						array_i.push_back(orbs_here[o1] + start_index);
						array_j.push_back(orbs_here[o2] + start_index);
						array_t.push_back(0.0);
					}
				}
		
			}
		}
	}
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
	a_inverse[0][1] = -a21/det_a;
	a_inverse[1][0] = -a12/det_a;
	a_inverse[1][1] = a11/det_a;
	
	
}

void Sheet::setReciprocal(){

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

}

double Sheet::crossProd(std::vector<double> x, std::vector<double> y, int dim){

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

double Sheet::getInverse(int i, int j){
	return a_inverse[i][j];
}

double Sheet::getReciprocal(int i, int j){
	return b[i][j];
}