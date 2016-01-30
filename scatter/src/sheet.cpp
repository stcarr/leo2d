/* 
 * File:   sheet.cpp
 * Author: Stephen
 * 
 * Created on January 4, 2016, 4:32 PM
 */

#include "sheet.h"
#include <math.h>

Sheet::Sheet(std::vector<std::vector<double> > _a, std::vector<int> _types, std::vector<std::vector<double> > _pos, std::vector<int> _min, std::vector<int> _max) {
    a = _a;
    min_shape = _min;
    max_shape = _max;
    atom_types = _types;
    atom_pos = _pos;
    setIndex();
           
}

Sheet::Sheet(Sdata input){
	a = input.a;
	min_shape = input.min_shape;
	max_shape = input.max_shape;
	atom_types = input.atom_types;
	atom_pos = input.atom_pos;
	setIndex();
	
}

Sheet::Sheet(const Sheet& orig) {
    a = orig.a;
    max_shape = orig.max_shape;
    min_shape = orig.min_shape;
    max_index = orig.max_index;
    atom_types = orig.atom_types;
    atom_pos = orig.atom_pos;
    grid_array = orig.grid_array;
    index_array = orig.index_array;
}

Sheet::~Sheet() {
}

void Sheet::setIndex(){
    
    int height = max_shape[0] - min_shape[0];
    int width = max_shape[1] - min_shape[1];
    int depth = atom_types.size();
  
    int k = -1;
  
    grid_array.resize(height);
    for (int i = 0; i < height; ++i) {
        grid_array[i].resize(width);

        for (int j = 0; j < width; ++j){
            grid_array[i][j].resize(depth);
      
            for (int l = 0; l < depth; ++l){
                double pos[3];
                int index_here[3] = {i,j,l};
                for (int x =0; x < 3; ++x){
                    pos[x] = posAtomGrid(index_here,x);
                }
                if(checkShape(pos)) { 
                    grid_array[i][j][l] = k;
                    ++k;
                } 
                else {
                    grid_array[i][j][l] = -1;
                }
            }
        }
    }
    
  
    max_index = k;
    index_array.resize(k);
    for (int y = 0; y < k; ++y){
        index_array[y].resize(3);
    }
    
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

double Sheet::posAtomIndex(int index, int dim){
    int grid_here[3];
    
    for (int x = 0; x < 3; ++x){
        grid_here[x] = index_array[index][x];
    }
    
    return posAtomGrid(grid_here, dim);
}
    
double Sheet::posAtomGrid(int (&grid_index)[3],int dim){
    
    
        int i = grid_index[0] + min_shape[0];
        int j = grid_index[1] + min_shape[1];
        int l = grid_index[2];
        
        double x = i*a[0][0] + j*a[1][0] + atom_pos[l][0];
        double y = i*a[0][1] + j*a[1][1] + atom_pos[l][1];
        double z = atom_pos[l][2];
        
        if (dim == 0)
            return x;
        if (dim == 1)
            return y;
        if (dim == 2)
            return z;
}

bool Sheet::checkShape(double (&pos)[3]){
        
    if (pow(pos[0],2) + pow(pos[1],2) < pow(max_shape[0],2)){
        return true;}
    else{
        return false;}
}

int Sheet::indexToGrid(int k,int dim){
    return index_array[k][dim];
}

int Sheet::gridToIndex(int (&grid_index)[3]){
    int i = grid_index[0];
    int j = grid_index[1];
    int l = grid_index[2];
        
    int k = grid_array[i][j][l];
    return k;
}

int Sheet::getMaxIndex(){
    return max_index;
}

double Sheet::getUnit(int vec, int dim){
    return a[vec][dim];
}