/* 
 * File:   hstruct.cpp
 * Author: Stephen
 * 
 * Created on January 13, 2016, 3:16 PM
 */

#include "hstruct.h"
#include<stdio.h>
#include <math.h>

Hstruct::Hstruct(std::vector<Sheet> sheets_in,std::vector<double> angles_in,std::vector<double> heights_in) {
    std::vector<double> blank_shift;
    for(int j = 0; j < 3; ++j)
        blank_shift.push_back(0);
    for(int i = 0; i < sheets_in.size(); ++i){
        sheets.push_back(sheets_in[i]);
        heights.push_back(heights_in[i]);
        angles.push_back(angles_in[i]);
        shifts.push_back(blank_shift);
    }
    max_sheets = sheets.size();
    setIndex();
        
    
}

Hstruct::Hstruct(const Hstruct& orig) {
}

Hstruct::~Hstruct() {
}

void Hstruct::setIndex(){
    int k_total = 0;
    
    for(int s = 0; s < max_sheets; ++s){
        k_total += sheets[s].getMaxIndex();
    }
    
    max_index = k_total;
	
	index_array.resize(max_index);
	
	for(int k = 0; k < max_index; ++k){
	
		index_array[k].resize(4);
		int i_here, j_here, l_here, s_here;
		
		int index_counter = 0;
		
		for(int s = 0; s < sheets.size(); ++s) {
			if (k < index_counter + sheets[s].getMaxIndex())
				s_here = s;
			else
				index_counter += sheets[s].getMaxIndex();
		}
		
		i_here = sheets[s_here].indexToGrid(k - index_counter, 0);
		j_here = sheets[s_here].indexToGrid(k - index_counter, 1);
		l_here = sheets[s_here].indexToGrid(k - index_counter, 2);
		
		index_array[k].push_back(i_here);
		index_array[k].push_back(j_here);
		index_array[k].push_back(l_here);
		index_array[k].push_back(s_here);
	
	}
	
}

int Hstruct::getMaxIndex(){
    return max_index;
}

void Hstruct::setShift(int sheet, std::vector<double> b){
    shifts[sheet] = b;
}

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
    double theta = angles[s];

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

int Hstruct::findNearest(double (&pos)[3],int s, int dim){

	double theta = angles[s];
	double x = pos[0];
	double y = pos[1];
	double z = pos[2];
	
	if (dim == 0){
		double x_new = x*cos(theta) + y*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,0) + shifts[s][1]*sheets[s].getUnit(1,0) + shifts[s][2]*sheets[s].getUnit(2,0);
		int i = std::min(std::max(int(floor(x_new)),sheets[s].getShape(0,0)),sheets[s].getShape(1,0)) - sheets[s].getShape(0,0);
		return i;
	}
	
	if (dim == 1){
		double y_new = y*cos(theta) - x*sin(theta) + shifts[s][0]*sheets[s].getUnit(0,1) + shifts[s][1]*sheets[s].getUnit(1,1) + shifts[s][2]*sheets[s].getUnit(2,1);
		int j = std::min(std::max(int(floor(y_new)),sheets[s].getShape(0,1)),sheets[s].getShape(1,1)) - sheets[s].getShape(0,1);
		return j;
	}
	
}

std::vector<std::vector<int> > Hstruct::getIndexArray(){
	return index_array;
}


std::vector<std::vector<int> > Hstruct::getPairs(){

	// !!!!! NEED TO ADD INTRALAYER PAIRS (only interlayer pairs at the moment) !!!!!

	int searchsize = 4;
	std::vector<std::vector<int> > pair_array;
	int pair_count = 0;
	
	for (int kh = 0; kh < max_index; ++kh){
		
		
		int ih = index_array[kh][0];
		int jh = index_array[kh][1];
		int lh = index_array[kh][2];
		int sh = index_array[kh][3];
		
		double pos_here[3];
		pos_here[0] = posAtomIndex(kh,0);
		pos_here[1] = posAtomIndex(kh,1);
		pos_here[2] = posAtomIndex(kh,2);
		
		if (sh > 0) {
		
			int i0 = findNearest(pos_here, sh - 1, 0);
			int j0 = findNearest(pos_here, sh - 1, 1);
			int s0 = sh - 1;
			
			for (int i = std::max(0, i0 - searchsize); i < std::min(sheets[s0].getShape(1,0)  - sheets[s0].getShape(0,0), i0 + searchsize); ++i) {
				for (int j = std::max(0,j0 - searchsize); j < std::min(sheets[s0].getShape(1,1) - sheets[s0].getShape(0,1), j0 + searchsize); ++j) {
					for (int l = 0; l < sheets[s0].getNumAtoms(); ++l) {
						int grid_2[3] = {i,j,l};
						int k2 = sheets[s0].gridToIndex(grid_2);
						if (k2 != -1){
							std::vector<int> pair_here;
							pair_here.push_back(kh);
							pair_here.push_back(k2);
							pair_array.push_back(pair_here);
							++pair_count;
						}
					}
				}
			}
		}
			
		if (sh < max_sheets - 1) {
			int i0 = findNearest(pos_here, sh + 1, 0);
			int j0 = findNearest(pos_here, sh + 1, 1);
			int s0 = sh + 1;
			
			for (int i = std::max(0, i0 - searchsize); i < std::min(sheets[s0].getShape(1,0)  - sheets[s0].getShape(0,0), i0 + searchsize); ++i) {
				for (int j = std::max(0,j0 - searchsize); j < std::min(sheets[s0].getShape(1,1) - sheets[s0].getShape(0,1), j0 + searchsize); ++j) {
					for (int l = 0; l < sheets[s0].getNumAtoms(); ++l) {
						int grid_2[3] = {i,j,l};
						int k2 = sheets[s0].gridToIndex(grid_2);
						if (k2 != -1){
							std::vector<int> pair_here;
							pair_here.push_back(kh);
							pair_here.push_back(k2);
							pair_array.push_back(pair_here);
							++pair_count;
						}
					}
				}
			}
		}
	}

	return pair_array;
 
}