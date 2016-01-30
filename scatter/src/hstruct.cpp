/* 
 * File:   hstruct.cpp
 * Author: Stephen
 * 
 * Created on January 13, 2016, 3:16 PM
 */

#include "hstruct.h"
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
