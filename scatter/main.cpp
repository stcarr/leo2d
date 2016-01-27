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

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    std::vector<std::vector<double> > unitCell;
    std::vector<double> a1,a2,a3;
    
    double unitCell_in[3][3] = 
    {
        {1,0,0},
        {0,1,0},
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
    
    int num_atoms = 1;
    
    vector<int> types;
    // 6 is carbon (atomic #)
    types.push_back(6);
    
    vector<vector<double> > pos;
    pos.resize(num_atoms);
    for (int i = 0; i < num_atoms; ++i)
        pos[i].resize(3);

    pos[0][0] = 0.0;
    pos[0][1] = 0.0;
    pos[0][2] = 0.0;

    
    std::vector<int> min;
    min.push_back(-25);
    min.push_back(-25);
    min.push_back(-25);
    
    std::vector<int> max;
    max.push_back(25);
    max.push_back(25);
    max.push_back(25);
	/*
    Sheet s1(unitCell,types,pos,min,max);
    Sheet s2(unitCell,types,pos,min,max);
    vector<Sheet> sheets;
    
    sheets.push_back(s1);
    sheets.push_back(s2);
    */
    
	vector<double> heights, angles;
    
	heights.push_back(0);
    heights.push_back(3);
    
	angles.push_back(0);
    angles.push_back(1);
    
    Locality loc();
	loc.setup();
	loc.initMPI(argc, argv);
	loc.constructGeom();
	loc.constructMatrix();
	loc.solveMatrix();
	loc.getRaw();
	#loc.getProcessed();
	loc.plot();
	#loc.save();
	loc.finMPI();
	
	#Hstruct h1(sheets, angles, heights);
    #printf("s1: %d, s2: %d, h1: %d \n", s1.getMaxIndex(),s2.getMaxIndex(),h1.getMaxIndex());
    
}