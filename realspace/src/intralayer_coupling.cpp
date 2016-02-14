/* 
 * File:   intralayer_coupling.cpp
 * Author: Stephen
 * 
 * Created on February 5, 2016, 3:33 PM
 */

 
#include "intralayer_coupling.h"
#include <stdio.h>

void intralayer_terms(int (&grid_index)[3],int mat,std::vector< std::vector<int> > &temp_grid, std::vector<double> &temp_t){
	
	if (mat == 0) // graphene
		intralayer_graphene(grid_index, temp_grid, temp_t);
	
	
	
	// throwaway return if mat is invalid
	
	std::vector<std::vector<double > > temp1;
	std::vector<double> temp2;
	temp2.push_back(0.0);
	temp1.push_back(temp2);
	return temp1;
}

void intralayer_graphene(int (&grid_index)[3],std::vector< std::vector<int> > &temp_grid, std::vector<double> &temp_t){

	int i = grid_index[0];
	int j = grid_index[1];
	int s = grid_index[2];
	double t1 = -2.892;
	double t2 =  0.243;
	double t3 = -0.266;
	double t4 =  0.024;
	double t5 =  0.052;
	double t6 = -0.021;
	double t7 = -0.015;
	double t8 = -0.021;
	double t[8] = {t1,t2,t3,t4,t5,t6,t7,t8};
	
	printf("first term: [%f,%f,s,%f]. \n",i,j,t1);
            
	if (s == 0) {
                int temp_array[39][4] = 
					{	{i  ,j  ,1,1},
						{i-1,j  ,1,1},
                				{i  ,j-1,1,1},
                				{i  ,j+1,0,2},
                				{i+1,j  ,0,2},
                				{i+1,j-1,0,2},
                				{i  ,j-1,0,2},
                				{i-1,j  ,0,2},
                				{i-1,j+1,0,2},
						{i-1,j-1,1,3},
                				{i+1,j-1,1,3},
                				{i-1,j+1,1,3},
                				{i  ,j+1,1,4},
                				{i+1,j  ,1,4},
                				{i+1,j-2,1,4},
                				{i  ,j-2,1,4},
                				{i-2,j  ,1,4},
                				{i-2,j+1,1,4},
                				{i-1,j+2,0,5},
						{i+1,j+1,0,5},
                				{i+2,j-1,0,5},
                				{i+1,j-2,0,5},
                				{i-1,j-1,0,5},
                				{i-2,j+1,0,5},
                				{i  ,j+2,0,6},
                				{i+2,j  ,0,6},
                				{i+2,j-2,0,6},
                				{i  ,j-2,0,6},
                				{i-2,j  ,0,6},
                				{i-2,j+2,0,6},
                				{i-2,j+2,1,7},
                				{i-1,j+2,1,7},
                				{i+2,j-1,1,7},
                				{i+2,j-2,1,7},
                				{i-2,j-1,1,7},
                				{i-1,j-2,1,7},
                				{i+1,j+1,1,8},
                				{i+1,j-3,1,8},
                				{i-3,j+1,1,8},
							};
							
				for (int x = 0; x < 39; ++x)
					std::vector<int> grid_here;
					for (int y = 0; y < 3; ++y)
						grid_here.push_back(temp_array[x][y]);
						temp_grid.push_back(grid_here);
						temp_t.push_back(t[temp_array[x][3] - 1]);
		}
                
            if (s == 1) {
				double temp_array[39][4] = 
					{	{i  ,j  ,0,t1},
						{i+1,j  ,0,t1},
                				{i  ,j+1,0,t1},
						{i  ,j+1,1,t2},
                				{i+1,j  ,1,t2},
                				{i+1,j-1,1,t2},
                				{i  ,j-1,1,t2},
                				{i-1,j  ,1,t2},
                				{i-1,j+1,1,t2},
                				{i+1,j+1,0,t3},
                				{i+1,j-1,0,t3},
                				{i-1,j+1,0,t3},
                				{i  ,j+2,0,t4},
                				{i-1,j+2,0,t4},
                				{i+2,j  ,0,t4},
                				{i+2,j-1,0,t4},
                				{i  ,j-1,0,t4},
                				{i-1,j  ,0,t4},
                				{i-1,j+2,1,t5},
                				{i+1,j+1,1,t5},
                				{i+2,j-1,1,t5},
                				{i+1,j-2,1,t5},
                				{i-1,j-1,1,t5},
                				{i-2,j+1,1,t5},
                				{i  ,j+2,1,t6},
                				{i+2,j  ,1,t6},
                				{i+2,j-2,1,t6},
                				{i  ,j-2,1,t6},
                				{i-2,j  ,1,t6},
                				{i-2,j+2,1,t6},
                				{i+1,j+2,0,t7},
                				{i+2,j+1,0,t7},
                				{i+2,j-2,0,t7},
                				{i+1,j-2,0,t7},
                				{i-2,j+1,0,t7},
                				{i-2,j+2,0,t7},
                				{i-1,j+3,0,t8},
						{i+3,j-1,0,t8},
                				{i-1,j-1,0,t8},
							};
				for (int x = 0; x < 39; ++x)
					for (int y = 0; y < 4; ++y)
						c_array[x][y] = temp_array[x][y];
			}
            
            for (int x = 0; x < 39; ++x){
			
				std::vector<double> temp;
				
				for (int y = 0; y < 4; ++y){
					temp.push_back(c_array[x][y]);
					
				}
				
				outarray.push_back(temp);
				
			}
			printf("c_array[5] = [%f,%f,%f,%f] \n", c_array[5][0], c_array[5][1], c_array[5][2], c_array[5][3]);
			printf("outarray[5] = [%f,%f,%f,%f] \n", outarray[5][0],outarray[5][1],outarray[5][2],outarray[5][3]);	
			return outarray;
			
}
