/* 
 * File:   intralayer_coupling.cpp
 * Author: Stephen
 * 
 * Created on February 5, 2016, 3:33 PM
 */

 
#include "intralayer_coupling.h"

double intralayer_term(double x1, double y1, double z1, double x2, double y2, double z2, int mat){
	
	if (mat == 0) // graphene
		intralayer_graphene(x1, y1, z1, x2, y2, z2);
		
}

double intralayer_graphene(double x1, double y1, double z1, double x2, double y2, double z2){

double delta = .1;
double t = 0;

double r_sq = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);

if ((r_sq > 2.0449 - delta) && (r_sq < 2.0449 + delta))
    t = -2.8922;
else if ((r_sq > 6.1345 - delta) && (r_sq < 6.1345 + delta))
    t =  0.2425;
else if ((r_sq > 14.3141 - delta) && (r_sq < 14.3141 + delta))
    t = -0.2656;
else if ((r_sq > 18.4041 - delta) && (r_sq < 18.4041 + delta))
    t =  0.0235;

// t5 = 0.0524
// t6 = -0.0209
// t7 = -0.0148
// t8 = -0.0211
// onsite = 0.3504


return t;	
}
