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
double alpha = 2.4768;
double alpha2 = alpha*alpha;

double r_sq = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);

double t_arr[8];
t_arr[0] = -2.8922;
t_arr[1] =  0.2425;
t_arr[2] = -0.2656;
t_arr[3] =  0.0235;
t_arr[4] =  0.0524;
t_arr[5] = -0.0209;
t_arr[6] = -0.0148;
t_arr[7] = -0.0211;

double r_sq_arr[8];
r_sq_arr[0] = (1.0/3.0)*alpha2;
r_sq_arr[1] = alpha2;
r_sq_arr[2] = (4.0/3.0)*alpha2;
r_sq_arr[3] = (7.0/3.0)*alpha2;
r_sq_arr[4] = 3.0*alpha2;
r_sq_arr[5] = 4.0*alpha2;
r_sq_arr[6] = (13.0/3.0)*alpha2;
r_sq_arr[7] = (16.0/3.0)*alpha2;

for (int i = 0; i < 8; ++i) {
	if ((r_sq > r_sq_arr[i] - delta) && (r_sq < r_sq_arr[i] + delta))
		t = t_arr[i];
}

// onsite = 0.3504

return t;	
}
