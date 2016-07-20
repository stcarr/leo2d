#include "interlayer_coupling.h"
#include <cmath>
#include <stdio.h>
// assumes sheet 1 to sheet 2 orientation, and sheet 2 is rotated by theta

const double pi2o3 	= M_PI*2/3;
const double pi1o3 	= M_PI/3;
const double pi2	= M_PI/2;
const double pi6 	= M_PI/6;
const double r_cut_graphene = 8;
const double r_cut2_graphene = 7;


// assuming (x,y) aims from sheet 1 to sheet 2 site
// theta rotates sheet counter-clockwise

double interlayer_term(double x1_in, double y1_in, double z1_in, double x2_in, double y2_in, double z2_in, int orbit1, int orbit2, double theta1, double theta2, int mat1, int mat2)
{

	if (mat1 == 0 && mat2 == 0) {
		return inter_graphene(x1_in,y1_in,x2_in,y2_in,orbit1,orbit2,theta1,theta2);
	}
	
	if (mat1 > 0 && mat1 < 5 && mat2 > 0 && mat2 < 5){
		return inter_tmdc(x1_in,y1_in,z1_in,x2_in,y2_in,z2_in,orbit1,orbit2,mat1,mat2);
	}

	return 0;

}

double inter_graphene(double x1_in, double y1_in, double x2_in, double y2_in, int orbit1, int orbit2, double theta1, double theta2) 
{

/*
double x1 = cos(-theta1)*x1_in - sin(-theta1)*y1_in;
double y1 = sin(-theta1)*x1_in + cos(-theta1)*y1_in;

double x2 = cos(-theta1)*x2_in - sin(-theta1)*y2_in;
double y2 = sin(-theta1)*x2_in + cos(-theta1)*y2_in;
*/

double x = x2_in - x1_in;
double y = y2_in - y1_in; 

double theta12 = 0;
double theta21 = 0;

double r = sqrt(x*x+y*y);
double t = 0;

//printf("input for inter term: x = %lf, y = %lf, r  = %lf, orbit1 = %d, orbit 2 = %d, theta1 = %lf, theta2 = %lf \n",x,y,r,orbit1,orbit2,theta1,theta2);

// deal with r = 0 case first (causes problems with theta computation)
if (r == 0.0)
	return 0.3155;

if (r < r_cut_graphene){

	double ac = acos(x/r);
	if ( y < 0 )
		ac = 2*M_PI - ac;

	// theta21 (angle to bond on sheet 1)

	theta21 = ac - theta1;
	if (orbit1 == 1){
		theta21 = theta21 + pi6;
	}
	if (orbit1 == 0){
		theta21 = theta21 - pi6;
	}

	// theta12 (angle to bond on sheet 2)
	
	theta12 = ac - theta2 + M_PI;
	if (orbit2 == 1) {
		theta12 = theta12 + pi6;
	}
	if (orbit2 == 0) {
		theta12 = theta12 - pi6;
	}

	//printf("r = %lf, ac = %lf, orbit2 = %d, theta12 = %lf, orbit1 = %d, theta21 = %lf \n", r, ac, orbit2, theta12, orbit1, theta21);


	double V0 = .3155*exp(-1.7543*(r/2.46)*(r/2.46))*cos(2.001*r/2.46);
	double V3 = -.0688*(r/2.46)*(r/2.46)*exp(-3.4692*(r/2.46 - .5212)*(r/2.46-.5212));
	double V6 = -.0083*exp(-2.8764*(r/2.46-1.5206)*(r/2.46-1.5206))*sin(1.5731*r/2.46);
	
	t = V0+V3*(cos(3*theta12)+cos(3*theta21)) + V6*(cos(6*theta12)+cos(6*theta21));
	if (r > r_cut2_graphene)
	{
		double inside_cut = r_cut2_graphene - r_cut_graphene;
		double inside_cut2 = r - r_cut_graphene;
		t = t*exp(1/(inside_cut*inside_cut)-1/(inside_cut2*inside_cut2));
	}

}

//printf("coupling = %f, from input = %f, %f, %d, %d, %f, %f \n", t, x, y, orbit1, orbit2, theta1, theta2);
//printf("%f inter-term computed. \n",t);

return t;
}

double inter_tmdc(double x1_in, double y1_in, double z1_in, double x2_in, double y2_in, double z2_in, int orbit1, int orbit2, int mat1, int mat2)
{
	double delta_z = z1_in - z2_in;
	
	double nu_sigma;
	double R_sigma;
	double eta_sigma;
	
	double nu_pi;
	double R_pi;
	double eta_pi;
	
	double c;
	double d_xx;
	
	// mat1/2 == 1: MoS_2 bilayer
	// This is for S-S interaction
	if (mat1 == 1 && mat2 == 1){
		nu_sigma 	=  2.627;
		R_sigma 	=  3.128;
		eta_sigma 	=  3.859;
		
		nu_pi 		= -0.708;
		R_pi 		=  2.923;
		eta_pi 		=  5.724;
		
		c 			=  12.29;
		d_xx		=  3.130;
	}
	// mat1/2 == 2: WSe_2 bilayer
	// This is for Se-Se interaction
	else if (mat1 == 2 && mat2 == 2){
		nu_sigma 	=  2.559;
		R_sigma 	=  3.337;
		eta_sigma 	=  4.114;
		
		nu_pi 		= -1.006;
		R_pi 		=  2.927;
		eta_pi 		=  5.185;
		
		c 			=  12.96;
		d_xx		=  3.350;
	}

	double XX_sep = (c/2.0) - d_xx;
	int p_i = 0;
	int p_j = 0;
	
	double r_vec[3] = {x2_in - x1_in, y2_in - y1_in, z2_in - z1_in};
	double r_sq = r_vec[0]*r_vec[0] + r_vec[1]*r_vec[1] + r_vec[2]*r_vec[2];
	double r = sqrt(r_sq);
	
	if (r > 7){
		return 0;
	}
	
	double temp_t = 0;
	
	if (delta_z < 0) {
		delta_z = -delta_z;
	}
	
	double delta = 0.05;
	if ((delta_z > XX_sep - delta) && (delta_z < XX_sep + delta)) {
	
		if (orbit1 == 5 || orbit1 == 8) {
			p_i = 0;
		}
		if (orbit1 == 6 || orbit1 == 9) {
			p_i = 1;
		}
		if (orbit1 == 7 || orbit1 == 10) {
			p_i = 2;
		}
		
		if (orbit2 == 5 || orbit2 == 8) {
			p_j = 0;
		}
		if (orbit2 == 6 || orbit2 == 9) {
			p_j = 1;
		}
		if (orbit2 == 7 || orbit2 == 10) {
			p_j = 2;
		}
		
		double V_sigma = nu_sigma*exp(-pow(r/R_sigma, eta_sigma));
		double V_pi    =    nu_pi*exp(-pow(r/R_pi,    eta_pi   ));
		
		temp_t += (V_sigma - V_pi)*((r_vec[p_i]*r_vec[p_j])/r_sq);
		if (p_i == p_j){
			temp_t += V_pi;
		}
		
		return temp_t;
		
	}

	return 0;
	
}
