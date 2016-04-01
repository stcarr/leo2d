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
// theta rotates sheet2 counter-clockwise

double inter_graphene(double x, double y, int orbit1, int orbit2, double theta)
{
double theta12 = 0;
double theta21 = 0;

double r = sqrt(x*x+y*y);
double t = 0;

// deal with r = 0 case first (causes problems with theta computation)
if (x == 0.0 && y == 0.0)
	return 0.3155;

if (r < r_cut_graphene)
{
double ac = acos(x/r);

if ((x < 0 && y < 0) || (x > 0 && y < 0))
	ac = 2*M_PI-ac;

if (orbit1 == 1)
	theta21 = ac + pi6;
else
	theta21 = ac - pi6;

while (theta21 >= pi2o3)
	theta21 -= pi2o3;
while (theta21 <= -pi2o3)
	theta21 += pi2o3;

ac = ac - theta;

if (orbit2 == 1)
	theta12 = ac + pi6;
else
	theta12 = ac - pi6;

while (theta12 >= pi2o3)
	theta12 -= pi2o3;
while (theta12 <= -pi2o3)
	theta12 += pi2o3;




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

//printf("interlayer_coupling input = %f, %f, %d, %d, %f \n", x, y, orbit1, orbit2, theta);
//printf("%f inter-term computed. \n",t);

return t;
}
