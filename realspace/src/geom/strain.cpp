/*
 * File:   hstruct.cpp
 * Author: Stephen
 *
 * Created on January 13, 2016, 3:16 PM
 */

#include "strain.h"
#include "tools/numbers.h"

#include <string>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <stdexcept>

using namespace numbers;

// ---------------------
// Use this constructor!
// ---------------------
StrainCalc::StrainCalc() {

}

StrainCalc::StrainCalc(const StrainCalc& orig) {
}

StrainCalc::~StrainCalc() {
}

void StrainCalc::loadConfigFile(std::string config_filename){

  // File should look like:
  /*
      | NUM_SHEETS
      |
  S1  | NUM_ORB GRID_X GRID_Y
      |
      | [ORB_1 DISP_X DATA]
      |
      | [ORB_1 DISP_Y DATA]
      |
      | [ORB_1 DISP_Z DATA]
      |
      | [ORB_2 DISP_X DATA]
      |
      ...
      |
  S2  | NUM_ORB GRID_X GRID_Y
      |
      | [ORB_1 DISP_X DATA]
      |
      ...
      |
  */

  //printf ("entering loadConfigFile ('%s')\n",config_filename.c_str());
  std::ifstream fin(config_filename.c_str());

  // First read the number of sheets
  fin >> num_sheets;

  // Allocate memory accordingly
  num_orb.resize(num_sheets);
  grid.resize(num_sheets);
  disp_x.resize(num_sheets);
  disp_y.resize(num_sheets);
  disp_z.resize(num_sheets);
  // Now go through the file, each sheet should have its own entry
  for (int s = 0; s < num_sheets; ++s){

    grid[s].resize(2);

    // get the header, which is a line consisting of three ints:
    fin >> num_orb[s];
    fin >> grid[s][0];
    fin >> grid[s][1];


    disp_x[s].resize(num_orb[s]);
    disp_y[s].resize(num_orb[s]);
    disp_z[s].resize(num_orb[s]);


    // Go through each orbital in that sheet
    for (int o = 0; o < num_orb[s]; ++o){
      disp_x[s][o].resize(grid[s][0]);
      disp_y[s][o].resize(grid[s][0]);
      disp_z[s][o].resize(grid[s][0]);

      // Load the x-displacement info
      for (int i_x = 0; i_x < grid[s][0]; ++i_x){

        disp_x[s][o][i_x].resize(grid[s][1]);

        for (int j_x = 0; j_x < grid[s][1]; ++j_x){

          fin >> disp_x[s][o][i_x][j_x];
        }

      }

      // Load the y-displacement info
      for (int i_y = 0; i_y < grid[s][0]; ++i_y){

        disp_y[s][o][i_y].resize(grid[s][1]);

        for (int j_y = 0; j_y < grid[s][1]; ++j_y){
          fin >> disp_y[s][o][i_y][j_y];
        }

      }

      // Load the z-displacement info
      for (int i_z = 0; i_z < grid[s][0]; ++i_z){

        disp_z[s][o][i_z].resize(grid[s][1]);

        for (int j_z = 0; j_z < grid[s][1]; ++j_z){
          fin >> disp_z[s][o][i_z][j_z];
        }

      }

    }
  }

	fin.close();

  /*
  printf("disp_z = \n");
  for (int y = 0; y < grid[0][1]; ++y){
    printf("        ");
    for (int x = 0; x < grid[0][0]; ++x){
      printf("%lf ",disp_z[0][0][x][y]);
    }
    printf("\n");
  }
  */

  //printf ("exiting loadConfigFile\n");

}

void StrainCalc::setOpts(Job_params opts_in){
	opts = opts_in;
}

std::vector<double> StrainCalc::fourierStrainDisp(std::vector<double> config_in, int sheet, int orb){


    // First we check inputs
    if (config_in.size() < 2){
      throw std::runtime_error("Variable 'config_in' has size less than 2 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (config_in[0] < 0){
      throw std::runtime_error("Variable 'config_in[0]' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (config_in[0] >= 1){
      throw std::runtime_error("Variable 'config_in[0]' >= 1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (config_in[1] < 0){
      throw std::runtime_error("Variable 'config_in[1]' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (config_in[1] >= 1){
      throw std::runtime_error("Variable 'config_in[1]' >= 1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (sheet < 0){
      throw std::runtime_error("Variable 'sheet' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (sheet > num_sheets-1){
      throw std::runtime_error("Variable 'sheet' > num_sheets-1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (orb < 0){
      throw std::runtime_error("Variable 'orb' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }
    if (orb > num_orb[sheet]-1){
      throw std::runtime_error("Variable 'orb' > num_orb-1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
    }

    double sign = 1.0;

    // Rework to force symmetry between layers
    // /*
    if (sheet == 1){
      sheet == 0;
      if (config_in[0] != 0) {
        config_in[0] = 1.0 - config_in[0];
      }
      if (config_in[1] != 0){
          config_in[1] = 1.0 - config_in[1];
      }
      sign = -1.0;
    }
    // */

    // output vector
    std::vector<double> disp_out;
    disp_out.resize(3);

    //double r_scale = 0.1072;
    double r_scale = 0.0524; // value for 26_25 super cell (~1.30 degrees)
    //double r_scale = 0.2;

    double x_c01 =  0.0;
    double x_c10 =  r_scale*sqrt(3.0)/2.0;;
    double x_c11 =  r_scale*sqrt(3.0)/2.0;;

    double y_c01 =  r_scale;
    double y_c10 = -r_scale/2.0;
    double y_c11 =  r_scale/2.0;

    double trig_x = config_in[0]*2.0*PI;
    double trig_y = config_in[1]*2.0*PI;

    disp_out[0] = sign*(x_c10*sin(trig_x) + x_c01*sin(trig_y) + x_c11*sin(trig_x + trig_y) );
    disp_out[1] = sign*(y_c10*sin(trig_x) + y_c01*sin(trig_y) + y_c11*sin(trig_x + trig_y) );
    disp_out[2] = 0.0;
    //printf("disp_out = [%lf,%lf,%lf]\n",disp_out[0],disp_out[1],disp_out[2]);
    return disp_out;

}

std::vector<double> StrainCalc::interpStrainDisp(std::vector<double> config_in, int sheet, int orb){

  // First we check inputs
  if (config_in.size() < 2){
    throw std::runtime_error("Variable 'config_in' has size less than 2 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (config_in[0] < 0){
    throw std::runtime_error("Variable 'config_in[0]' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (config_in[0] >= 1){
    throw std::runtime_error("Variable 'config_in[0]' >= 1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (config_in[1] < 0){
    throw std::runtime_error("Variable 'config_in[1]' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (config_in[1] >= 1){
    throw std::runtime_error("Variable 'config_in[1]' >= 1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (sheet < 0){
    throw std::runtime_error("Variable 'sheet' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (sheet > num_sheets-1){
    throw std::runtime_error("Variable 'sheet' > num_sheets-1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (orb < 0){
    throw std::runtime_error("Variable 'orb' < 0 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }
  if (orb > num_orb[sheet]-1){
    throw std::runtime_error("Variable 'orb' > num_orb-1 on entry to StrainCalc::interpStrainDisp! (in src/geom/strain.cpp) \n");
  }

  int sign = 1;

  // Rework to force symmetry between layers
  // /*
  if (sheet == 1){
    sheet == 0;
    if (config_in[0] != 0) {
      config_in[0] = 1.0 - config_in[0];
    }
    if (config_in[1] != 0){
        config_in[1] = 1.0 - config_in[1];
    }
    sign = -1;
  }
  // */

  // Force symmetry between orbitals
  orb = 0;

  // output vector
  std::vector<double> disp_out;
  disp_out.resize(3);

  // set-up grid information
  std::vector<int> scaled_beg;
  std::vector<int> scaled_end;
  std::vector<double> scaled_res;
  scaled_beg.resize(2);
  scaled_end.resize(2);
  scaled_res.resize(2);

  // figure out where we are on the grid (with residue and boundary conditions loaded properly)
  // for both dimensions (d == 0 -> i, d == 1 -> j)
  for (int d = 0; d < 2; ++d){

    // truncate the input config onto the grid
    scaled_beg[d]  = int(config_in[d]*grid[sheet][d]);

    // compute the ending grid position, with periodic boundary conditions
    if (scaled_beg[d]+1 == grid[sheet][d]){
      scaled_end[d] = 0;
    } else {
      scaled_end[d] = scaled_beg[d]+1;
    }

    // save residue for interpolation
    scaled_res[d] = config_in[d]*grid[sheet][d] - scaled_beg[d];
  }

  // set up the interpolation in x
  double v1_x = disp_x[sheet][orb][scaled_beg[0]][scaled_beg[1]];
  double v2_x = disp_x[sheet][orb][scaled_end[0]][scaled_beg[1]];
  double v3_x = disp_x[sheet][orb][scaled_beg[0]][scaled_end[1]];
  double v4_x = disp_x[sheet][orb][scaled_end[0]][scaled_end[1]];

  //printf("config_in = [%lf,%lf] => [%d,%d] (& [%d,%d]) + [%lf,%lf]\n",config_in[0],config_in[1],scaled_beg[0],scaled_beg[1],scaled_end[0],scaled_end[1],scaled_res[0],scaled_res[1]);
  //printf("vi_x = [%lf,%lf,%lf,%lf]\n",v1_x,v2_x,v3_x,v4_x);
  disp_out[0] = sign*interp_4point(scaled_res[0],scaled_res[1],v1_x,v2_x,v3_x,v4_x);

  // set up the interpolation in y
  double v1_y = disp_y[sheet][orb][scaled_beg[0]][scaled_beg[1]];
  double v2_y = disp_y[sheet][orb][scaled_end[0]][scaled_beg[1]];
  double v3_y = disp_y[sheet][orb][scaled_beg[0]][scaled_end[1]];
  double v4_y = disp_y[sheet][orb][scaled_end[0]][scaled_end[1]];

  disp_out[1] = sign*interp_4point(scaled_res[0],scaled_res[1],v1_y,v2_y,v3_y,v4_y);

  // set up the interpolation in z
  double v1_z = disp_z[sheet][orb][scaled_beg[0]][scaled_beg[1]];
  double v2_z = disp_z[sheet][orb][scaled_end[0]][scaled_beg[1]];
  double v3_z = disp_z[sheet][orb][scaled_beg[0]][scaled_end[1]];
  double v4_z = disp_z[sheet][orb][scaled_end[0]][scaled_end[1]];
  //printf("vi_z = [%lf,%lf,%lf,%lf]\n",v1_z,v2_z,v3_z,v4_z);

  disp_out[2] = sign*interp_4point(scaled_res[0],scaled_res[1],v1_z,v2_z,v3_z,v4_z);

  //printf("disp_out = [%lf,%lf,%lf]\n",disp_out[0],disp_out[1],disp_out[2]);
  return disp_out;

}

double StrainCalc::interp_4point(double x, double y, double v1, double v2, double v3, double v4) {

  // x,y in [0,1] and v1,v2,v3,v4 are the four corners following:
  //
  /*

    1   v3        v4
    |   |
    y   |--(x,y)
    |   |    |
    0   v1---|----v2

        0----x----1


  */
  //

  double value = v1*(1-x)*(1-y) + v2*x*(1-y)+ v3*(1-x)*y + v4*x*y;
  return value;
}

std::vector<double> StrainCalc::supercellDisp(std::vector<double> pos_in, int sheet, int orb){

	std::vector<double> disp;
	disp.resize(3);

  /*
	double amp = 20.0;
	double freq = 1.0;

	// u_x:
	disp[0] = amp*sin(2.0*PI*freq*pos_in[0]);
	// u_y:
	disp[1] = amp*cos(2.0*PI*freq*pos_in[1]);
	// u_z:
	disp[2] = 0.0;
  */

  double sign = 1.0;

  // Rework to force symmetry between layers
  // /*
  if (sheet == 1){
    sheet == 0;
    if (pos_in[0] != 0) {
      pos_in[0] = 1.0 - pos_in[0];
    }
    if (pos_in[1] != 0){
        pos_in[1] = 1.0 - pos_in[1];
    }
    sign = -1.0;
  }
  // */

  //double r_scale = 0.1072;
  double r_scale = 0.0524; // value for 26_25 super cell (~1.30 degrees)
  //double r_scale = 0.2;

  double x_c01 =  0.0;
  double x_c10 =  r_scale*sqrt(3.0)/2.0;;
  double x_c11 =  r_scale*sqrt(3.0)/2.0;;

  double y_c01 =  r_scale;
  double y_c10 = -r_scale/2.0;
  double y_c11 =  r_scale/2.0;

  double trig_x = pos_in[0]*2.0*PI;
  double trig_y = pos_in[1]*2.0*PI;

  disp[0] = sign*(x_c10*sin(trig_x) + x_c01*sin(trig_y) + x_c11*sin(trig_x + trig_y) );
  disp[1] = sign*(y_c10*sin(trig_x) + y_c01*sin(trig_y) + y_c11*sin(trig_x + trig_y) );
  disp[2] = 0.0;

	return disp;

}

std::vector< std::vector<double> > StrainCalc::supercellStrain(std::vector<double> pos_in, int sheet, int orb){

	std::vector< std::vector<double> > strain_here;
	// u_xx, u_xy, u_yx, u_yy;
	strain_here.resize(2);
	strain_here[0].resize(2);
	strain_here[1].resize(2);

  /*
	double amp = 20.0;
	double freq = 1.0;

	// u_xx:
	strain_here[0][0] = -2.0*PI*freq*amp*cos(2.0*PI*freq*pos_in[0]);
	// u_xy:
	strain_here[0][1] =  0.0;
	// u_yx:
	strain_here[1][0] =  0.0;
	// u_yy:
	strain_here[1][1] =  2.0*PI*freq*amp*sin(2.0*PI*freq*pos_in[1]);
  */

  double sign = 1.0;

  // Rework to force symmetry between layers
  // /*
  if (sheet == 1){
    sheet == 0;
    if (pos_in[0] != 0) {
      pos_in[0] = 1.0 - pos_in[0];
    }
    if (pos_in[1] != 0){
        pos_in[1] = 1.0 - pos_in[1];
    }
    sign = -1.0;
  }
  // */

  //double r_scale = 0.1072;
  double r_scale = 0.0524; // value for 26_25 super cell (~1.30 degrees)
  //double r_scale = 0.2;

  double x_c01 =  0.0;
  double x_c10 =  r_scale*sqrt(3.0)/2.0;;
  double x_c11 =  r_scale*sqrt(3.0)/2.0;;

  double y_c01 =  r_scale;
  double y_c10 = -r_scale/2.0;
  double y_c11 =  r_scale/2.0;

  double trig_x = pos_in[0]*2.0*PI;
  double trig_y = pos_in[1]*2.0*PI;

  //disp[0] = sign*(x_c10*sin(trig_x) + x_c01*sin(trig_y) + x_c11*sin(trig_x + trig_y) );
  //disp[1] = sign*(y_c10*sin(trig_x) + y_c01*sin(trig_y) + y_c11*sin(trig_x + trig_y) );
  //disp[2] = 0.0;

  // u_xx:
	strain_here[0][0] = 2.0*PI*sign*(x_c10*cos(trig_x) + x_c11*cos(trig_x + trig_y) );
	// u_xy:
	strain_here[0][1] = 2.0*PI*sign*(x_c01*cos(trig_y) + x_c11*cos(trig_x + trig_y) );
	// u_yx:
	strain_here[1][0] = 2.0*PI*sign*(y_c10*cos(trig_x) + y_c11*cos(trig_x + trig_y) );
	// u_yy:
	strain_here[1][1] = 2.0*PI*sign*(y_c01*cos(trig_y) + y_c11*cos(trig_x + trig_y) );

	// Now we rescale for the supercell:

	std::vector< std::vector<double> > sc = opts.getDoubleMat("supercell");

	std::vector< std::vector<double> > sc_inv;
	sc_inv.resize(2);
	sc_inv[0].resize(2);
	sc_inv[1].resize(2);

	double det = sc[0][0]*sc[1][1] - sc[0][1]*sc[1][0];
	sc_inv[0][0] =  sc[1][1]/det;
	sc_inv[0][1] = -sc[1][0]/det;
	sc_inv[1][0] = -sc[0][1]/det;
	sc_inv[1][1] =  sc[0][0]/det;

	std::vector< std::vector<double> > strain_out;
	strain_out.resize(2);
	strain_out[0].resize(2);
	strain_out[1].resize(2);

	for (int i = 0; i < 2; ++i){
		for (int j = 0; j < 2; ++j){
		strain_out[i][j] = sc_inv[j][0]*strain_here[i][0] + sc_inv[j][1]*strain_here[i][1];
		}
	}

	return strain_out;
}

std::vector<double> StrainCalc::realspaceDisp(std::vector<double> pos_in, int sheet, int orb){

	double lambda = opts.getDouble("strain_lambda");
	double shift = opts.getDouble("strain_shift");
	double a = 2.4768;

	std::vector<double> disp;
	disp.resize(3);


	double x = pos_in[0] - shift;
  double y = pos_in[1];

	double disp_vec[2];

	disp_vec[0] = 0.5*a;
	disp_vec[1] = SQRT3_6*a;

  /*
	if (sheet == 0){
		// u_x:
		disp[0] = disp_vec[0]*atan(x/lambda)/(PI_2);; // = a/PI * atan(x/lambda)
		// u_y:
		disp[1] = disp_vec[1]*atan(x/lambda)/(PI_2);; // = a*SQRT3/(3*PI) * atan(x/lambda)
		// u_z:
		disp[2] = 0.0;
	} else {
		disp[0] = 0.0;
		disp[1] = 0.0;
		disp[2] = 0.0;
	}
  */

  // Make two AB-BA solitons
  int disp_sign = 1;
  if (sheet == 1){
    disp_sign = -1;
  }

  if (y <= 25){
    disp[0] = disp_sign*0.5*disp_vec[0];
    disp[1] = disp_sign*0.5*disp_vec[1];
  } else if (y <= 75){
    double d = (y-25.0)/50.0;
    disp[0] = disp_sign*(0.5-d)*disp_vec[0];
    disp[1] = disp_sign*(0.5-d)*disp_vec[1];
  } else if (y <= 125){
    disp[0] = -disp_sign*0.5*disp_vec[0];
    disp[1] = -disp_sign*0.5*disp_vec[1];
  } else  if (y <= 175){
    double d = (y-125.0)/50.0;
    disp[0] = -disp_sign*(0.5-d)*disp_vec[0];
    disp[1] = -disp_sign*(0.5-d)*disp_vec[1];
  } else {
    disp[0] = disp_sign*0.5*disp_vec[0];
    disp[1] = disp_sign*0.5*disp_vec[1];
  }

	//printf("disp = [%lf, %lf] \n",disp[0],disp[1]);
	return disp;

}

std::vector< std::vector<double> > StrainCalc::realspaceStrain(std::vector<double> pos_in, int sheet, int orb){

	double lambda = opts.getDouble("strain_lambda");
	double shift = opts.getDouble("strain_shift");
	double a = 2.4768;

	std::vector< std::vector<double> > strain_here;
	// u_xx, u_xy, u_yx, u_yy;
	strain_here.resize(2);
	strain_here[0].resize(2);
	strain_here[1].resize(2);

	double x = pos_in[0] - shift;
  double y = pos_in[0] - shift;

	if (sheet == 0){
		// u_xx:
		strain_here[0][0] = (a/(PI*lambda)) / (1+(x/lambda)*(x/lambda));
		// u_xy:
		strain_here[0][1] = (a*SQRT3_6/(PI*lambda)) / (1+(x/lambda)*(x/lambda));
		// u_yx:
		strain_here[1][0] =  strain_here[0][1] ;
		// u_yy:
		strain_here[1][1] =  0.0;

	} else {
		// u_xx:
		strain_here[0][0] =  0.0;
		// u_xy:
		strain_here[0][1] =  0.0;
		// u_yx:
		strain_here[1][0] =  0.0;
		// u_yy:
		strain_here[1][1] =  0.0;
	}

	return strain_here;
}

double StrainCalc::gsfeHeight(std::vector<double> config_in){

  // Direct form by Srolovitz paper
  /*
  double x = config_in[0]  + config_in[1]/2.0;
  double y = config_in[0]*0.0 + config_in[1]*SQRT3_2;

  double psi = y + 1.0/SQRT3;
  double phi = x;

  double pre = 2*PI;
  // GSFE d = 3.39, our d = 3.35
  double c0 =     3.47889 - 0.04;
  double c1 =    -0.02648;
  double c2 =    -0.00352;
  double c3 =     0.00037;
  double c4 =     SQRT3*c1;
  double c5 =    -SQRT3*c3;


  double F =  c0 +
              c1*(cos(pre*(phi + psi/SQRT3)) + cos(pre*(phi - psi/SQRT3)) + cos(pre*2*psi/SQRT3) ) +
              c2*(cos(pre*(phi + psi*SQRT3)) + cos(pre*(phi - psi*SQRT3)) + cos(pre*2*phi) ) +
              c3*(cos(pre*(2*phi + 2*psi/SQRT3)) + cos(pre*(2*phi - 2*psi/SQRT3)) + cos(pre*4*psi/SQRT3) ) +
              c4*(sin(pre*(phi - psi/SQRT3)) - sin(pre*(phi + psi/SQRT3)) + sin(pre*2*psi/SQRT3) ) +
              c5*(sin(pre*(2*phi - 2*psi/SQRT3)) - sin(pre*(2*phi + 2*psi/SQRT3)) + sin(pre*4*psi/SQRT3) );
  return F - 3.35;
  */

  // Better form for config. space
  // GSFE d = 3.39, our d = 3.35
  double c0 =     3.47889 - 0.04;
  double c1 =    -0.02648;
  double c2 =    -0.00352;
  double c3 =     0.00037;
  double c4 =     SQRT3*c1;
  double c5 =    -SQRT3*c3;

  double b0 = c0;
  double b1 = -2.0*c1;
  double b2 = c2;
  double b3 = -2.0*c3;

  double x = 2.0*PI*config_in[0];
  double y = 2.0*PI*config_in[1];

  double F =  b0 +
              b1*(cos(x)          + cos(y)      + cos(x + y)) +
              b2*(cos(x + 2.0*y)  + cos(x - y)  + cos(2.0*x + y)) +
              b3*(cos(2.0*x)      + cos(2.0*y)  + cos(2.0*x + 2.0*y));

  return F - 3.35;


}
