/*
 * File:   hstruct.cpp
 * Author: Stephen
 *
 * Created on January 13, 2016, 3:16 PM
 */

#include "strain.h"

#include <string>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <stdexcept>

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

  printf ("entering loadConfigFile ('%s')\n",config_filename.c_str());
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

  printf ("exiting loadConfigFile\n");

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
  disp_out[0] = interp_4point(scaled_res[0],scaled_res[1],v1_x,v2_x,v3_x,v4_x);

  // set up the interpolation in y
  double v1_y = disp_y[sheet][orb][scaled_beg[0]][scaled_beg[1]];
  double v2_y = disp_y[sheet][orb][scaled_end[0]][scaled_beg[1]];
  double v3_y = disp_y[sheet][orb][scaled_beg[0]][scaled_end[1]];
  double v4_y = disp_y[sheet][orb][scaled_end[0]][scaled_end[1]];

  disp_out[1] = interp_4point(scaled_res[0],scaled_res[1],v1_y,v2_y,v3_y,v4_y);

  // set up the interpolation in z
  double v1_z = disp_z[sheet][orb][scaled_beg[0]][scaled_beg[1]];
  double v2_z = disp_z[sheet][orb][scaled_end[0]][scaled_beg[1]];
  double v3_z = disp_z[sheet][orb][scaled_beg[0]][scaled_end[1]];
  double v4_z = disp_z[sheet][orb][scaled_end[0]][scaled_end[1]];
  //printf("vi_z = [%lf,%lf,%lf,%lf]\n",v1_z,v2_z,v3_z,v4_z);

  disp_out[2] = interp_4point(scaled_res[0],scaled_res[1],v1_z,v2_z,v3_z,v4_z);

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
