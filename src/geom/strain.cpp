/*
 * File:   hstruct.cpp
 * Author: Stephen
 *
 * Created on January 13, 2016, 3:16 PM
 */

#include "strain.h"
#include "tools/numbers.h"

#include <cstring> // for strtok
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


// Loads four text files (fourier relations on SUPERCELL data, from k-dot-p codebase).
// These need to be converted into configuration space OR the moire pattern needs to be figured out
void StrainCalc::loadFourierConfigFile_interp(std::string thetas_filename, std::string x_filename, std::string y_filename, std::string z_filename){

  // load thetas.txt
  std::ifstream fin;
  fin.open(thetas_filename.c_str(), std::ifstream::in);

  double temp_var;
  while(!fin.eof()) {
    fin >> temp_var;
    theta_list.push_back(temp_var);
  }
  fin.close();

  // load coeffs_x.txt
  fin.open(x_filename.c_str(), std::ifstream::in);

	std::string line = "";
	while (getline(fin, line))
	{
    std::vector<double> vec;

    char * str = (char*)line.c_str();
    char * pch;
    pch = strtok(str," ,"); // delimiters for token
    while (pch != NULL)
    {
      vec.push_back(std::strtod(pch,NULL));
      pch = strtok (NULL, " ,"); // move location of token
    }

		coeffs_x.push_back(vec);
    //printf("first values of coeffs_x line: [[%lf, %lf, %lf] \n", vec[0], vec[1], vec[2]);

	}
	// Close the File
	fin.close();

  // load coeffs_y.txt
  fin.open(y_filename.c_str(), std::ifstream::in);

  while (getline(fin, line))
  {
    std::vector<double> vec;

    char * str = (char*)line.c_str();
    char * pch;
    pch = strtok(str," ,"); // delimiters for token
    while (pch != NULL)
    {
      vec.push_back(std::strtod(pch,NULL));
      pch = strtok (NULL, " ,"); // move location of token
    }

    coeffs_y.push_back(vec);
    //printf("first values of coeffs_y line: [[%lf, %lf, %lf] \n", vec[0], vec[1], vec[2]);

  }
  // Close the File
  fin.close();


  // load coeffs_z.txt
  fin.open(z_filename.c_str(), std::ifstream::in);

  while (getline(fin, line))
  {
    std::vector<double> vec;

    char * str = (char*)line.c_str();
    char * pch;
    pch = strtok(str," ,"); // delimiters for token
    while (pch != NULL)
    {
      vec.push_back(std::strtod(pch,NULL));
      pch = strtok (NULL, " ,"); // move location of token
    }

    coeffs_z.push_back(vec);
    //printf("first values of coeffs_x line: [[%lf, %lf, %lf] \n", vec[0], vec[1], vec[2]);

  }
  // Close the File
  fin.close();


  // MATLAB code, for later implementation
  /*
  d_theta = theta_list(2) - theta_list(1);
  %theta_list(tar_theta);
  [interp_diff, nearest_theta] = min(abs(tar_theta - theta_list));
  if (tar_theta - theta_list(nearest_theta) < 0 )
    nearest_theta = nearest_theta - 1;
    interp_diff = -interp_diff+d_theta;
  end
  interp_diff = interp_diff/d_theta;
  if (tar_theta - theta_list(nearest_theta) < 0 )
    nearest_theta = nearest_theta - 1;
    interp_diff = -interp_diff+d_theta;
  end
  interp_diff = interp_diff/d_theta;
  fprintf("theta = %f \n",theta_list(nearest_theta))
  coeffs_1 = [coeffs_x(:,nearest_theta) coeffs_y(:,nearest_theta) coeffs_z(:,nearest_theta)];
  coeffs_2 = [coeffs_x(:,nearest_theta+1) coeffs_y(:,nearest_theta+1) coeffs_z(:,nearest_theta+1)];
  coeffs = coeffs_1*(1-interp_diff) + coeffs_2*interp_diff;
  */

}

void StrainCalc::loadFourierConfigFile(std::string config_filename){

  // File should look like:
  /*
      \ MAX_K
      \ U_X COEFFS
      \ U_Y COEFFS
      \ U_Z COEFFS
  */

  //printf ("entering loadConfigFile ('%s')\n",config_filename.c_str());
  std::ifstream fin(config_filename.c_str());

  // First read the number of sheets
  fin >> max_k;
  int num_coeffs = max_k*(max_k+1)/2 + 1;

  // Allocate memory accordingly
  coeffs.resize(3);
  for (int i = 0; i < 3; ++i){
    coeffs[i].resize(num_coeffs);
    for (int j = 0; j < num_coeffs; ++j){
      fin >> coeffs[i][j];
    }
  }

}

void StrainCalc::setOpts(Job_params opts_in){
	opts = opts_in;
}

std::vector<double> StrainCalc::fourierStrainDisp_sc(double* r, double* b1, double* b2, int s){

  std::vector<double> temp_disp = fourierStrainDisp(r,b1,b2);

  std::vector<double> disp_out;
  disp_out.resize(3);

  // strain_shift is single-layer strain, so multiply by 2 for this form: h( b + 2u(b) )
  double disp_sign = 1.0; // sign is 1 if s == 0
  if (s == 1){
    disp_sign = -1.0; // change sign if s == 1
  }

  // un-rotate the relaxation and include the sign
  disp_out[0] =  disp_sign*temp_disp[1];
  disp_out[1] = -disp_sign*temp_disp[0];
  disp_out[2] =  disp_sign*temp_disp[2];

  return disp_out;

}

// returns fourier strain in config space (for realspace non-supercell methods), not complete!
std::vector<double> StrainCalc::fourierStrainDisp_config(std::vector<double> r, double* b1, double* b2, int sheet){

    // output vector
    std::vector<double> disp_out;
    disp_out.resize(3);

    double disp_x = 0.0;
    double disp_y = 0.0;
    double disp_z = 0.0;

    double k[2];
    double c[3];
    double k_base[2];
    double c_base[3];

    int idx = 0;

    // loop over [i,j] Fourier components
    for (int i = 0; i < max_k+1; ++i){
        for (int j = 0; j < max(i,1); ++j){
            //printf("[%d, %d] \n",i,j);
            // fixing factor, to account for exp(ikr) vs const term
            double ff = 2.0;
            int rot_max = 3;
            if (i == 0) {
               rot_max = 1;
               ff = 1.0;
            }

            // define Fourier component direction
            k_base[0] =  i*b1[0] + j*b2[0];
            k_base[1] =  i*b1[1] + j*b2[1];

            // define Fourier coefficients for this direction
            c_base[0] =  coeffs[0][idx]; // coeffs are in SUPERCELL basis
            c_base[1] =  coeffs[1][idx]; // later we put them in CONFIG basis
            c_base[2] =  coeffs[2][idx];

            for (int r_idx = 0; r_idx < rot_max; ++r_idx){
                double r_theta = r_idx*PI/3.0;
                k[0] = cos(r_theta)*k_base[0] - sin(r_theta)*k_base[1];
                k[1] = sin(r_theta)*k_base[0] + cos(r_theta)*k_base[1];

                //c = rot^(r_idx)*[coeffs[idx,]];
                c[0] = cos(r_theta)*c_base[0] - sin(r_theta)*c_base[1];
                c[1] = sin(r_theta)*c_base[0] + cos(r_theta)*c_base[1];
                c[2] = c_base[2];

                // the coeffs are for exp(i*k*r), so we mulitply by 2 for sin/cos
                disp_x = disp_x + ff*c[0]*sin(k[0]*r[0] + k[1]*r[1]);
                disp_y = disp_y + ff*c[1]*sin(k[0]*r[0] + k[1]*r[1]);
                disp_z = disp_z + ff*c[2]*cos(k[0]*r[0] + k[1]*r[1]);
            }

            idx++;

        }
    }


    double basis_rot = atan2(b1[1],b1[0]); // above data assumes b1 is along x axis
    double rot_disp_x = cos(basis_rot)*disp_x - sin(basis_rot)*disp_y;
    double rot_disp_y = sin(basis_rot)*disp_x + cos(basis_rot)*disp_y;
    // rotate by ~pi/2 to put into CONFIG basis
    disp_out[0] = -rot_disp_y;
    disp_out[1] = rot_disp_x;
    disp_out[2] = disp_z;
    if (sheet == 1){
      disp_out[0] = -disp_out[0];
      disp_out[1] = -disp_out[1];
      disp_out[2] = -disp_out[2];
    }
    //printf("disp_out = [%lf,%lf,%lf]\n",disp_out[0],disp_out[1],disp_out[2]);
    return disp_out;

}

std::vector<double> StrainCalc::fourierStrainDisp(double* r, double* b1, double* b2){

    // output vector
    std::vector<double> disp_out;
    disp_out.resize(3);

    double disp_x = 0.0;
    double disp_y = 0.0;
    double disp_z = 0.0;

    double k[2];
    double c[3];
    double k_base[2];
    double c_base[3];

    int idx = 0;

    // loop over [i,j] Fourier components
    for (int i = 0; i < max_k+1; ++i){
        for (int j = 0; j < max(i,1); ++j){
            //printf("[%d, %d] \n",i,j);
            // fixing factor, to account for exp(ikr) vs const term
            double ff = 2.0;
            int rot_max = 3;
            if (i == 0) {
               rot_max = 1;
               ff = 1.0;
            }

            // define Fourier component direction
            k_base[0] =  i*b1[0] + j*b2[0];
            k_base[1] =  i*b1[1] + j*b2[1];

            // define Fourier coefficients for this direction
            c_base[0] =  coeffs[0][idx]; // coeffs are in SUPERCELL basis
            c_base[1] =  coeffs[1][idx]; // later we put them in CONFIG basis
            c_base[2] =  coeffs[2][idx];

            for (int r_idx = 0; r_idx < rot_max; ++r_idx){
                double r_theta = r_idx*PI/3.0;
                k[0] = cos(r_theta)*k_base[0] - sin(r_theta)*k_base[1];
                k[1] = sin(r_theta)*k_base[0] + cos(r_theta)*k_base[1];

                //c = rot^(r_idx)*[coeffs[idx,]];
                c[0] = cos(r_theta)*c_base[0] - sin(r_theta)*c_base[1];
                c[1] = sin(r_theta)*c_base[0] + cos(r_theta)*c_base[1];
                c[2] = c_base[2];

                // the coeffs are for exp(i*k*r), so we mulitply by 2 for sin/cos
                disp_x = disp_x + ff*c[0]*sin(k[0]*r[0] + k[1]*r[1]);
                disp_y = disp_y + ff*c[1]*sin(k[0]*r[0] + k[1]*r[1]);
                disp_z = disp_z + ff*c[2]*cos(k[0]*r[0] + k[1]*r[1]);
            }

            idx++;

        }
    }


    double basis_rot = atan2(b1[1],b1[0]); // above data assumes b1 is along x axis
    double rot_disp_x = cos(basis_rot)*disp_x - sin(basis_rot)*disp_y;
    double rot_disp_y = sin(basis_rot)*disp_x + cos(basis_rot)*disp_y;
    // rotate by ~pi/2 to put into CONFIG basis
    disp_out[0] = -rot_disp_y;
    disp_out[1] = rot_disp_x;
    disp_out[2] = disp_z;
    //printf("disp_out = [%lf,%lf,%lf]\n",disp_out[0],disp_out[1],disp_out[2]);
    return disp_out;

}

std::vector< std::vector<double> > StrainCalc::fourierStrainDisp_vectorized(std::vector< std::vector<double> > r, double* b1, double* b2){

    // length of input
    int num_atoms = r.size();

    // output vector
    std::vector< std::vector<double> > disp_out;
    disp_out.resize(num_atoms);
    for (int n = 0; n < num_atoms; ++n){
      disp_out[n].resize(3);
      for (int d = 0; d < 3; ++d){
        disp_out[n][d] = 0.0;
      }
    }

    double k[2];
    double c[3];
    double k_base[2];
    double c_base[3];

    int idx = 0;

    // loop over [i,j] Fourier components
    for (int i = 0; i < max_k+1; ++i){
        for (int j = 0; j < max(i,1); ++j){
            //printf("[%d, %d] fourier comp. \n",i,j);

            // fixing factor, to account for exp(ikr) vs const term
            double ff = 2.0;
            int rot_max = 3;
            if (i == 0) {
               rot_max = 1;
               ff = 1.0;
            }

            // define Fourier component direction
            k_base[0] =  i*b1[0] + j*b2[0];
            k_base[1] =  i*b1[1] + j*b2[1];

            // define Fourier coefficients for this direction
            c_base[0] =  coeffs[0][idx]; // coeffs are in SUPERCELL basis
            c_base[1] =  coeffs[1][idx]; // later we put them in CONFIG basis
            c_base[2] =  coeffs[2][idx];

            for (int r_idx = 0; r_idx < rot_max; ++r_idx){
                double r_theta = r_idx*PI/3.0;
                k[0] = cos(r_theta)*k_base[0] - sin(r_theta)*k_base[1];
                k[1] = sin(r_theta)*k_base[0] + cos(r_theta)*k_base[1];

                //c = rot^(r_idx)*[coeffs[idx,]];
                c[0] = cos(r_theta)*c_base[0] - sin(r_theta)*c_base[1];
                c[1] = sin(r_theta)*c_base[0] + cos(r_theta)*c_base[1];
                c[2] = c_base[2];
                for (int n = 0; n < num_atoms; ++n){
                  // the coeffs are for exp(i*k*r), so we mulitply by 2 for sin/cos
                  disp_out[n][0] = disp_out[n][0] + ff*c[0]*sin(k[0]*r[n][0] + k[1]*r[n][1]);
                  disp_out[n][1] = disp_out[n][1] + ff*c[1]*sin(k[0]*r[n][0] + k[1]*r[n][1]);
                  disp_out[n][2] = disp_out[n][2] + ff*c[2]*cos(k[0]*r[n][0] + k[1]*r[n][1]);
                }
            }

            idx++;

        }
    }


    double basis_rot = atan2(b1[1],b1[0]); // above data assumes b1 is along x axis
    for (int n = 0; n < num_atoms; ++n){
      double rot_disp_x = cos(basis_rot)*disp_out[n][0] - sin(basis_rot)*disp_out[n][1];
      double rot_disp_y = sin(basis_rot)*disp_out[n][0] + cos(basis_rot)*disp_out[n][1];
      // rotate by ~pi/2 to put into CONFIG basis
      disp_out[n][0] = -rot_disp_y;
      disp_out[n][1] = rot_disp_x;
    }

    return disp_out;

}

std::vector<double> StrainCalc::fourierStrainDisp_old(std::vector<double> config_in, int sheet, int orb){


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
    //double r_scale = 0.0524;  // value for 26_25 super cell (~1.30 degrees)
    //double r_scale = 0.0666;    // value for ~1.12 degrees
    //double r_scale = 0.04237;   // value for ~1.47 degrees
    //double r_scale = 0.02404;   // value for ~2.00 degrees
    //double r_scale = 0.2;

    double r_scale = opts.getDouble("gsfe_r_scale");

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

std::vector<double> StrainCalc::supercellDisp(double* r, double* b1, double* b2, int s, double theta, int type){

  // hardcoded for graphene sandwich system
  // theta should be supplied IN DEGREES

	std::vector<double> disp;
	disp.resize(3);

  double in_plane_fac = 1.0;
  double out_of_plane_fac = 1.0;

  coeffs.resize(3);
  double y2_term = 0.0;
  double y3_term = 0.0;

  double z1_term = 0.0;
  double z2_term = 0.0;
  double z3_term = 0.0;
  double z4_term = 0.0;

  if (type == 6){ // graphene sandwich, from Zoe Zhu

    y2_term = -0.0151*pow(theta,2) + 0.0681*theta - 0.0875;
    y3_term = -0.0040*pow(theta,2) + 0.0148*theta - 0.0140;

  //printf("theta = %lf, y2_term = %lf \n",theta,y2_term);

    z1_term =  0.0161*pow(theta,2) - 0.0765*theta + 0.0158;
    z2_term =  0.0061*pow(theta,2) - 0.0255*theta + 0.0023;
    z3_term = -0.0017*pow(theta,2) + 0.0096*theta - 0.0146;
    z4_term = -0.0002*pow(theta,2) + 0.0034*theta - 0.0072;

    // relaxation symmetry for sandwich
    if (s == 2){
      out_of_plane_fac = -1.0;

    } else if (s == 1) {
      in_plane_fac = -2.0;
      out_of_plane_fac = 0.0;

    }

  } else if (type == 7) { // simple double bilayer, from Dorri Halbertal

    y2_term = -0.0319;
    y3_term = -0.0013;

    z1_term = -0.0622;
    z2_term = -0.0232;
    z3_term = -0.0060;
    z4_term = -0.0032;

    // relaxation symmetry for simple double bilayer
    if (s > 1){
      out_of_plane_fac = -1.0;
      in_plane_fac = -1.0;
    }

  } else if (type == 8) { // interior/exterior double bilayer, from Dorri Halbertal

    if (s == 0 || s == 3){ // exterior layers
      y2_term = -0.0106;
      y3_term = -0.0003;

      z1_term = -0.0017 + -0.0510;
      z2_term = -0.0004 + -0.0199;
      z3_term =  0.0003 + -0.0075;
      z4_term =  0.0002 + -0.0042;
    } else { // interior layers
      y2_term = -0.0457;
      y3_term = -0.0038;

      z1_term = -0.0510;
      z2_term = -0.0199;
      z3_term = -0.0075;
      z4_term = -0.0042;
    }

    // relaxation symmetry for double bilayer
    if (s > 1){ // if 3rd or 4th layer, invert strength
      out_of_plane_fac = -1.0;
      in_plane_fac = -1.0;
    }

  }


  std::vector<double> coeffs_x;
  coeffs_x.resize(4);
  coeffs_x[0] = in_plane_fac*0.0;
  coeffs_x[1] = in_plane_fac*0.0;
  coeffs_x[2] = in_plane_fac*0.0;
  coeffs_x[3] = in_plane_fac*0.0;

  std::vector<double> coeffs_y;
  coeffs_y.resize(4);
  coeffs_y[0] = in_plane_fac*0.0;
  coeffs_y[1] = in_plane_fac*y2_term;
  coeffs_y[2] = in_plane_fac*y3_term;
  coeffs_y[3] = in_plane_fac*0.0;

  std::vector<double> coeffs_z;
  coeffs_z.resize(4);
  coeffs_z[0] = out_of_plane_fac*z1_term;
  coeffs_z[1] = out_of_plane_fac*z2_term;
  coeffs_z[2] = out_of_plane_fac*z3_term;
  coeffs_z[3] = out_of_plane_fac*z4_term;

  coeffs[0] = coeffs_x;
  coeffs[1] = coeffs_y;
  coeffs[2] = coeffs_z;

  std::vector<double> k_base;
  k_base.resize(2);

  std::vector<double> k;
  k.resize(2);

  std::vector<double> c_base;
  c_base.resize(3);

  std::vector<double> c;
  c.resize(3);

  double disp_x = 0.0;
  double disp_y = 0.0;
  double disp_z = 0.0;


  int idx = 0;

  max_k = 2;
  // loop over [i,j] Fourier components
  for (int i = 0; i < max_k+1; ++i){
      for (int j = 0; j < max(i,1); ++j){
          //printf("[%d, %d] \n",i,j);
          // fixing factor, to account for exp(ikr) vs const term
          double ff = 2.0;
          int rot_max = 3;
          if (i == 0) {
             rot_max = 1;
             ff = 1.0;
          }

          // define Fourier component direction
          k_base[0] =  i*b1[0] + j*b2[0];
          k_base[1] =  i*b1[1] + j*b2[1];

          // define Fourier coefficients for this direction
          c_base[0] =  coeffs[0][idx]; // coeffs are in SUPERCELL basis
          c_base[1] =  coeffs[1][idx]; // later we put them in CONFIG basis
          c_base[2] =  coeffs[2][idx];

          for (int r_idx = 0; r_idx < rot_max; ++r_idx){
              double r_theta = r_idx*PI/3.0;
              k[0] = cos(r_theta)*k_base[0] - sin(r_theta)*k_base[1];
              k[1] = sin(r_theta)*k_base[0] + cos(r_theta)*k_base[1];

              //c = rot^(r_idx)*[coeffs[idx,]];
              c[0] = cos(r_theta)*c_base[0] - sin(r_theta)*c_base[1];
              c[1] = sin(r_theta)*c_base[0] + cos(r_theta)*c_base[1];
              c[2] = c_base[2];

              // the coeffs are for exp(i*k*r), so we mulitply by 2 for sin/cos
              disp_x = disp_x + ff*c[0]*sin(k[0]*r[0] + k[1]*r[1]);
              disp_y = disp_y + ff*c[1]*sin(k[0]*r[0] + k[1]*r[1]);
              disp_z = disp_z + ff*c[2]*cos(k[0]*r[0] + k[1]*r[1]);
          }

          idx++;

      }
  }


  double basis_rot = atan2(b1[1],b1[0]); // above data assumes b1 is along x axis
  double rot_disp_x = cos(basis_rot)*disp_x - sin(basis_rot)*disp_y;
  double rot_disp_y = sin(basis_rot)*disp_x + cos(basis_rot)*disp_y;

  disp[0] = rot_disp_x;
  disp[1] = rot_disp_y;
  disp[2] = disp_z;

	return disp;

}

std::vector< std::vector<double> > StrainCalc::supercellStrain(std::vector<double> pos_in, int sheet, int orb){

	std::vector< std::vector<double> > strain_here;
	// u_xx, u_xy, u_yx, u_yy;
	strain_here.resize(2);
	strain_here[0].resize(2);
	strain_here[1].resize(2);
  /*
  strain_here[0][0] = 0.0;
  strain_here[0][1] = 0.0;
  strain_here[1][0] = 0.0;
  strain_here[1][1] = 0.0;
  return strain_here;
  */

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
  //double r_scale = 0.0524; // value for 26_25 super cell (~1.30 degrees)
  //double r_scale = 0.2;
  double r_scale = opts.getDouble("gsfe_r_scale");

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
