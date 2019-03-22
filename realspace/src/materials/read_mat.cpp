/*
 * File:   read_mat.cpp
 * Author:  Stephen Carr
 *
 * Created on May 21, 2018, 5:46 PM
 */

#include <string>
#include <fstream>
#include <sstream>

#include "materials/read_mat.h"

using namespace std;

array<array<double, 2>, 2>  ReadMat::getLattice(LoadedMat& mat, int sheet){
  // for now we hard code sheet = 0, as we only have bilayer systems with only one copy of monolayer information
  sheet = 0;
  return mat.intra_data[sheet].lattice;
}


int ReadMat::n_orbitals(LoadedMat& mat, int sheet){
  sheet = 0;
  return mat.intra_data[sheet].num_orbs;
}


double ReadMat::intra_search_radius(LoadedMat& mat){
  return 3.0;
}

double ReadMat::inter_search_radius(LoadedMat& mat){
  return 3.0;
}

double ReadMat::orbital_pos(LoadedMat& mat, int idx, int dim, int sheet){
  sheet = 0;
  return mat.intra_data[sheet].orb_pos[idx][dim];
}

double ReadMat::intralayer_term(int orbital_row, int orbital_col, std::array<int, 2>& vector, LoadedMat& mat, int sheet){

  sheet = 0;
  int max_R = mat.intra_data[sheet].intralayer_terms.size();
  int c = (max_R - 1)/2;
  // return 0.0 if part of vector is bigger than intralyer_terms span
  if ( abs(vector[0]) > c || abs(vector[1]) > c )
    return 0.0;

  /* debugging
  if (orbital_row == 5 && orbital_col == 5){
    printf("vector = [%d, %d], max_R = %d, c = %d, t = % lf\n",vector[0],vector[1],max_R,c, mat.intra_data[sheet].intralayer_terms[vector[0]+c][vector[1]+c][orbital_row][orbital_col]);
  }
  */
  
  return mat.intra_data[sheet].intralayer_terms[vector[0]+c][vector[1]+c][orbital_row][orbital_col];

}

double ReadMat::interlayer_term_xy_sym(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col, LoadedMat& mat){

  // swap to other direction for easy hermiticity:
  // if vector is pointing down, pass tranposed term instead (ensure Hermiticity):
  if (vector[2] < 0.0){
    std::array<double, 3> new_vector;
    for (int d = 0; d < 3; ++d){
      new_vector[d] = -vector[d];
    }
    //printf("passing reversed vector to ReadMat::interlayer_term \n");
    //printf("[%d, %d, %lf, %lf, %lf] \n",orbital_row,orbital_col,new_vector[0],new_vector[1],new_vector[2]);

    return interlayer_term_xy_sym(orbital_col, orbital_row, new_vector, angle_col, angle_row, mat);
  }

  // Decompose each orbital into correct [x,y] channels for proper symmetrization via the interlayer coupling functional
  std::vector< std::vector<int> > xy_pairs;
  xy_pairs.resize(4);
  xy_pairs[0].resize(2);
  xy_pairs[1].resize(2);
  xy_pairs[2].resize(2);
  xy_pairs[3].resize(2);

  /*
  // wannier90 convention should be:
  0 : dz^2
  1 : dxz
  2 : dyz
  3 : dx^2-y^2
  4 : dxy
  5 : top pz
  6 : top px
  7 : top py
  8 : bot pz
  9 : bot px
  10: bot py
  */

  // bottom Chalcogenide
  xy_pairs[0][0] = 6;  // top p_x
  xy_pairs[0][1] = 7;  // top p_y

  // top Chalcogenide
  xy_pairs[1][0] = 9;  // bot p_x
  xy_pairs[1][1] = 10; // bot p_y

  // d orbitals :
  xy_pairs[2][0] = 1;  // d_xz
  xy_pairs[2][1] = 2;  // d_yz

  // d orbitals :
  xy_pairs[3][0] = 3;  // d_{x^2 - y^2}
  xy_pairs[3][1] = 4;  // d_xy

  std::vector<int> row_orbs{orbital_row};
  std::vector<int> col_orbs{orbital_col};
  std::vector<double> row_facs{1.0};
  std::vector<double> col_facs{1.0};

  // Project rotated orbitals into unrotated components
  for (int idx = 0; idx < xy_pairs.size(); ++idx){

    if (xy_pairs[idx][0] == orbital_row){
      row_orbs = xy_pairs[idx];
      row_facs.resize(2);
      row_facs[0] = cos(angle_row);
      row_facs[1] = sin(angle_row);

    }

    if (xy_pairs[idx][1] == orbital_row){
      row_orbs = xy_pairs[idx];
      row_facs.resize(2);
      row_facs[0] = -sin(angle_row);
      row_facs[1] =  cos(angle_row);
    }

    if (xy_pairs[idx][0] == orbital_col){
      col_orbs = xy_pairs[idx];
      col_facs.resize(2);
      col_facs[0] = cos(angle_col);
      col_facs[1] = sin(angle_col);

    }

    if (xy_pairs[idx][1] == orbital_col){
      col_orbs = xy_pairs[idx];
      col_facs.resize(2);
      col_facs[0] = -sin(angle_col);
      col_facs[1] =  cos(angle_col);
    }

  }

  double t = 0;

  // loop over all unrotated channels and add them up!
  for(int r_idx = 0; r_idx < row_orbs.size(); ++r_idx){
    for(int c_idx = 0; c_idx < col_orbs.size(); ++c_idx){
      t += row_facs[r_idx]*col_facs[c_idx]*interlayer_term(row_orbs[r_idx], col_orbs[c_idx], vector, 0.0, 0.0, mat);
    }
  }

  return t;

}

double ReadMat::interlayer_term(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col, LoadedMat& mat){


    // if vector is pointing down, pass tranposed term instead (ensure Hermiticity):
    if (vector[2] < 0.0){
      std::array<double, 3> new_vector;
      for (int d = 0; d < 3; ++d){
        new_vector[d] = -vector[d];
      }
      //printf("passing reversed vector to ReadMat::interlayer_term \n");
      //printf("[%d, %d, %lf, %lf, %lf] \n",orbital_row,orbital_col,new_vector[0],new_vector[1],new_vector[2]);

      return interlayer_term(orbital_col, orbital_row, new_vector, angle_col, angle_row, mat);
    }


    // mat_row and mat_col are identical now
    // need to find a good system of dealing with general # of interlayer coupling functionals
    int inter_idx = 0; // needs to be generalized eventually for modeling beyond bilayer systems

    double max_r = mat.inter_data[inter_idx].max_r;
    int gridsize_x = mat.inter_data[inter_idx].gridsize_x;
    int gridsize_y = mat.inter_data[inter_idx].gridsize_y;

    std::vector<int> saved_orbs1 = mat.inter_data[inter_idx].saved_orbs1;
    std::vector<int> saved_orbs2 = mat.inter_data[inter_idx].saved_orbs2;

    double x_h = vector[0];
    double y_h = vector[1];

    // see if vector is within coupling cutoff square.
    if (x_h > max_r || x_h < -max_r || y_h > max_r || y_h < -max_r){
      return 0.0;
    }

    // see if we have a coupling pair for the two supplied orbital indices
    // if so, save the local orbital indices for the interlayer_terms dataset
    int found_pair = 0;
    int o1, o2;
    int o1_here, o2_here;

    for (int o1_idx = 0; o1_idx < (int)saved_orbs1.size(); ++o1_idx){
      o1_here = saved_orbs1[o1_idx]-1;
      for (int o2_idx = 0; o2_idx < (int)saved_orbs2.size(); ++o2_idx){
        o2_here = saved_orbs2[o2_idx]-1;
        if(o1_here == orbital_row && o2_here == orbital_col){
          o1 = o1_idx;
          o2 = o2_idx;
          found_pair = 1;
          o1_idx = 999;
          o2_idx = 999;
        }
      }
    }

    if (found_pair == 0){
      return 0.0;
    }

    // spacing between grid-points in x and y
    double dx = 2.0*max_r/( (double)(gridsize_x-1) );
    double dy = 2.0*max_r/( (double)(gridsize_y-1) );

    // find number of dx's that vector[0] is away from the start of the grid, -max_r:
    double x_tar = (x_h - (-max_r) ) / dx;
    double y_tar = (y_h - (-max_r) ) / dy;

    int x_floor = floor(x_tar);
    int y_floor = floor(y_tar);

    double x_fac = x_tar - x_floor;
    double y_fac = y_tar - y_floor;

    /*
    2d linear interpolant:
            c---d
            |   |
            a---b
    [x_floor][y_floor] is given by "a"
    */


    double a = mat.inter_data[inter_idx].interlayer_terms[o1][o2][x_floor  ][y_floor  ];
    double b = mat.inter_data[inter_idx].interlayer_terms[o1][o2][x_floor+1][y_floor  ];
    double c = mat.inter_data[inter_idx].interlayer_terms[o1][o2][x_floor  ][y_floor+1];
    double d = mat.inter_data[inter_idx].interlayer_terms[o1][o2][x_floor+1][y_floor+1];


    // 2D linear interpolation algorithm (bilinear interpolant):
    double t = ( a*(1.0-x_fac)*(1.0-y_fac) + b*x_fac*(1.0-y_fac) + c*(1.0-x_fac)*y_fac + d*x_fac*y_fac  );
    return t;

    /*
      const double nu_sigma   =  2.627;   const double nu_pi      = -0.708;
      const double R_sigma    =  3.128;   const double R_pi       =  2.923;
      const double eta_sigma  =  3.859;   const double eta_pi     =  5.724;

      const double XX_sep     = (12.29/2.0) - 3.130;

      // Corrected method for unrotating the system:
      // Assume the layer from the row has 0 twist:
      double rot_vector[3];
      // Rotate the displacement vector "backwards" so that is now defined in row_layer's coordinate system
      rot_vector[0] = cos(-theta_row)*vector[0] - sin(-theta_row)*vector[1];
      rot_vector[1] = sin(-theta_row)*vector[0] + cos(-theta_row)*vector[1];
      rot_vector[2] = vector[2];

      double tot_theta = theta_col - theta_row;
      int p_col_start = (index(orbit_col) - 5) % 3;
      double p_col[3];

      if (p_col_start == 0){ // if p_col is the x orbital
        p_col[0] = 1.0*cos(tot_theta) - 0.0*sin(tot_theta);
        p_col[1] = 1.0*sin(tot_theta) + 0.0*cos(tot_theta);
        p_col[2] = 0.0;
      } else if (p_col_start == 1){ // if p_col is the y orbital
        p_col[0] = 0.0*cos(tot_theta) - 1.0*sin(tot_theta);
        p_col[1] = 0.0*sin(tot_theta) + 1.0*cos(tot_theta);
        p_col[2] = 0.0;
      } else if (p_col_start == 2){ // if p_col is the z orbital
        p_col[0] = 0.0;
        p_col[1] = 0.0;
        p_col[2] = 1.0;
      }

      double r_sq = rot_vector[0]*rot_vector[0] + rot_vector[1]*rot_vector[1] + rot_vector[2]*rot_vector[2];
      if ( (r_sq < TMDC::inter_cutoff_radius * TMDC::inter_cutoff_radius)
                          && (std::abs(std::abs(vector[2]) - XX_sep) < 0.05) )
      {
          assert( (atom(orbit_row) == Atom::X_A && atom(orbit_col) == Atom::X_B)
                  || (atom(orbit_row) == Atom::X_B && atom(orbit_col) == Atom::X_A) ); // This assures that the inner layers only ever couple

          double r = std::sqrt(r_sq);
          // Determine character of row p orbit /
          int p_row = (index(orbit_row) - 5) % 3;

          double V_sigma = nu_sigma*std::exp(-std::pow(r/R_sigma, eta_sigma));
          double V_pi    =    nu_pi*std::exp(-std::pow(r/R_pi,    eta_pi   ));
          double sum_t = 0.0;
          for (int idx = 0; idx < 3; ++idx){
            sum_t = sum_t + p_col[idx]*((V_sigma - V_pi)*(rot_vector[p_row] * rot_vector[idx] / r_sq) + (p_row == idx ? V_pi : 0));
          }
          return sum_t;
      }
      else
          return 0;
  }
  */
}

LoadedMat ReadMat::loadMat(std::string filename){

  ifstream fileIn;
  fileIn.open(filename);

  fileIn.ignore(256, '\n');  // header info line

  LoadedMat outMat;
  std::vector< LoadedIntraData > intra_data;
  std::vector< LoadedInterData > inter_data;

  int num_intra_data;
  int num_inter_data;

  string temp;

  // Get number of data types in this file
  fileIn >> temp;
  num_intra_data = atoi(temp.c_str());
  intra_data.resize(num_intra_data);

  fileIn >> temp;
  num_inter_data = atoi(temp.c_str());
  inter_data.resize(num_inter_data);

  // loop over intra data (comes  before inter data)
  for (int idx = 0; idx < num_intra_data; ++idx){

    string name;
    array< array<double, 2>, 2> lattice;

    int num_orbs;
    vector< array<double, 3> > orb_pos;

    int num_mono_couplings;
    int big_R;
    int max_R;
    vector< vector< vector< vector<double> > > > intralayer_terms;

    fileIn >> name;   // name
    //printf("name = %s \n",name.c_str());
    fileIn.ignore(256, '\n');  // end of line
    fileIn.ignore(256, '\n');  // blank line
    fileIn.ignore(256, '\n');  // begin unit cell

    fileIn >> temp;
    lattice[0][0] = atof(temp.c_str());
    fileIn >> temp;
    lattice[0][1] = atof(temp.c_str());
    fileIn.ignore(256, '\n');

    fileIn >> temp;
    lattice[1][0] = atof(temp.c_str());
    fileIn >> temp;
    lattice[1][1] = atof(temp.c_str());
    fileIn.ignore(256, '\n');

    //printf("lattice = [%lf, %lf ; %lf, %f] \n",lattice[0][0],lattice[0][1],lattice[1][0],lattice[1][1]);

    fileIn.ignore(256, '\n');  // c axis
    fileIn.ignore(256, '\n');  // end unit cell
    fileIn.ignore(256, '\n');  // blank line

    fileIn.ignore(256, '\n');  // begin orb list

    fileIn >> num_orbs;
    //printf("num_orbs = %d \n",num_orbs);
    fileIn.ignore(256, '\n');  // end of line

    // Now read the positions of every orbital
    orb_pos.resize(num_orbs);
    for (int i = 0; i < num_orbs; ++i){
      fileIn >> orb_pos[i][0];
      fileIn >> orb_pos[i][1];
      fileIn >> orb_pos[i][2];
      //printf("orb_pos[%d] = [%lf, %lf, %lf] \n",i,orb_pos[i][0],orb_pos[i][1],orb_pos[i][2]);
    }

    fileIn.ignore(256, '\n');  // end of line
    fileIn.ignore(256, '\n');  // end orb list
    fileIn.ignore(256, '\n');  // blank line

    fileIn.ignore(256, '\n');  // begin  mono couplings
    fileIn >> num_mono_couplings;
    fileIn >> big_R;
    fileIn.ignore(256, '\n');  // end of line
    max_R = big_R*2 + 1;

    // resize intralyer_terms array
    intralayer_terms.resize(max_R);
    for (int r_i = 0; r_i < max_R; ++r_i){
      intralayer_terms[r_i].resize(max_R);
      for (int r_j = 0; r_j < max_R; ++r_j){
        intralayer_terms[r_i][r_j].resize(num_orbs);
        for (int o1 = 0; o1 < num_orbs; ++o1){
          intralayer_terms[r_i][r_j][o1].resize(num_orbs);
        }
      }
    }

    int R_i;
    int R_j;
    int o1;
    int o2;
    double t;

    for (int i = 0; i < num_mono_couplings; ++i){

      fileIn >> R_i;
      fileIn >> R_j;
      fileIn >> o1;
      fileIn >> o2;
      fileIn >> t;

      intralayer_terms[R_i + big_R][R_j + big_R][o1-1][o2-1] = t;

    }

    intra_data[idx].name = name;
    intra_data[idx].num_orbs = num_orbs;
    intra_data[idx].lattice = lattice;
    intra_data[idx].orb_pos = orb_pos;
    intra_data[idx].intralayer_terms = intralayer_terms;

    //printf("done with ReadMat::loadMat (%d orbs, %d couplings) \n",num_orbs, num_mono_couplings);
    fileIn.ignore(256, '\n');  // end of line

  }

  fileIn.ignore(256, '\n');  // end mono couplings
  fileIn.ignore(256, '\n');  // blank line

  // Now loop over all inter data
  for (int idx = 0; idx < num_inter_data; ++idx){

        string name;

        int num_orbs1;
        int num_orbs2;

        vector< array<double, 3> > orb_pos;

        int num_saved_orbs1;
        int num_saved_orbs2;

        std::vector<int> saved_orbs1;
        std::vector<int> saved_orbs2;

        double max_r;
        int gridsize_x;
        int gridsize_y;
        vector< vector< vector< vector<double> > > > interlayer_terms;

        fileIn >> name;   // name
        //printf("name = %s \n",name.c_str());
        fileIn.ignore(256, '\n');  // end of line
        fileIn.ignore(256, '\n');  // blank line
        fileIn.ignore(256, '\n');  // begin orb list

        fileIn >> num_orbs1;
        fileIn >> temp; // skip "x"
        fileIn >> num_orbs2;
        fileIn.ignore(256, '\n'); // end of line

        // Now read the positions of every orbital
        orb_pos.resize(num_orbs1 + num_orbs2);
        for (int i = 0; i < num_orbs1+num_orbs2; ++i){
          fileIn >> orb_pos[i][0];
          fileIn >> orb_pos[i][1];
          fileIn >> orb_pos[i][2];
          //printf("orb_pos[%d] = [%lf, %lf, %lf] \n",i,orb_pos[i][0],orb_pos[i][1],orb_pos[i][2]);
        }

        fileIn.ignore(256, '\n');  // end of line
        fileIn.ignore(256, '\n');  // end orb list
        fileIn.ignore(256, '\n');  // blank line

        fileIn.ignore(256, '\n');  // start saved_orb list

        fileIn >> num_saved_orbs1;
        fileIn >> temp; // skip "x"
        fileIn >> num_saved_orbs2;
        fileIn.ignore(256, '\n'); // end of line

        saved_orbs1.resize(num_saved_orbs1);
        saved_orbs2.resize(num_saved_orbs2);

        for (int i = 0; i < num_saved_orbs1; ++i){
          fileIn >> saved_orbs1[i];
        }
        fileIn.ignore(256, '\n'); // end of line

        for (int i = 0; i < num_saved_orbs2; ++i){
          fileIn >> saved_orbs2[i];
        }
        fileIn.ignore(256, '\n'); // end of line
        fileIn.ignore(256, '\n'); // end saved_orb list

        fileIn.ignore(256, '\n');  // blank line
        fileIn.ignore(256, '\n');  // begin inter couplings

        fileIn >> max_r;
        fileIn >> gridsize_x;
        fileIn >> gridsize_y;
        fileIn.ignore(256, '\n'); // end of line
        fileIn.ignore(256, '\n'); // blank line

        interlayer_terms.resize(num_saved_orbs1);

        //printf("max_r = %lf, gridsize_x = %d, gridsize_y = %d \n",max_r, gridsize_x, gridsize_y);

        int temp_o1;
        int temp_o2;

        for (int o1 = 0; o1 < num_saved_orbs1; ++o1){
          interlayer_terms[o1].resize(num_saved_orbs2);
          for (int o2 = 0; o2 < num_saved_orbs2; ++o2){
            interlayer_terms[o1][o2].resize(gridsize_x);

            fileIn >> temp_o1;
            fileIn >> temp_o2;

            if (temp_o1 != saved_orbs1[o1]){
              throw std::runtime_error("ReadMat::LoadMat Interlayer coupling orbital index mismatch in orbital 1! \n");
            }

            if (temp_o2 != saved_orbs2[o2]){
              throw std::runtime_error("ReadMat::LoadMat Interlayer coupling orbital index mismatch in orbital 2! \n");
            }

            for (int gx = 0; gx < gridsize_x; ++gx){
              interlayer_terms[o1][o2][gx].resize(gridsize_y);
              for (int gy = 0; gy < gridsize_y; ++gy){
                fileIn >> interlayer_terms[o1][o2][gx][gy];
              }
            }

            fileIn.ignore(256, '\n'); // end of line
            fileIn.ignore(256, '\n'); // blank line

        }
      }

      inter_data[idx].name = name;
      inter_data[idx].num_orbs1 = num_orbs1;
      inter_data[idx].num_orbs2 = num_orbs2;
      inter_data[idx].num_saved_orbs1 = num_saved_orbs1;
      inter_data[idx].num_saved_orbs2 = num_saved_orbs2;
      inter_data[idx].saved_orbs1 = saved_orbs1;
      inter_data[idx].saved_orbs2 = saved_orbs2;
      inter_data[idx].max_r = max_r;
      inter_data[idx].gridsize_x = gridsize_x;
      inter_data[idx].gridsize_y = gridsize_y;
      inter_data[idx].orb_pos = orb_pos;
      inter_data[idx].interlayer_terms = interlayer_terms;



      //printf("done with ReadMat::loadMat (%d orbs, %d couplings) \n",num_orbs, num_mono_couplings);
  }

  fileIn.close();

  outMat.num_intra_data = num_intra_data;
  outMat.num_inter_data = num_inter_data;
  outMat.intra_data = intra_data;
  outMat.inter_data = inter_data;
  return outMat;

}
