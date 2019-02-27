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

array<array<double, 2>, 2>  ReadMat::getLattice(LoadedMat mat, int sheet){
  return mat.intra_data[sheet].lattice;
}


int ReadMat::n_orbitals(LoadedMat mat, int sheet){
  return mat.intra_data[sheet].num_orbs;
}


double ReadMat::intra_search_radius(LoadedMat mat){
  return 3.0;
}

double ReadMat::inter_search_radius(LoadedMat mat){
  return 10.0;
}

double ReadMat::orbital_pos(LoadedMat mat, int idx, int dim, int sheet){
  return mat.intra_data[sheet].orb_pos[idx][dim];
}

double ReadMat::intralayer_term(int orbital_row, int orbital_col, std::array<int, 2>& vector, LoadedMat mat, int sheet){

  int max_R = mat.intra_data[sheet].intralayer_terms.size();
  int c = (max_R - 1)/2;
  // return 0.0 if part of vector is bigger than intralyer_terms span
  if ( abs(vector[0]) > c || abs(vector[1] > c) )
    return 0.0;

  return mat.intra_data[sheet].intralayer_terms[vector[0]+c][vector[1]+c][orbital_row][orbital_col];

}

double ReadMat::interlayer_term(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col, LoadedMat mat_row, LoadedMat mat_col){

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

      intralayer_terms[R_i + big_R][R_j + big_R][o1][o2] = t;

    }

    intra_data[idx].name = name;
    intra_data[idx].num_orbs = num_orbs;
    intra_data[idx].lattice = lattice;
    intra_data[idx].orb_pos = orb_pos;
    intra_data[idx].intralayer_terms = intralayer_terms;

    //printf("done with ReadMat::loadMat (%d orbs, %d couplings) \n",num_orbs, num_mono_couplings);

  }

  // Now loop over all inter data
  for (int idx = 0; idx < num_inter_data; ++idx){

  }

  fileIn.close();

  outMat.num_intra_data = num_intra_data;
  outMat.num_inter_data = num_inter_data;
  outMat.intra_data = intra_data;
  outMat.inter_data = inter_data;
  return outMat;

}
