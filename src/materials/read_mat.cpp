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

array<array<double, 2>, 2>  ReadMat::getLattice(LoadedMat& mat, Sdata& sheet_data){

  for (int i = 0; i < mat.intra_data.size(); ++i){
    if (sheet_data.lmat_name.compare(mat.intra_data[i].name) == 0){
      return mat.intra_data[i].lattice;
    }
  }

  throw std::runtime_error("ReadMat::getLattice ERROR!! Material not found in loaded material file: " + sheet_data.lmat_name);

}


int ReadMat::n_orbitals(LoadedMat& mat,  Sdata& sheet_data){

  for (int i = 0; i < mat.intra_data.size(); ++i){
    if (sheet_data.lmat_name.compare(mat.intra_data[i].name) == 0){
      return mat.intra_data[i].num_orbs;
    }
  }

  throw std::runtime_error("ReadMat::n_orbitals ERROR!! Material not found in loaded material file: " + sheet_data.lmat_name);

}


double ReadMat::intra_search_radius(LoadedMat& mat){
  return 3.0;
}

double ReadMat::inter_search_radius(LoadedMat& mat){
  return 3.0;
}

double ReadMat::orbital_pos(LoadedMat& mat, int idx, int dim, Sdata& sheet_data){

  for (int i = 0; i < mat.intra_data.size(); ++i){
    if (sheet_data.lmat_name.compare(mat.intra_data[i].name) == 0){
      return mat.intra_data[i].orb_pos[idx][dim];
    }
  }

  throw std::runtime_error("ReadMat::orbital_pos ERROR!! Material not found in loaded material file: " + sheet_data.lmat_name);

}

double ReadMat::intralayer_term(int orbital_row, int orbital_col, std::array<int, 2>& vector, LoadedMat& mat,  Sdata& sheet_data){

  int sheet = -1;
  for (int i = 0; i < mat.intra_data.size(); ++i){
    if (sheet_data.lmat_name.compare(mat.intra_data[i].name) == 0){
      sheet = i;
    }
  }

  if (sheet == -1){
    throw std::runtime_error("ReadMat::intralayer_term ERROR!! Material not found in loaded material file: " + sheet_data.lmat_name);
  }

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

double ReadMat::interlayer_term_basic_c3_sym(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col,
                                       LoadedMat& mat, Sdata& sheet_data1, Sdata& sheet_data2){

  double t_temp = 0;
  std::array<double, 3> vector_here{0.0, 0.0, 0.0};
  vector_here[2] = vector[2];

  int r_max = 3;//3;
  int c2_max = 2;//2;
  double inter_r_cutoff = 5.0;

  double r_sq = vector[0]*vector[0] + vector[1]*vector[1];

  double t_here;
  // 3 fold rotational symm
  for (int r_idx = 0; r_idx < r_max; ++r_idx){

    double rot_theta = ((double)r_idx)*2.0*M_PI/3.0;
    vector_here[0] = cos(rot_theta)*vector[0] - sin(rot_theta)*vector[1];
    vector_here[1] = sin(rot_theta)*vector[0] + cos(rot_theta)*vector[1];

    // c2 sym, swaps A<->B orbitals
    for (int c2_idx = 0; c2_idx < c2_max; ++c2_idx){

      if (c2_idx == 0){
        t_here = interlayer_term(orbital_row, orbital_col, vector_here, angle_row, angle_col, mat, sheet_data1, sheet_data2);
      } else {


        // c2 operation
        vector_here[0] = -vector_here[0];
        vector_here[1] = -vector_here[1];

        // swap the orbitals
        int o_row_here;
        int o_col_here;

        if (orbital_row == 0)
          o_row_here = 1;
        else
          o_row_here = 0;
        if (orbital_col == 0)
          o_col_here = 1;
        else
          o_col_here = 0;

        t_here = interlayer_term(o_row_here, o_col_here, vector_here, angle_row, angle_col, mat, sheet_data1, sheet_data2);



      }

      t_temp += t_here;

      if (abs(t_here) > .2){
      //printf("[%d,%d] [%lf, %lf, %lf]: adding %lf \n",r_idx,m_idx,vector[0],vector[1],vector[2],t_here);
      }
    }
  }

  if (r_sq < inter_r_cutoff*inter_r_cutoff) {
    //printf("[%lf, %lf, %lf]: returning %lf \n",vector[0],vector[1],vector[2],t_temp/6.0);
  } else{
    t_temp = 0.0;
  }
  return t_temp/((double) r_max*c2_max);


}


double ReadMat::interlayer_term_xy_sym(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col,
                                       LoadedMat& mat, Sdata& sheet_data1, Sdata& sheet_data2){

  // swap to other direction for easy hermiticity:
  // if vector is pointing down, pass tranposed term instead (ensure Hermiticity):
  if (vector[2] < 0.0){
    std::array<double, 3> new_vector;
    for (int d = 0; d < 3; ++d){
      new_vector[d] = -vector[d];
    }
    //printf("passing reversed vector to ReadMat::interlayer_term \n");
    //printf("[%d, %d, %lf, %lf, %lf] \n",orbital_row,orbital_col,new_vector[0],new_vector[1],new_vector[2]);

    return interlayer_term_xy_sym(orbital_col, orbital_row, new_vector, angle_col, angle_row, mat, sheet_data2, sheet_data1);
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
      t += row_facs[r_idx]*col_facs[c_idx]*interlayer_term(row_orbs[r_idx], col_orbs[c_idx], vector, 0.0, 0.0, mat, sheet_data1, sheet_data2);
    }
  }

  return t;

}

double ReadMat::interlayer_term(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col,
                                LoadedMat& mat, Sdata& sheet_data1, Sdata& sheet_data2){

    // if vector is pointing down, pass tranposed term instead (ensure Hermiticity):
    if (vector[2] < 0.0){
      std::array<double, 3> new_vector;
      for (int d = 0; d < 3; ++d){
        new_vector[d] = -vector[d];
      }
      //printf("passing reversed vector to ReadMat::interlayer_term \n");
      //printf("[%d, %d, %lf, %lf, %lf] \n",orbital_row,orbital_col,new_vector[0],new_vector[1],new_vector[2]);

      return interlayer_term(orbital_col, orbital_row, new_vector, angle_col, angle_row, mat, sheet_data2, sheet_data1);
    }

    // look for bilayer info
    std::string bilayer_name = sheet_data1.lmat_name + '-' + sheet_data2.lmat_name;

    int inter_idx = -1;

    for (int i = 0; i < mat.inter_data.size(); ++i){
      if (bilayer_name.compare(mat.inter_data[i].name) == 0){
        inter_idx = i;
      }
    }

    if (inter_idx == -1){
      throw std::runtime_error("ReadMat::interlayer_term ERROR!! Bilayer data not found in loaded material file: " + bilayer_name);
    }

    //int inter_idx = 0;


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

LoadedMat ReadMat::loadMat(std::string filename, std::vector<Sdata> sdata, double intra_rsq_cutoff){

  // do some parsing of sdata to get list of monolayer info and interayer info desired

  ifstream fileIn;
  fileIn.open(filename);

  fileIn.ignore(256, '\n');  // header info line

  LoadedMat outMat;
  std::vector< LoadedIntraData > intra_data;
  std::vector< LoadedInterData > inter_data;

  int num_intra_data;
  int num_inter_data;

  string temp;

  string intra_header("MonolayerTBH_Data");
  string inter_header("InterlayerTBH_Data");


  // Get number of data types in this file
  /*
  fileIn >> temp;
  num_intra_data = atoi(temp.c_str());
  intra_data.resize(num_intra_data);

  fileIn >> temp;
  num_inter_data = atoi(temp.c_str());
  inter_data.resize(num_inter_data);
  */

  /*
  // loop over intra data (comes  before inter data)
  for (int idx = 0; idx < num_intra_data; ++idx){
  */

  while(!fileIn.eof()){

    fileIn >> temp;

    if (temp.compare(intra_header) == 0){

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
        fileIn >> o1; // FROM o1
        fileIn >> o2; // TO o2
        fileIn >> t;

        double delta_x = R_i*lattice[0][0] + R_j*lattice[1][0] + orb_pos[o2-1][0] - orb_pos[o1-1][0];
        double delta_y = R_i*lattice[0][1] + R_j*lattice[1][1] + orb_pos[o2-1][1] - orb_pos[o1-1][1];
        double delta_rsq = delta_x*delta_x + delta_y*delta_y;

        if (delta_rsq < intra_rsq_cutoff){
          //printf("[%d, %d, %d, %d] = %lf added \n",R_i, R_j, o1,o2,t);
          intralayer_terms[R_i + big_R][R_j + big_R][o1-1][o2-1] = t;
        }

      }

      LoadedIntraData temp_intradata;

      temp_intradata.name = name;
      temp_intradata.num_orbs = num_orbs;
      temp_intradata.lattice = lattice;
      temp_intradata.orb_pos = orb_pos;
      temp_intradata.intralayer_terms = intralayer_terms;
      intra_data.push_back(temp_intradata);

      /*
      intra_data[idx].name = name;
      intra_data[idx].num_orbs = num_orbs;
      intra_data[idx].lattice = lattice;
      intra_data[idx].orb_pos = orb_pos;
      intra_data[idx].intralayer_terms = intralayer_terms;
      */

      //printf("done with ReadMat::loadMat (%d orbs, %d couplings) \n",num_orbs, num_mono_couplings);
      //fileIn.ignore(256, '\n');  // end of line

      // end monolayer info parsing
    } else if (temp.compare(inter_header) == 0){

      // Now loop over all inter data
      //for (int idx = 0; idx < num_inter_data; ++idx){

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
            throw std::runtime_error("ReadMat::LoadMat ERROR!! Interlayer coupling orbital index mismatch in orbital 1! \n");
          }

          if (temp_o2 != saved_orbs2[o2]){
            throw std::runtime_error("ReadMat::LoadMat ERROR!! Interlayer coupling orbital index mismatch in orbital 2! \n");
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

    LoadedInterData temp_interdata;

    temp_interdata.name = name;
    temp_interdata.num_orbs1 = num_orbs1;
    temp_interdata.num_orbs2 = num_orbs2;
    temp_interdata.num_saved_orbs1 = num_saved_orbs1;
    temp_interdata.num_saved_orbs2 = num_saved_orbs2;
    temp_interdata.saved_orbs1 = saved_orbs1;
    temp_interdata.saved_orbs2 = saved_orbs2;
    temp_interdata.max_r = max_r;
    temp_interdata.gridsize_x = gridsize_x;
    temp_interdata.gridsize_y = gridsize_y;
    temp_interdata.orb_pos = orb_pos;
    temp_interdata.interlayer_terms = interlayer_terms;

    inter_data.push_back(temp_interdata);
    /*
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
    */

    //printf("done with ReadMat::loadMat (%d orbs, %d couplings) \n",num_orbs, num_mono_couplings);
    } // end interlayer info parsing

  } // end fileIn.eof() loop

  fileIn.close();

  // check if we got everything we need!
  // first monolayer data
  for (int i = 0; i < sdata.size(); ++ i){

    string name_here = sdata[i].lmat_name;

    //printf("lmat_name: %s \n",name_here.c_str());

    int mat_found = 0;

    for (int j = 0; j < intra_data.size(); ++j){
      if (name_here.compare(intra_data[j].name) == 0){
        mat_found = 1;
        break;
      }
    }

    if (mat_found == 0){
      throw std::runtime_error("ReadMat::LoadMat ERROR!! Material not found in loaded material file: " + name_here);
    }

  }

  // interlayer (bilayer) data
  for (int i = 0; i < sdata.size() - 1; ++ i){

    string name_here_bot = sdata[i].lmat_name;
    string name_here_top = sdata[i+1].lmat_name;
    string bilayer_name_here = name_here_bot + '-' + name_here_top;

    //printf("bilayer_name_here: %s \n",bilayer_name_here.c_str());

    int mat_found = 0;

    for (int j = 0; j < inter_data.size(); ++j){
      if (bilayer_name_here.compare(inter_data[j].name) == 0){
        mat_found = 1;
        break;
      }
    }

    if (mat_found == 0){
      throw std::runtime_error("ReadMat::LoadMat ERROR!! Interlayer data not found in loaded material file: " + bilayer_name_here);
    }
  }


  outMat.num_intra_data = intra_data.size();
  outMat.num_inter_data = inter_data.size();
  outMat.intra_data = intra_data;
  outMat.inter_data = inter_data;
  return outMat;

}

void ReadMat::monolayerSym(LoadedMat& mat, int target_index){


  std::vector< std::vector< std::vector< std::vector<double> > > > intralayer_terms = mat.intra_data[target_index].intralayer_terms;
  int max_r = intralayer_terms.size();
  //max_R = big_R*2 + 1;
  int big_R = (max_r - 1) / 2;

  std::vector< std::vector< std::vector< std::vector<double> > > > intralayer_terms_sym;
  std::vector< std::vector< std::vector< std::vector<double> > > > intralayer_terms_sym_temp;
  intralayer_terms_sym = intralayer_terms;
  intralayer_terms_sym_temp = intralayer_terms;

  //intralayer_terms[R_x+big_R][R_y+big_R][o1][o2] = t;
  int R_x_size = intralayer_terms.size();
  int R_y_size = intralayer_terms[0].size();
  int o1_size  = intralayer_terms[0][0].size();
  int o2_size  = intralayer_terms[0][0][0].size();

  int num_sym_ops = 2;

  // Row-major, A_C3[0][1] = first row, second col
  std::vector <std::vector<int> > A_C3 = {{ -1 , -1},{ 1 , 0}};
  std::vector <std::vector<int> > sitemap_C3 = { {0,0,0,1}, {1,-1,0,1}};
  int C3_order = 2;

  std::vector <std::vector<int> >  A_M = {{ 0 , -1},{ -1 , 0}};
  std::vector <std::vector<int> >  sitemap_M = { {1,0,0,1}, {0,0,0,1}};
  int M_order = 1;

  int tot_order = 1 + C3_order + M_order;
  //double rescale = 1.0/( (double) tot_order);

  /*
  intralayer_terms_sym.resize(R_x_size);

  for (int R_x_idx = 0; R_x_idx < R_x_size; R_x_idx++){
    intralayer_terms_sym[R_x_idx].resize(R_y_size);
    for (int R_y_idx = 0; R_y_idx < R_y_size; R_y_idx++){
      intralayer_terms_sym[R_x_idx][R_y_idx].resize(o1_size);
      for (int o1 = 0; o1 < o1_size; o1++){
        intralayer_terms_sym[R_x_idx][R_y_idx][o1].resize(o2_size);

        for (int o2 = 0; o2 < o2_size; o2++){

          double t_here = intralayer_terms[R_x_idx][R_y_idx][o1][o2];
          intralayer_terms_sym[R_x_idx][R_y_idx][o1][o2] = t_here;

        }

      }
    }
  }
  */


  // go over each symmetry condition
  for (int sym_idx = 0; sym_idx < num_sym_ops; sym_idx++){

    std::vector <std::vector<int> >  A;
    std::vector <std::vector<int> >  sitemap;
    int order;

    if (sym_idx == 0) {
      // C3 symmetry
      A = A_C3;
      sitemap = sitemap_C3;
      order = C3_order;

    } else {
      // Mirror symmetry
      A = A_M;
      sitemap = sitemap_M;
      order = M_order;

    }

    for (int R_x_idx = 0; R_x_idx < R_x_size; R_x_idx++){
      for (int R_y_idx = 0; R_y_idx < R_y_size; R_y_idx++){
        for (int o1 = 0; o1 < o1_size; o1++){
          for (int o2 = 0; o2 < o2_size; o2++){

            int sym_R_x;
            int sym_R_y;
            int sym_o1;
            int sym_o2;

            double t_here = intralayer_terms_sym[R_x_idx][R_y_idx][o1][o2];
            int R_x_h = R_x_idx - big_R;
            int R_y_h = R_y_idx - big_R;
            int o1_h = o1;
            int o2_h = o2;

            int R_x_new;
            int R_y_new;
            int o1_new;
            int o2_new;

            for (int idx = 0; idx < order; ++idx){

              //printf("idx = %d \n", idx);

              // hopping: FROM o1, TO o2
              R_x_new = A[0][0]*R_x_h + A[0][1]*R_y_h + sitemap[o2_h][1] - sitemap[o1_h][1];
              R_y_new = A[1][0]*R_x_h + A[1][1]*R_y_h + sitemap[o2_h][2] - sitemap[o1_h][2];
              o1_new = sitemap[o1_h][0];
              o2_new = sitemap[o2_h][0];

              //printf("new loc = [%d][%d][%d][%d] \n", R_x_new+big_R, R_y_new+big_R, o1_new, o2_new);

              if (R_x_new + big_R < R_x_size && R_x_new + big_R > 0 && R_y_new + big_R < R_y_size && R_y_new + big_R > 0) {

                double t_new = intralayer_terms_sym[R_x_new + big_R][R_y_new + big_R][o1_new][o2_new];
                if (abs(t_new) > 2){
                  //printf("[%d, %d]: averaging %lf with %lf: [%d, %d, %d, %d] to [%d, %d, %d, %d] \n",sym_idx, idx, t_here, t_new, R_x_h, R_y_h, o1_h+1, o2_h+1, R_x_new, R_y_new, o1_new+1, o2_new+1);
                }

                t_here += t_new;

              }

              R_x_h = R_x_new;
              R_y_h = R_y_new;
              o1_h = o1_new;
              o2_h = o2_new;

            }

            intralayer_terms_sym_temp[R_x_idx][R_y_idx][o1][o2] = t_here/( (double) 1 + order);

          }


        }
      }
    }

    intralayer_terms_sym = intralayer_terms_sym_temp;

  }

  /*
  for (int R_x_idx = 0; R_x_idx < R_x_size; R_x_idx++){
    for (int R_y_idx = 0; R_y_idx < R_y_size; R_y_idx++){
      for (int o1 = 0; o1 < o1_size; o1++){
        for (int o2 = 0; o2 < o2_size; o2++){

          double t_here = intralayer_terms_sym[R_x_idx][R_y_idx][o1][o2];
          if (abs(t_here) > 1e-10){
            printf("symmetrized [%d, %d, %d, %d] = %.20f \n", R_x_idx - big_R, R_y_idx - big_R, o1+1,o2+1, t_here);
          }

        }
      }
    }
  }
  */

  mat.intra_data[target_index].intralayer_terms = intralayer_terms_sym;

}
