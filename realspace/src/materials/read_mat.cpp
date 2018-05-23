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
  return mat.lattice;
}


int ReadMat::n_orbitals(LoadedMat mat, int sheet){
  return mat.num_orbs;
}


double ReadMat::intra_search_radius(LoadedMat mat){
  return 3.0;
}

double ReadMat::inter_search_radius(LoadedMat mat){
  return 10.0;
}

double ReadMat::orbital_pos(LoadedMat mat, int idx, int dim, int sheet){
  return mat.orb_pos[idx][dim];
}

double ReadMat::intralayer_term(int orbital_row, int orbital_col, std::array<int, 2>& vector, LoadedMat mat, int sheet){

  int max_R = mat.intralayer_terms.size();
  int c = (max_R - 1)/2;
  // return 0.0 if part of vector is bigger than intralyer_terms span
  if ( abs(vector[0]) > c || abs(vector[1] > c) )
    return 0.0;

  return mat.intralayer_terms[vector[0]+c][vector[1]+c][orbital_row][orbital_col];

}

double ReadMat::interlayer_term(int orbital_row, int orbital_col, std::array<double, 3>& vector, double angle_row, double angle_col, LoadedMat mat_row, LoadedMat mat_col){
  return -1.0;
}

LoadedMat ReadMat::loadMat(std::string filename){

  ifstream fileIn;
  fileIn.open(filename);

  LoadedMat outMat;

  string name;
  array< array<double, 2>, 2> lattice;

  int num_orbs;
  vector< array<double, 3> > orb_pos;

  int num_mono_couplings;
  int big_R;
  int max_R;
  vector< vector< vector< vector<double> > > > intralayer_terms;


  string temp;

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

  outMat.name = name;
  outMat.num_orbs = num_orbs;
  outMat.lattice = lattice;
  outMat.orb_pos = orb_pos;
  outMat.intralayer_terms = intralayer_terms;

  printf("done with ReadMat::loadMat (%d orbs, %d couplings) \n",num_orbs, num_mono_couplings);

  fileIn.close();
  return outMat;

}
