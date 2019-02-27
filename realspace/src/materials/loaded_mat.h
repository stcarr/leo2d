/*
 * File:   loaded_mat.h
 * Author:  Stephen Carr
 *
 * Created on May 21, 2018, 5:46 PM
 */

#ifndef loaded_mat_h
#define loaded_mat_h

#include <array>
#include <vector>

#include "tools/numbers.h"

struct LoadedIntraData{

  std::string name;
  std::array<std::array<double, 2>, 2> lattice;

  int num_orbs;
  std::vector< std::array<double, 3> > orb_pos;

  // intralyer terms, intralyer_terms[R_i][R_j][o1][o2]
  std::vector< std::vector< std::vector< std::vector<double> > > > intralayer_terms;

};

struct LoadedInterData{

};

struct LoadedMat{

          int num_intra_data;
          int num_inter_data;

          std::vector< LoadedIntraData > intra_data;
          std::vector< LoadedInterData > inter_data;

};

#endif
