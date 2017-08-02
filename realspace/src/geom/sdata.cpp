/*
 * File:   sdata.cpp
 * Author: Stephen Carr
 *
 * Created on January 29, 2016, 4:28 PM
 */

 #include "sdata.h"
 #include <vector>


Sdata::Sdata(Materials::Mat _mat, std::vector<int> _min, std::vector<int> _max, int _boundary_condition, int _solver_space, int _strain_type, std::string _strain_file) {
  mat = _mat;
  a.resize(2);

  for(int i = 0; i < 2; ++i){
    a[i].resize(2);
    for(int j = 0; j < 2; ++j){
      a[i][j] = Materials::lattice(mat)[i][j];
    }
  }


  min_shape = _min;
  max_shape = _max;
  boundary_condition = _boundary_condition;
  solver_space = _solver_space;
  strain_type = _strain_type;
  strain_file = _strain_file;
}

Sdata::Sdata() {

}

Sdata::Sdata(const Sdata& orig) {
  mat = orig.mat;
  a = orig.a;
  max_shape = orig.max_shape;
  min_shape = orig.min_shape;
  boundary_condition = orig.boundary_condition;
  solver_space = orig.solver_space;
  strain_type = orig.strain_type;
  strain_file = orig.strain_file;

  supercell = orig.supercell;
  supercell_stride = orig.supercell_stride;
}

Sdata::~Sdata() {

}
