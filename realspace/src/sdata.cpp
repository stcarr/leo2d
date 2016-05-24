/* 
 * File:   sdata.cpp
 * Author: Stephen Carr
 *
 * Created on January 29, 2016, 4:28 PM
 */
 
 #include "sdata.h"
 #include <vector>

Sdata::Sdata(std::vector<std::vector<double> > _a, std::vector<int> _types, std::vector<std::vector<double> > _pos, std::vector<int> _min, std::vector<int> _max, int _mat) {
    a = _a;
    min_shape = _min;
    max_shape = _max;
    atom_types = _types;
    atom_pos = _pos;
	mat = _mat;
           
}

Sdata::Sdata() {

}

Sdata::Sdata(const Sdata& orig) {
    a = orig.a;
    max_shape = orig.max_shape;
    min_shape = orig.min_shape;
    atom_types = orig.atom_types;
    atom_pos = orig.atom_pos;
	mat = orig.mat;
}

Sdata::~Sdata() {

}