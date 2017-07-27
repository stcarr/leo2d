/*
 * File:   sdata.h
 * Author: Stephen Carr
 *
 * Created on January 29, 2016, 4:28 PM
 */


#ifndef SDATA_H
#define SDATA_H
#include <vector>
#include <string>
#include "materials/materials.h"


class Sdata {

  public:

    Sdata(Materials::Mat, std::vector<int>, std::vector<int>, int, int, int, std::string);
    Sdata(const Sdata& orig);
		Sdata();
    ~Sdata();
    Materials::Mat mat;
    std::vector<std::vector<double> > a;
    std::vector<int> max_shape, min_shape;
		int boundary_condition;
    std::vector< std::vector<double> > supercell;
    std::vector< std::vector<int> > supercell_stride;
		int solver_space;
		int strain_type;
		std::string strain_file;

};


#endif /* SDATA_H */