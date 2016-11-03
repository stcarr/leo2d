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


class Sdata {
        
    public:
        Sdata(std::vector<std::vector<double> >, std::vector<int>, std::vector<std::vector<double> >, std::vector<int>, std::vector<int>, int, int, int);
        Sdata(const Sdata& orig);
		Sdata();
        ~Sdata();
		std::vector<std::vector<double> > a;
        std::vector<int> max_shape, min_shape;
		int boundary_condition;
        std::vector<int> atom_types;
		std::vector<std::vector<double> > atom_pos; //atom_pos[i][x] = orbital i's displacement from the unit cell origin in the xth-dimension direction
		int mat;
		int solver_space;

};


#endif /* SDATA_H */