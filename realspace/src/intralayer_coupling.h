/* 
 * File:   intralayer_coupling.h
 * Author: Stephen
 * 
 * Created on February 5, 2016, 3:33 PM
 */

#include <math.h>
#include <vector>

#ifndef INTRALAYER_COUPLING_H
#define INTRALAYER_COUPLING_H

std::vector<std::vector<double> > intralayer_terms(int (&grid_index)[3], int);
std::vector<std::vector<double> > intralayer_graphene(int (&grid_index)[3]);

#endif /* INTRALAYER_COUPLING_H */