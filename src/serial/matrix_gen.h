
/*
 * File:   matrix_gen.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:39 PM
 */

#ifndef MATRIX_GEN_H
#define MATRIX_GEN_H

#include "serial/locality_serial.h"
#include "matrix/dmatrix.h"

// IMPORTANT, only use getLeoMatrix as a copy constructor! i.e. only do:
//     Dmatrix out_matrix = getLeoMatrix(hstruct_input_file);
// If you use assignment (i.e. declare out_matrix then assign it at a later line)
// then we will get problems due to a lack of a proper assignment operator in DMatrix!
DMatrix getLeoMatrix(string hstruct_input_file);

#endif /* MATRIX_GEN_H */
