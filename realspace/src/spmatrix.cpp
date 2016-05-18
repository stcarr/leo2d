/* 
 * File:   hstruct.cpp
 * Author: Stephen Carr
 * Based on MATKIT from the FILTLAN package
 * 
 * Created on May 18, 2016, 12:35 PM
 */

#include "spmatrix.h"

#include <stdio.h>
#include <stdlib.h>
 
#ifdef USE_MKL
	#include "mkl.h"
//#elif USE_ESSL
#endif
 
// constructors
// a constructor for an empty matrix
SpMatrix::SpMatrix() {
    nrows = 0;
    ncols = 0;
    val = NULL;
    rowIndex = NULL;
    colPointer = NULL;
    maxnnz = 0;
}


// a constructor for an nr-by-nc (sparse) matrix of zeros
SpMatrix::SpMatrix(int nr, int nc) {
    nrows = nr;
    ncols = nc;
    val = NULL;
    rowIndex = NULL;
    colPointer = new int[nc+1];
    for (int j = 0; j <= nc; j++)
        colPointer[j] = 0;
    maxnnz = 0;
}

// a constructor for an nr-by-nc (sparse) matrix stored in CSC format in store0[], rowIndex0[], colPointer0[]
// note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
// if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
// memory of size maxnnz0 is allocated, with the address pointed to by val0 and so is it for rowIndex0
SpMatrix::SpMatrix(int nr, int nc, double *val0, int *rowIndex0, int *colPointer0, int maxnnz0) {
    // set dimensions
    nrows = nr;
    ncols = nc;

    val = val0;
    rowIndex = rowIndex0;
    colPointer = colPointer0;

    maxnnz = colPointer[ncols]-colPointer[0];
    if (maxnnz < maxnnz0)
        maxnnz = maxnnz0;
    // the default maxnnz0 is 0, in which case we have maxnnz = colPointer[ncols]-colPointer[0]
}

// a destructor
SpMatrix::~SpMatrix() {
    delete [] val;
    delete [] rowIndex;
    delete [] colPointer;
}

// Sparse Matrix - Vector Multiplication

void SpMatrix::vectorMultiply(double *vec_in, double *vec_out, double alpha, double beta) {

    if (vec_out == NULL) {
        // allocate memory
        vec_out = new double[ncols];
		for (int i = 0; i < ncols; ++i){
			vec_out[i] = 0;
		}
    }

#ifdef USE_MKL
	char mv_type = 'N';
	char matdescra[6] = {'G',' ',' ','C',' ',' '};

	// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
	mkl_dcscmv(
		&mv_type,		// Specifies operator, transpose or not	
		&nrows,			// Number of rows in matrix
		&ncols,			// Number of cols in matrix
		&alpha,			// scalar ALPHA
		matdescra,		// Specifies type of matrix, *char
		val,			// nonzero elements
		rowIndex,		// row indicies
		colPointer,		// begin of col ptr
		colPointer + 1,	// end of col ptr
		vec_in,			// input vector
		&beta,			// scalar BETA
		vec_out			// output vector
		);
		
//#elif USE_ESSL
	
#else
	for (int r = 0; r < nrows; ++r){
		
		double temp_sum = 0;
		int begin = colPointer[r];
		int end = colPointer[r+1];
		
		for (int i = begin; i < end; ++i){
			temp_sum += vec_in[rowIndex[i]]*val[i];
		}
		
		vec_out[r] = alpha*temp_sum + beta*vec_out[r];
	}
#endif
}
