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
#elif USE_ESSL
	#include <essl.h> 
#endif
 
// constructors
// a constructor for an empty matrix
SpMatrix::SpMatrix() {
    nrows = 0;
    ncols = 0;
    val = NULL;
    colIndex = NULL;
    rowPointer = NULL;
    maxnnz = 0;
}


// a constructor for an nr-by-nc (sparse) matrix of zeros
SpMatrix::SpMatrix(int nr, int nc) {
    nrows = nr;
    ncols = nc;
    val = NULL;
    colIndex = NULL;
    rowPointer = new int[nc+1];
    for (int j = 0; j <= nc; j++)
        rowPointer[j] = 0;
    maxnnz = 0;
}

// a constructor for an nr-by-nc (sparse) matrix stored in CSC format in store0[], rowIndex0[], colPointer0[]
// note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
// if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
// memory of size maxnnz0 is allocated, with the address pointed to by val0 and so is it for rowIndex0
SpMatrix::SpMatrix(int nr, int nc, double *val0, int *colIndex0, int *rowPointer0, int maxnnz0) {
    // set dimensions
    nrows = nr;
    ncols = nc;

    val = val0;
    colIndex = colIndex0;
    rowPointer = rowPointer0;

    maxnnz = rowPointer[nrows]-rowPointer[0];
    if (maxnnz < maxnnz0)
        maxnnz = maxnnz0;
    // the default maxnnz0 is 0, in which case we have maxnnz = colPointer[ncols]-colPointer[0]
}

// a destructor
SpMatrix::~SpMatrix() {
    delete [] val;
    delete [] colIndex;
    delete [] rowPointer;
}

// Sparse Matrix - Vector Multiplication

void SpMatrix::vectorMultiply(double *vec_in, double *vec_out, double alpha, double beta) {
	/*
    if (vec_out == NULL) {
        // allocate memory
        vec_out = new double[nrows];
		for (int i = 0; i < nrows; ++i){
			vec_out[i] = 0;
		}
    }
	*/

#ifdef USE_MKL
	char mv_type = 'N';
	char matdescra[6] = {'G',' ',' ','C',' ',' '};

	// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
	mkl_dcsrmv(
		&mv_type,		// Specifies operator, transpose or not	
		&nrows,			// Number of rows in matrix
		&ncols,			// Number of cols in matrix
		&alpha,			// scalar ALPHA
		matdescra,		// Specifies type of matrix, *char
		val,			// nonzero elements
		colIndex,		// row indicies
		rowPointer,		// begin of col ptr
		rowPointer + 1,	// end of col ptr
		vec_in,			// input vector
		&beta,			// scalar BETA
		vec_out			// output vector
		);
		
#elif USE_ESSL

		int nz = nrows;
		double* ac;
		int* ka;
		
		double* orig_vec[nrows];
		for (int i = 0; i < nrows; ++i){
			orig_vec[i] = vec_out[i];
		}

		// converts csr format to compressed matrix storage format
		dsrsm(
			0, 			// 0 = general sparse matrix, 1 = only upper triangle (symmetric)
			val, 		// array of the values of the matrix
			colIndex, 	// array of the column indices
			rowPointer,	// row pointer array
			nrows, 		// number of rows
			&nz,		// maximum number of nonzero elements in each row 
			ac, 		// values of the converted sparse matrix
			ka, 		// column indices of the converted sparse matrix
			nrows		// size of leading dimension of the arrays
			);
			
		dsmmx(
			nrows, 		// number of rows in matrix
			nz, 		// maximum number of nonzero elements in each row
			ac, 		// values of the sparse matrix
			ka, 		// column indices of the sparse matrix
			nrows, 		// leading dimension of the arrays (>= nrows)
			vec_in, 
			vec_out
			);
			
		for (int i = 0; i < nrows; ++i){
			vec_out[i] = alpha*vec_out[i] + beta*orig_vec[i];
		}


	
#else
	for (int r = 0; r < nrows; ++r){
		
		double temp_sum = 0;
		int begin = rowPointer[r];
		int end = rowPointer[r+1];
		
		for (int i = begin; i < end; ++i){
			temp_sum += vec_in[colIndex[i]]*val[i];
		}
		
		vec_out[r] = alpha*temp_sum + beta*vec_out[r];
	}
#endif
}
