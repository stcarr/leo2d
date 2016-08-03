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
	type = 0;
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
	type = 0;
}

// a constructor for an nr-by-nc (sparse) matrix stored in CSC format in store0[], rowIndex0[], colPointer0[]
// note that it performs a shallow copy of store0, rowIndex0, and colPointer0, and the memory will be freed by the destructor
// if maxnnz0 is greater than the number of nonzero elements, then it is assumed that
// memory of size maxnnz0 is allocated, with the address pointed to by val0 and so is it for rowIndex0
SpMatrix::SpMatrix(int nr, int nc, double *val0, int *colIndex0, int *rowPointer0, int maxnnz0) {
    
	type = 0;


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
   
#ifdef USE_ESSL
	nz = nrows+3;

	double ac[nrows][nrows];
	int ka[nrows][nrows];
	//ac = new double[nrows*nrows];
	//ka = new int[nrows*nrows];
	printf("calling essl_dsrsm. \n");
	// converts csr format to compressed matrix storage format
	printf("nz in = %d \n",nz);	
	dsrsm(
			0, 		// 0 = general sparse matrix, 1 = only upper triangle (symmetric)
			val, 		// array of the values of the matrix
			colIndex, 	// array of the column indices
			rowPointer,	// row pointer array
			nrows, 		// number of rows
			nz,		// maximum number of nonzero elements in each row 
			ac, 		// values of the converted sparse matrix
			ka, 		// column indices of the converted sparse matrix
			nrows		// size of leading dimension of the arrays
			);
	printf("nz out = %d \n",nz);
	for(int i = 0; i < nrows; ++i){
		for(int j = 0;j < nrows; ++j){
			printf("ka[%d][%d] = %d: ac[%d][%d] = %lf ([%d][%d] = %lf) \n",j,i,ka[j][i],j,i,ac[j][i],i,j,val[i*nrows + j]);
		}
	}

#endif

}

// Complex constructor
SpMatrix::SpMatrix(int nr, int nc, std::complex<double> *val_c0, int *colIndex0, int *rowPointer0, int maxnnz0) {
    
	type = 1;

	// set dimensions
	nrows = nr;
	ncols = nc;

	val_c = val_c0;
	colIndex = colIndex0;
	rowPointer = rowPointer0;

	maxnnz = rowPointer[nrows]-rowPointer[0];
	if (maxnnz < maxnnz0)
		maxnnz = maxnnz0;
    	// the default maxnnz0 is 0, in which case we have maxnnz = colPointer[ncols]-colPointer[0]
		
}


// a destructor
SpMatrix::~SpMatrix() {

	if (type == 0){
		delete [] val;
	}
	
	if (type == 1){
		delete [] val_c;
	}
	
	delete [] colIndex;
	delete [] rowPointer;
}

// Real setup
void SpMatrix::setup(int nr, int nc, double *val0, int *colIndex0, int *rowPointer0, int maxnnz0){   
	
	type = 0;


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
   
	#ifdef USE_ESSL
		nz = nrows+3;

		double ac[nrows][nrows];
		int ka[nrows][nrows];
		//ac = new double[nrows*nrows];
		//ka = new int[nrows*nrows];
		printf("calling essl_dsrsm. \n");
		// converts csr format to compressed matrix storage format
		printf("nz in = %d \n",nz);	
		dsrsm(
				0, 		// 0 = general sparse matrix, 1 = only upper triangle (symmetric)
				val, 		// array of the values of the matrix
				colIndex, 	// array of the column indices
				rowPointer,	// row pointer array
				nrows, 		// number of rows
				nz,		// maximum number of nonzero elements in each row 
				ac, 		// values of the converted sparse matrix
				ka, 		// column indices of the converted sparse matrix
				nrows		// size of leading dimension of the arrays
				);
		printf("nz out = %d \n",nz);
		for(int i = 0; i < nrows; ++i){
			for(int j = 0;j < nrows; ++j){
				printf("ka[%d][%d] = %d: ac[%d][%d] = %lf ([%d][%d] = %lf) \n",j,i,ka[j][i],j,i,ac[j][i],i,j,val[i*nrows + j]);
			}
		}
	#endif
}

// Complex setup
void SpMatrix::setup(int nr, int nc, std::complex<double> *val_c0, int *colIndex0, int *rowPointer0, int maxnnz0) {
    
	type = 1;

	// set dimensions
	nrows = nr;
	ncols = nc;

	val_c = val_c0;
	colIndex = colIndex0;
	rowPointer = rowPointer0;

	maxnnz = rowPointer[nrows]-rowPointer[0];
	if (maxnnz < maxnnz0)
		maxnnz = maxnnz0;
    	// the default maxnnz0 is 0, in which case we have maxnnz = colPointer[ncols]-colPointer[0]
		
}


// Sparse Matrix - Vector Multiplication
void SpMatrix::vectorMultiply(double *vec_in, double *vec_out, double alpha, double beta) {

	if (vec_out == NULL) {
		// allocate memory
		vec_out = new double[nrows];
		for (int i = 0; i < nrows; ++i){
			vec_out[i] = 0.0;
			}
	}
		
	/*
	for (int i = 0; i < nrows; ++i){
		printf("vec_in[%d] = %lf \n",i,vec_in[i]);
		printf("and vec_out[%d] = %lf \n",i,vec_out[i]);
		int end = rowPointer[nrows];
		for (int j = 0; j < end; ++j){
			printf("val[%d] = %lf \n",j,val[j]);
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
		/*
		for (int i = 0; i < nrows; ++i)
			printf("vec_out[%d] = %lf \n",i,vec_out[i]);		
		*/
	#elif USE_ESSL
					
		double orig_vec[nrows];
		for (int i = 0; i < nrows; ++i){
			orig_vec[i] = vec_out[i];
		}
		printf("checking input for dsmmx: \n");
		printf("nrows = %d \n", nrows);
		printf("nz = %d \n", nz);
		printf("ac[7][7] = %lf \n", ac[7][7]);
		printf("ka[7][7] = %d \n", ka[7][7]);
		printf("vec_in[7] = %lf \n", vec_in[7]);
		printf("vec_out[7] = %lf \n", vec_out[7]);
		printf("calling dsmmx. \n");	
		dsmmx(
			nrows, 		// number of rows in matrix
			nz, 		// maximum number of nonzero elements in each row
			ac, 		// values of the sparse matrix
			ka, 		// column indices of the sparse matrix
			nrows, 		// leading dimension of the arrays (>= nrows)
			vec_in, 
			vec_out
			);
		printf("giving matrix-vector product output. \n");		
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

void SpMatrix::vectorMultiply(std::complex<double> *vec_in, std::complex<double> *vec_out, std::complex<double> alpha, std::complex<double> beta) {

	#ifdef USE_MKL
	char mv_type = 'N';
	char matdescra[6] = {'G',' ',' ','C',' ',' '};

	// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
	mkl_zcsrmv(
		&mv_type,		// Specifies operator, transpose or not	
		&nrows,			// Number of rows in matrix
		&ncols,			// Number of cols in matrix
		&alpha,			// scalar ALPHA
		matdescra,		// Specifies type of matrix, *char
		val_c,			// nonzero elements
		colIndex,		// row indicies
		rowPointer,		// begin of col ptr
		rowPointer + 1,	// end of col ptr
		vec_in,			// input vector
		&beta,			// scalar BETA
		vec_out			// output vector
		);
		
	#else

	for (int r = 0; r < nrows; ++r){
			
				std::complex<double> temp_sum;
				temp_sum = std::complex<double>(0,0);
				int begin = rowPointer[r];
				int end = rowPointer[r+1];
				
				for (int i = begin; i < end; ++i){
					temp_sum += vec_in[colIndex[i]]*val_c[i];
				}
				
				vec_out[r] = alpha*temp_sum + beta*vec_out[r];
	}
	#endif
}


