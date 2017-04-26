/* 
 * File:   dmatrix.cpp
 * Author: Stephen Carr
 * Based on MATKIT from the FILTLAN package
 * 
 * Created on April 12, 2017, 3:38 PM
 */

#include "dmatrix.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
 
#ifdef USE_MKL
	#include "mkl.h"
#elif USE_ESSL
	#include <essl.h> 
#endif
 
// constructors
// a constructor for an empty matrix
DMatrix::DMatrix() {
    nrows = 0;
    ncols = 0;
    val = NULL;
	val_c = NULL;
    nval = 0;
	type = 0;
}


// a constructor for an nr-by-nc (sparse) matrix of zeros
DMatrix::DMatrix(int nr, int nc) {
    nrows = nr;
    ncols = nc;
    val = NULL;
    nval = nr*nc;
	type = 0;
}

// a constructor for an nr-by-nc matrix
DMatrix::DMatrix(int nr, int nc, double *val0) {
    
	type = 0;


	// set dimensions
	nrows = nr;
	ncols = nc;

	val = val0;
	nval = nr*nc;

}

// Complex constructor
DMatrix::DMatrix(int nr, int nc, std::complex<double> *val_c0) {
    
	type = 1;

	// set dimensions
	nrows = nr;
	ncols = nc;

	val_c = val_c0;
	nval = nr*nc;
}


// a destructor
DMatrix::~DMatrix(){

	if (type == 0){
		delete [] val;
	}
	
	if (type == 1){
		delete [] val_c;
	}

}



void DMatrix::setup(int nr, int nc){
	
	nrows = nr;
	ncols = nc;
	nval  = nr*nc;
	
}

// Real setup
void DMatrix::setup(int nr, int nc, double *val0){   
	
	type = 0;


	// set dimensions
	nrows = nr;
	ncols = nc;

	val = val0;
	nval = nr*nc;
}

double* DMatrix::allocRealVal(){
	val = new double[nval];
	for (int i = 0; i < nval; ++i){
		val[i] = 0;
	}
	type = 0;
	return val;
}

std::complex<double>* DMatrix::allocCpxVal(){
	val_c = new std::complex<double>[nval];
	for (int i = 0; i < nval; ++i){
		val_c[i] = std::complex<double>(0.0, 0.0);
	}
	type = 1;
	return val_c;
}

// Complex setup
void DMatrix::setup(int nr, int nc, std::complex<double> *val_c0) {
    
	type = 1;

	// set dimensions
	nrows = nr;
	ncols = nc;

	val_c = val_c0;
	nval = nr*nc;
}

double* DMatrix::getValPtr(){
	return val;
}

void DMatrix::setVal(double* ptr){
	for (int i = 0; i < nval; ++i){
		val[i] = ptr[i];
	}
}

void DMatrix::squareAllVals(){
	for (int i = 0; i < nval; ++i){
		val[i] = val[i]*val[i];
	}
}


// Matrix-Vector Multiplication
void DMatrix::vectorMultiply(double *vec_in, double *vec_out, double alpha, double beta) {

	if (vec_out == NULL) {
		// allocate memory
		vec_out = new double[ncols];
		for (int i = 0; i < ncols; ++i){
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
		int v_stride = 1;
		// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
		dgemv(
			&mv_type,		// Specifies operator, transpose or not	
			&nrows,			// Number of rows in matrix
			&ncols,			// Number of cols in matrix
			&alpha,			// scalar ALPHA
			val,			// nonzero elements
			&nrows,			// leading dimension (we do not skip any elements)
			vec_in,			// input vector
			&v_stride,		// increment/iterator of vec_in (how to move in memory, starting at the pointer)
			&beta,			// scalar BETA
			vec_out,		// output vector
			&v_stride		// increment/iterator of vec_out (see above)	
			);
		/*
		for (int i = 0; i < nrows; ++i)
			printf("vec_out[%d] = %lf \n",i,vec_out[i]);		
		*/
	#else
	// Column-major matrix-vector product algorithm
	for (int r = 0; r < nrows; ++r){
		vec_out[r] = beta*vec_out[r];
	}
	
	for (int c = 0; c < ncols; ++c){
		
		for (int r = 0; r < nrows; ++r){
			vec_out[r] += alpha*val[c*nrows + r]*vec_in[c];
		}
	}
	#endif

}

void DMatrix::vectorMultiply(std::complex<double> *vec_in, std::complex<double> *vec_out, std::complex<double> alpha, std::complex<double> beta) {
	
	if (vec_out == NULL) {
		// allocate memory
		vec_out = new std::complex<double>[ncols];
		for (int i = 0; i < ncols; ++i){
			vec_out[i] = std::complex<double>(0.0,0.0);
			}
	}
	#ifdef USE_MKL
	char mv_type = 'N';
	int v_stride = 1;
	// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
	zgemv(
		&mv_type,		// Specifies operator, transpose or not	
		&nrows,			// Number of rows in matrix
		&ncols,			// Number of cols in matrix
		&alpha,			// scalar ALPHA
		val_c,			// nonzero elements
		&nrows,			// leading dimension (we do not skip any elements)
		vec_in,			// input vector
		&v_stride,		// increment/iterator of vec_in (how to move in memory, starting at the pointer)
		&beta,			// scalar BETA
		vec_out,		// output vector
		&v_stride		// increment/iterator of vec_out (see above)	
		);
		
	#else
	// Column-major matrix-vector product algorithm
	for (int r = 0; r < nrows; ++r){
		vec_out[r] = beta*vec_out[r];
	}
	
	for (int c = 0; c < ncols; ++c){
		
		for (int r = 0; r < nrows; ++r){
			vec_out[r] += alpha*val_c[c*nrows + r]*vec_in[c];
		}
	}
	#endif
}

void DMatrix::matrixMultiply(DMatrix &C, DMatrix &B, double alpha, double beta) {
	
	// A = this matrix
	// C := alpha*A*B + beta*C,

	int ncols_A = ncols;
	int nrows_A = nrows;
	
	int ncols_B = B.getNumCols();
	int nrows_B = B.getNumRows();
	
	if (ncols_A != nrows_B){
		throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
	}
	
	double* val_B;
	val_B = B.getValPtr();
	
	double* val_C;
	int nval_C;
	
	int ncols_C;
	int nrows_C;
	
	if (C.getValPtr() == NULL){
	
		ncols_C = ncols_B;
		nrows_C = nrows_A;
	
		C.setup(nrows_C,ncols_C);
		nval_C = nrows_C*ncols_C;
		val_C = C.allocRealVal();
	
	} else {
		ncols_C = C.getNumCols();
		nrows_C= C.getNumRows();
		
		if (nrows_A != ncols_C || ncols_B != nrows_C){
			throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
		}

		nval_C = ncols_C*nrows_C;
		val_C = C.getValPtr();
	}
	
	#ifdef USE_MKL
		char mm_type = 'N';
		// C := alpha*op( A )*op( B ) + beta*C,
		dgemm(
			&mm_type,		// Specifies operator on A, transpose or not	
			&mm_type,		// Specifies operator on B, transpose or not
			&nrows_A,		// Number of rows in matrix A (rows of C)
			&ncols_C,		// Number of cols in matrix C (cols of B)
			&ncols_A,		// Number of internal cols/rows of A/B (summed over)
			&alpha,			// scalar ALPHA
			val,			// elements of matrix A
			&nrows_A,		// leading dimension (we do not skip any elements)
			val_B,			// elements of matrix B
			&nrows_B,		// leading dimension (we do not skip any elements)
			&beta,			// scalar BETA
			val_C,			// output matrix
			&nrows_C		// increment/iterator of vec_out (see above)	
			);
			
	#else
		for (int c = 0; c < ncols_B; ++c){
			for (int r = 0; r < nrows_A; ++r){
				val_C[c*nrows_C + r] = beta*val_C[c*nrows_C + r];
			}
		}
		
		for (int c = 0; c < ncols_B; ++c){
			for (int c_A = 0; c_A < ncols_A; ++c_A){
				for (int r = 0; r < nrows_A; ++r){
					val_C[c*nrows_C + r] += alpha*val[c_A*nrows_A + r]*val_B[c*nrows_B + c_A];
				}
			}
		}
	#endif

}

void DMatrix::matrixMultiply(DMatrix &C, DMatrix &B, double alpha, double beta, char A_type, char B_type) {
	
	// A = this matrix
	// C := alpha*op( A )*op( B ) + beta*C,

	int ncols_A = ncols;
	int nrows_A = nrows;
	
	int ncols_B = B.getNumCols();
	int nrows_B = B.getNumRows();

	if (A_type == 'N'){
		if (B_type == 'N'){
			if (ncols_A != nrows_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		} else if (B_type == 'T'){
			if (ncols_A != ncols_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		}
	} else if (A_type == 'T'){
		if (B_type == 'N'){
			if (nrows_A != nrows_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		} else if (B_type == 'T'){
			if (nrows_A != ncols_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		}
	}	

	double* val_B;
	val_B = B.getValPtr();
	
	double* val_C;
	int nval_C;
	
	int ncols_C;
	int nrows_C;
	
	int k_internal;
	
	if (C.getValPtr() == NULL){
	
		if (B_type == 'N'){
			ncols_C = ncols_B;
		} else if (B_type == 'T'){
			ncols_C = nrows_B;
		}
		
		if (A_type == 'N'){
			nrows_C = nrows_A;
			k_internal = ncols_A;
			//printf("A = 'N' \n");
		} else if (A_type == 'T'){
			nrows_C = ncols_A;
			k_internal = nrows_A;
			//printf("A = 'T' \n");

		}
	
		C.setup(nrows_C,ncols_C);
		nval_C = nrows_C*ncols_C;
		val_C = C.allocRealVal();
	
	} else {
		ncols_C = C.getNumCols();
		nrows_C= C.getNumRows();
		
		if (A_type == 'N'){
			if (B_type == 'N'){
				if (nrows_A != nrows_C || ncols_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			} else if (B_type == 'T'){
				if (nrows_A != nrows_C || nrows_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			}
		} else if (A_type == 'T'){
			if (B_type == 'N'){
				if (ncols_A != nrows_C || ncols_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			} else if (B_type == 'T'){
				if (ncols_A != nrows_C || nrows_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			}
		}

		nval_C = ncols_C*nrows_C;
		val_C = C.getValPtr();
	}

	
	#ifdef USE_MKL
		// A_type, B_type can be 'N' or 'T'
		// 'N' For normal, 'T' for transpose

		// C := alpha*op( A )*op( B ) + beta*C,
		dgemm(
			&A_type,		// Specifies operator on A, transpose or not	
			&B_type,		// Specifies operator on B, transpose or not
			&nrows_A,		// Number of rows in matrix A (rows of C)
			&ncols_B,		// Number of cols in matrix B (cols of C)
			&ncols_A,		// Number of internal cols/rows of A/B (summed over)
			&alpha,			// scalar ALPHA
			val,			// elements of matrix A
			&nrows_A,		// leading dimension (we do not skip any elements)
			val_B,			// elements of matrix B
			&nrows_B,		// leading dimension (we do not skip any elements)
			&beta,			// scalar BETA
			val_C,			// output matrix
			&nrows_C		// increment/iterator of vec_out (see above)	
			);
			
	#else
		for (int c = 0; c < ncols_C; ++c){
			for (int r = 0; r < nrows_C; ++r){
				val_C[c*nrows_C + r] = beta*val_C[c*nrows_C + r];
			}
		}
		
		//printf("ncols_C = %d, k_internal = %d, nrows_C = %d \n", ncols_C, k_internal, nrows_C);
		for (int c = 0; c < ncols_C; ++c){
			for (int k = 0; k < k_internal; ++k){
				for (int r = 0; r < nrows_C; ++r){
					//printf("[%d, %d, %d] \n",c,k,r);
					if (A_type == 'N'){
						if (B_type == 'N'){
							// Access C(r,c) += A(r,k)*B(k,c)
							val_C[c*nrows_C + r] += alpha*val[k*nrows_A + r]*val_B[c*nrows_B + k];
						} else if (B_type == 'T'){ 
							// Access C(r,c) += A(r,k)*B(c,k)
							val_C[c*nrows_C + r] += alpha*val[k*nrows_A + r]*val_B[k*nrows_B + c];
						}
					} else if (A_type == 'T'){
						if (B_type == 'N'){
							// Access C(r,c) += A(k,r)*B(k,c)
							val_C[c*nrows_C + r] += alpha*val[r*nrows_A + k]*val_B[c*nrows_B + k];
						} else if (B_type == 'T'){ 
							// Access C(r,c) += A(k,r)*B(c,k)
							val_C[c*nrows_C + r] += alpha*val[r*nrows_A + k]*val_B[k*nrows_B + c];
						}
					} 
				}
			}
		}

	#endif

}

void DMatrix::eleMatrixMultiply(DMatrix &C, DMatrix &B, double alpha, double beta){
	
	// A = this matrix
	// C := alpha*A.*B + beta*C,

	int ncols_A = ncols;
	int nrows_A = nrows;
	
	int ncols_B = B.getNumCols();
	int nrows_B = B.getNumRows();
	
	if (ncols_A != ncols_B || nrows_A != nrows_B){
		throw std::invalid_argument("Matrix Size Mismatch for A,B dimensions in C = A.*B \n");
	}
	
	double* val_B;
	val_B = B.getValPtr();
	
	double* val_C;
	int nval_C;
	
	int ncols_C;
	int nrows_C;
	
	if (C.getValPtr() == NULL){
	
		ncols_C = ncols_A;
		nrows_C = nrows_A;
	
		C.setup(nrows_C,ncols_C);
		nval_C = nrows_C*ncols_C;
		val_C = C.allocRealVal();
	
	} else {
		ncols_C = C.getNumCols();
		nrows_C= C.getNumRows();
		
		if (nrows_C != nrows_A || ncols_C != ncols_A){
			throw std::invalid_argument("Matrix Size Mismatch for C in C = A.*B \n");
		}

		nval_C = ncols_C*nrows_C;
		val_C = C.getValPtr();
	}
	
	for (int i = 0; i < nval_C; ++i){
		val_C[i] = alpha*val[i]*val_B[i] + beta*val_C[i];
	}

}

int DMatrix::getNumRows(){
	return nrows;
}

int DMatrix::getNumCols(){
	return ncols;
}

double DMatrix::getFirstVal(){
	if (type == 0)
		return val[0];
	else
		return -999;
}
