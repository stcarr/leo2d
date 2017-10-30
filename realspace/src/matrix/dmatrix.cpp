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
#else
	#include <Eigen/Dense>
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


DMatrix::DMatrix(const DMatrix& orig){

	nrows = orig.getNumRows();
	ncols = orig.getNumCols();
	nval = nrows*ncols;
	val = new double[nval];
	orig.getValCopy(val);

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

void DMatrix::debugPrint(){


	for (int r = 0; r < nrows; ++r){
		printf("         ");
		for (int c = 0; c < ncols; ++c){
			if (type == 0){
				printf("%lf ", val[c*nrows + r]);
			} else if (type == 1){
				printf("[%lf,%lf] ", val_c[c*nrows + r].real(),val_c[c*nrows + r].imag());
			}
		}
		printf("\n");
	}
	printf("\n");

}

void DMatrix::setup(int nr, int nc){

	nrows = nr;
	ncols = nc;
	nval  = nr*nc;

}

void DMatrix::setup(int nr, int nc, int t){

	nrows = nr;
	ncols = nc;
	nval  = nr*nc;
	type = t;

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

// Complex setup
void DMatrix::setup(int nr, int nc, std::complex<double> *val_c0) {

	type = 1;

	// set dimensions
	nrows = nr;
	ncols = nc;

	val_c = val_c0;
	nval = nr*nc;
}

std::complex<double>* DMatrix::allocCpxVal(){
	val_c = new std::complex<double>[nval];
	for (int i = 0; i < nval; ++i){
		val_c[i] = std::complex<double>(0.0, 0.0);
	}
	type = 1;
	return val_c;
}


double* DMatrix::getValPtr(){
	return val;
}

std::complex<double>* DMatrix::getCpxValPtr(){
	return val_c;
}

void DMatrix::getValCopy(double* val_copy) const{

	for (int i = 0; i < nval; ++i){
		val_copy[i] = val[i];
	}

}

void DMatrix::getValCopy(std::complex<double>* val_copy) const{

	for (int i = 0; i < nval; ++i){
		val_copy[i] = val_c[i];
	}

}

void DMatrix::getCPPValCopy(std::vector< std::vector<double> >& val_copy) const{

	if (type == 0){

		for (int i = 0; i < nrows; ++i){
			std::vector<double> temp;

			for (int j = 0; j < ncols; ++j) {
				temp.push_back(val[i + j*nrows]);
			}

			val_copy.push_back(temp);

		}
	} else if (type == 1){
		for (int i = 0; i < nrows; ++i){
			std::vector<double> temp;

			for (int j = 0; j < ncols; ++j) {
				temp.push_back(val_c[i + j*nrows].real());
			}

			val_copy.push_back(temp);

		}
	}

}
void DMatrix::getCPPValCopy(std::vector< std::vector< std::complex<double> > >& val_copy) const{

	if (type == 0){

		for (int i = 0; i < nrows; ++i){
			std::vector< std::complex<double> > temp;

			for (int j = 0; j < ncols; ++j) {
				temp.push_back(val[i + j*nrows]);
			}

			val_copy.push_back(temp);

		}
	} else if (type == 1){
		for (int i = 0; i < nrows; ++i){
			std::vector< std::complex<double> > temp;

			for (int j = 0; j < ncols; ++j) {
				temp.push_back(val_c[i + j*nrows]);
			}

			val_copy.push_back(temp);

		}
	}

}

void DMatrix::setVal(double* ptr){
	for (int i = 0; i < nval; ++i){
		val[i] = ptr[i];
	}
}

void DMatrix::setVal(std::complex<double>* ptr){
	for (int i = 0; i < nval; ++i){
		val_c[i] = ptr[i];
	}
}

void DMatrix::squareAllVals(){
	if (type == 0){
		for (int i = 0; i < nval; ++i){
			val[i] = val[i]*val[i];
		}
	} else if (type == 1){
		for (int i = 0; i < nval; ++i){
			val_c[i] = val_c[i]*(std::conj(val_c[i]));
		}
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

	matrixMultiply(C,B,alpha,beta,'N','N');

}

void DMatrix::matrixMultiply(DMatrix &C, DMatrix &B, double alpha, double beta, char A_type, char B_type) {

	// A = this matrix
	// C := alpha*op( A )*op( B ) + beta*C,
	// *_type can be:
	//								'N': Normal
	//								'T': Transpose

	if (type != 0){
		throw std::invalid_argument("Using real matrixMultiply routine for non-real matrix! \n");
	}

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
			&nrows_C,		// Number of rows in matrix A (rows of C)
			&ncols_C,		// Number of cols in matrix C (cols of B)
			&k_internal,		// Number of internal cols/rows of A/B (summed over)
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

void DMatrix::matrixMultiply(DMatrix &C, DMatrix &B, std::complex<double> alpha, std::complex<double> beta){

	matrixMultiply(C,B,alpha,beta,'N','N');

}

void DMatrix::matrixMultiply(DMatrix &C, DMatrix &B, std::complex<double> alpha, std::complex<double> beta, char A_type, char B_type) {

	// A = this matrix
	// C := alpha*op( A )*op( B ) + beta*C,
	// *_type can be:
	//								'N': Normal
	//								'T': Transpose
	//								'C': Conjugate Transpose (Hermitian Conjugate)



	if (type != 1){
		throw std::invalid_argument("Using complex matrixMultiply routine for non-complex matrix! \n");
	}

	int ncols_A = ncols;
	int nrows_A = nrows;

	int ncols_B = B.getNumCols();
	int nrows_B = B.getNumRows();

	if (A_type == 'N'){
		if (B_type == 'N'){
			if (ncols_A != nrows_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		} else if (B_type == 'T' || B_type == 'C'){
			if (ncols_A != ncols_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		}
	} else if (A_type == 'T' || A_type == 'C'){
		if (B_type == 'N'){
			if (nrows_A != nrows_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		} else if (B_type == 'T' || B_type == 'C'){
			if (nrows_A != ncols_B){
				throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
			}
		}
	}

	std::complex<double>* val_B_c;
	val_B_c = B.getCpxValPtr();

	std::complex<double>* val_C_c;
	int nval_C;

	int ncols_C;
	int nrows_C;

	int k_internal;

	if (C.getCpxValPtr() == NULL){

		if (B_type == 'N'){
			ncols_C = ncols_B;
		} else if (B_type == 'T' || B_type == 'C'){
			ncols_C = nrows_B;
		}

		if (A_type == 'N'){
			nrows_C = nrows_A;
			k_internal = ncols_A;
			//printf("A = 'N' \n");
		} else if (A_type == 'T' || A_type == 'C'){
			nrows_C = ncols_A;
			k_internal = nrows_A;
			//printf("A = 'T' \n");

		}

		C.setup(nrows_C,ncols_C,1);
		nval_C = nrows_C*ncols_C;
		val_C_c = C.allocCpxVal();

	} else {
		ncols_C = C.getNumCols();
		nrows_C= C.getNumRows();

		if (A_type == 'N'){
			if (B_type == 'N'){
				if (nrows_A != nrows_C || ncols_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			} else if (B_type == 'T' || B_type == 'C'){
				if (nrows_A != nrows_C || nrows_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			}
		} else if (A_type == 'T' || A_type == 'C'){
			if (B_type == 'N'){
				if (ncols_A != nrows_C || ncols_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			} else if (B_type == 'T' || B_type == 'C'){
				if (ncols_A != nrows_C || nrows_B != ncols_C){
					throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
				}
			}
		}

		nval_C = ncols_C*nrows_C;
		val_C_c = C.getCpxValPtr();
	}


	#ifdef USE_MKL
		// A_type, B_type can be 'N' or 'T'
		// 'N' For normal, 'T' for transpose

		// C := alpha*op( A )*op( B ) + beta*C,
		zgemm(
			&A_type,		// Specifies operator on A, transpose or not
			&B_type,		// Specifies operator on B, transpose or not
			&nrows_C,		// Number of rows in matrix A (rows of C)
			&ncols_C,		// Number of cols in matrix C (cols of B)
			&k_internal,	// Number of internal cols/rows of A/B (summed over)
			&alpha,			// scalar ALPHA
			val_c,			// elements of matrix A
			&nrows_A,		// leading dimension (we do not skip any elements)
			val_B_c,		// elements of matrix B
			&nrows_B,		// leading dimension (we do not skip any elements)
			&beta,			// scalar BETA
			val_C_c,		// output matrix
			&nrows_C		// increment/iterator of vec_out (see above)
			);

	#else
		for (int c = 0; c < ncols_C; ++c){
			for (int r = 0; r < nrows_C; ++r){
				val_C_c[c*nrows_C + r] = beta*val_C_c[c*nrows_C + r];
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
							val_C_c[c*nrows_C + r] += alpha*val_c[k*nrows_A + r]*val_B_c[c*nrows_B + k];
						} else if (B_type == 'T'){
							// Access C(r,c) += A(r,k)*B(c,k)
							val_C_c[c*nrows_C + r] += alpha*val_c[k*nrows_A + r]*val_B_c[k*nrows_B + c];
						} else if (B_type == 'C'){
							// Access C(r,c) += A(r,k)*B(c,k)
							val_C_c[c*nrows_C + r] += alpha*val_c[k*nrows_A + r]*std::conj(val_B_c[k*nrows_B + c]);
						}
					} else if (A_type == 'T'){
						if (B_type == 'N'){
							// Access C(r,c) += A(k,r)*B(k,c)
							val_C_c[c*nrows_C + r] += alpha*val_c[r*nrows_A + k]*val_B_c[c*nrows_B + k];
						} else if (B_type == 'T'){
							// Access C(r,c) += A(k,r)*B(c,k)
							val_C_c[c*nrows_C + r] += alpha*val_c[r*nrows_A + k]*val_B_c[k*nrows_B + c];
						} else if (B_type == 'C'){
							// Access C(r,c) += A(k,r)*B(c,k)
							val_C_c[c*nrows_C + r] += alpha*val_c[r*nrows_A + k]*std::conj(val_B_c[k*nrows_B + c]);
						}
					} else if (A_type == 'C'){
						if (B_type == 'N'){
							// Access C(r,c) += A(k,r)*B(k,c)
							val_C_c[c*nrows_C + r] += alpha*std::conj(val_c[r*nrows_A + k])*val_B_c[c*nrows_B + k];
						} else if (B_type == 'T'){
							// Access C(r,c) += A(k,r)*B(c,k)
							val_C_c[c*nrows_C + r] += alpha*std::conj(val_c[r*nrows_A + k])*val_B_c[k*nrows_B + c];
						} else if (B_type == 'C'){
							// Access C(r,c) += A(k,r)*B(c,k)
							val_C_c[c*nrows_C + r] += alpha*std::conj(val_c[r*nrows_A + k])*std::conj(val_B_c[k*nrows_B + c]);
						}
					}
				}
			}
		}

	#endif

}

void DMatrix::eleMatrixMultiply(DMatrix &C, DMatrix &B, double alpha, double beta, char A_type, char B_type){

	// A = this matrix
	// C := alpha*A.*B + beta*C,
	// *_type can be 'N' (normal) or 'C' (conjugate, but NOT transpose).


	int type_A = type;
	int type_B = B.getType();
	int type_C = C.getType();

	if (type_C == 0 && (type_A == 1 || type_B == 1)){
		throw std::invalid_argument("C cannot be real if A or B are complex in C = A.*B \n");
	}

	int ncols_A = ncols;
	int nrows_A = nrows;

	int ncols_B = B.getNumCols();
	int nrows_B = B.getNumRows();

	if (ncols_A != ncols_B || nrows_A != nrows_B){
		throw std::invalid_argument("Matrix Size Mismatch for A,B dimensions in C = A.*B \n");
	}

	double* val_B;
	std::complex<double>* val_B_c;
	if (type_B == 0){
		val_B = B.getValPtr();
	} else if (type_B == 1){
		val_B_c = B.getCpxValPtr();
	}

	double* val_C;
	std::complex<double>* val_C_c;

	int nval_C;

	int ncols_C;
	int nrows_C;

	if (type_C == 0){
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
	} else if (type_C == 1){
		if (C.getCpxValPtr() == NULL){

			ncols_C = ncols_A;
			nrows_C = nrows_A;

			C.setup(nrows_C,ncols_C);
			nval_C = nrows_C*ncols_C;
			val_C_c= C.allocCpxVal();

		} else {
			ncols_C = C.getNumCols();
			nrows_C= C.getNumRows();

			if (nrows_C != nrows_A || ncols_C != ncols_A){
				throw std::invalid_argument("Matrix Size Mismatch for C in C = A.*B \n");
			}

			nval_C = ncols_C*nrows_C;
			val_C_c = C.getCpxValPtr();
		}
	}

	if (type_C == 0){
		for (int i = 0; i < nval_C; ++i){
			val_C[i] = alpha*val[i]*val_B[i] + beta*val_C[i];
		}
	} else if (type_C == 1){
		if (type_B == 0){
			if (type_A == 0){
				for (int i = 0; i < nval_C; ++i){
					val_C_c[i] = alpha*val[i]*val_B[i] + beta*val_C_c[i];
				}
			} else if (type_A == 1){
				if (A_type == 'N'){
					for (int i = 0; i < nval_C; ++i){
						val_C_c[i] = alpha*val_c[i]*val_B[i] + beta*val_C_c[i];
					}
				} else if (A_type == 'C'){
					for (int i = 0; i < nval_C; ++i){
						val_C_c[i] = alpha*std::conj(val_c[i])*val_B[i] + beta*val_C_c[i];
					}
				}
			}
		} else if (type_B == 1){
			if (type_A == 0){
				if (B_type == 'N') {
					for (int i = 0; i < nval_C; ++i){
						val_C_c[i] = alpha*val[i]*val_B_c[i] + beta*val_C_c[i];
					}
				} else if (B_type == 'C'){
					for (int i = 0; i < nval_C; ++i){
						val_C_c[i] = alpha*val[i]*std::conj(val_B_c[i]) + beta*val_C_c[i];
					}
				}
			} else if (type_A == 1){
				if (B_type == 'N') {
					if (A_type == 'N') {
						for (int i = 0; i < nval_C; ++i){
							val_C_c[i] = alpha*val_c[i]*val_B_c[i] + beta*val_C_c[i];
						}
					} else if (A_type == 'C'){
						for (int i = 0; i < nval_C; ++i){
							val_C_c[i] = alpha*std::conj(val_c[i])*val_B_c[i] + beta*val_C_c[i];
						}
					}
				} else if (B_type == 'C'){
					if (A_type == 'N') {
						for (int i = 0; i < nval_C; ++i){
							val_C_c[i] = alpha*val_c[i]*std::conj(val_B_c[i]) + beta*val_C_c[i];
						}
					} else if (A_type == 'C'){
						for (int i = 0; i < nval_C; ++i){
							val_C_c[i] = alpha*std::conj(val_c[i])*std::conj(val_B_c[i]) + beta*val_C_c[i];
						}
					}
				}
			}
		}

	}

}

void DMatrix::eigenSolve(std::vector<double> &eigvals, DMatrix &eigvecs){

	if (nrows != ncols){
		printf("!! LEO2D Error !!: Trying to diagonalize non-square matrix in DMatrix.eigenSolve() \n");
		exit(1);
	}

	if (type == 0){

		eigvecs.setup(nrows, nrows);
		double* eigvecs_ptr;
		eigvecs_ptr = eigvecs.allocRealVal();

		#ifdef USE_MKL

			MKL_INT info;
			MKL_INT isuppz[2*nrows];
			int num_eigen;
			double abstol = -1;
			double vl = 0.0;
			double vu = 1.0;
			int il = 0;
			int iu = nrows;
			info = LAPACKE_dsyevr(		LAPACK_COL_MAJOR,
																'V',    			// jobz, 'N' For just vals, 'V' for vecs too
																'A',					// range, 'A' for all vals, 'V' val between  vl and vu, 'I' val indices il to iu
																'U',					// uplo, 'U' for upper triangular, 'L' for lower triangular
																nrows,				// n, order of the matrix a
																val,					// a, ptr to the matrix. Overwritten on output
																nrows,				// lda, leading dimension of the matrix, as stored in memory
																vl,						// vl, lower value for 'V' in range (second parameter)
																vu,						// vu, upper value for 'V' in range (second parameter)
																il,						// il, lower index for 'I' in range (second parameter)
																iu,						// iu, upper index for 'I' in range (second parameter)
																abstol,				// abstol, tolerance for eigensolve
																&num_eigen,		// m, output total number of eigenvalues found
																&eigvals[0],	// w, output eigenvalues
																eigvecs_ptr,	// z, output orthonormal eigenvectors
																nrows,				// ldz, leading dimension of output z (eigenvectors)
																isuppz);			// isuppz, the index of support for output z (eigenvectors)

			/* Check for convergence */
      if( info > 0 ) {
              printf( "LAPACK_dsyevr() failed to converge (eigenSolve in dmatrix.cpp). \n" );
              exit( 1 );
      }

		#else

			Eigen::MatrixXd Mat_for_eigen = Eigen::MatrixXd::Zero(nrows,nrows);

			for (int c = 0; c < ncols; ++c){
				for (int r = 0; r < nrows; ++r){
					 Mat_for_eigen(r,c) = val[c*nrows + r];
				}
			}

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Mat_for_eigen);

			Eigen::VectorXd::Map(&eigvals[0], nrows) = es.eigenvalues();
			Eigen::MatrixXd::Map(&eigvecs_ptr[0], nrows, nrows) = es.eigenvectors();

		#endif

	} else if (type == 1){

		eigvecs.setup(nrows, nrows, 1);
		std::complex<double>* eigvecs_ptr;
		eigvecs_ptr = eigvecs.allocCpxVal();

		#ifdef USE_MKL

			MKL_INT info;
			MKL_INT isuppz[2*nrows];
			int num_eigen;
			double abstol = -1;
			double vl = 0.0;
			double vu = 1.0;
			int il = 0;
			int iu = nrows;
			info = LAPACKE_zheevr(		LAPACK_COL_MAJOR,
																'V',    			// jobz, 'N' For just vals, 'V' for vecs too
																'A',					// range, 'A' for all vals, 'V' val between  vl and vu, 'I' val indices il to iu
																'U',					// uplo, 'U' for upper triangular, 'L' for lower triangular
																nrows,				// n, order of the matrix a
																val_c,				// a, ptr to the matrix. Overwritten on output
																nrows,				// lda, leading dimension of the matrix, as stored in memory
																vl,						// vl, lower value for 'V' in range (second parameter)
																vu,						// vu, upper value for 'V' in range (second parameter)
																il,						// il, lower index for 'I' in range (second parameter)
																iu,						// iu, upper index for 'I' in range (second parameter)
																abstol,				// abstol, tolerance for eigensolve
																&num_eigen,		// m, output total number of eigenvalues found
																&eigvals[0],	// w, output eigenvalues
																eigvecs_ptr,	// z, output orthonormal eigenvectors
																nrows,				// ldz, leading dimension of output z (eigenvectors)
																isuppz);			// isuppz, the index of support for output z (eigenvectors)

			/* Check for convergence */
			if( info > 0 ) {
							printf( "LAPACK_dsyevr() failed to converge (eigenSolve in dmatrix.cpp). \n" );
							exit( 1 );
			}

		#else

			Eigen::MatrixXcd Mat_for_eigen = Eigen::MatrixXcd::Zero(nrows,nrows);

			for (int c = 0; c < ncols; ++c){
				for (int r = 0; r < nrows; ++r){
					 Mat_for_eigen(r,c) = val_c[c*nrows + r];
				}
			}

			Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(Mat_for_eigen);

			Eigen::VectorXd::Map(&eigvals[0], nrows) = es.eigenvalues();
			Eigen::MatrixXcd::Map(&eigvecs_ptr[0], nrows, nrows) = es.eigenvectors();

		#endif
	}

}

int DMatrix::getNumRows() const{
	return nrows;
}

int DMatrix::getNumCols() const{
	return ncols;
}

int DMatrix::getType() const{
	return type;
}

double DMatrix::getFirstVal(){
	if (type == 0)
		return val[0];
	else
		return -9999;
}
