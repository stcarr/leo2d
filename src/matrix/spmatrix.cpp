/*
 * File:   spmatrix.cpp
 * Author: Stephen Carr
 * Based on MATKIT from the FILTLAN package
 *
 * Created on May 18, 2016, 12:35 PM
 */

#include "spmatrix.h"

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



void SpMatrix::setup(int maxnnz0, int nr, int nc){

	nrows = nr;
	ncols = nc;
	maxnnz = maxnnz0;

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

void SpMatrix::debugPrint(){

	for (int r = 0; r < nrows; ++r){
		printf("      ");

		int begin = rowPointer[r];
		int end = rowPointer[r+1];

		for (int i = begin; i < end; ++i){
			if (type == 0){
				printf("%lf(%d,%d) ",val[i],r,colIndex[i]);
			} else if (type == 1){
				printf("[%lf,%lf](%d,%d) ",val_c[i].real(), val_c[i].imag(),r,colIndex[i]);
			}
		}
		printf("\n");

	}
	printf("\n");

}

int* SpMatrix::allocColIndx(){
	colIndex = new int[maxnnz];
	return colIndex;
}

int* SpMatrix::allocRowPtr(){
	rowPointer = new int[nrows+1];
	return rowPointer;
}

double* SpMatrix::allocRealVal(){
	val = new double[maxnnz];
	for (int idx = 0; idx < maxnnz; ++idx){
		val[idx] = 0.0;
	}
	type = 0;
	return val;
}

std::complex<double>* SpMatrix::allocCpxVal(){
	val_c = new std::complex<double>[maxnnz];
	for (int idx = 0; idx < maxnnz; ++idx){
		val_c[idx] = std::complex<double>(0.0,0.0);
	}
	type = 1;
	return val_c;
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
		// !! STILL BUGGED !!
		// ESSL is annoying and only used at ALCF so this development is paused for now

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

// For real matrices
void SpMatrix::vectorMultiply(std::complex<double> *vec_in, std::complex<double> *vec_out, std::complex<double> alpha, std::complex<double> beta) {

	if (vec_out == NULL) {
		// allocate memory
		vec_out = new std::complex<double>[ncols];
		for (int i = 0; i < ncols; ++i){
			vec_out[i] = std::complex<double>(0.0);
			}
	}

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
		colIndex,		// row indices
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

// Sparse Matrix - Dense Matrix Multiplication
void SpMatrix::denseMatrixMultiply(DMatrix &C, DMatrix &B, double alpha, double beta) {

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

	/*

	// !! DO NOT USE mkl_dcsrmm HERE !!
	// If one uses 0-based indexing (i.e. matdescra[4] == 'C', i.e. C,C++) then the dense matrix MUST be in row major order
	// For column-major order dense matrices, one must use 1-based indexing (matdescra[4] == 'F', i.e. Fortran)
	//
	// We want to do lots of dense x dense matrix multiplication, but only one sparse x dense matrix multiplication
	// Unfortunately, I think this means we just ignore the MKL package for this method -Stephen


	#ifdef USE_MKL
		char mm_type = 'N';
		char matdescra[6] = {'G',' ',' ','C',' ',' '};

		// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
		mkl_dcsrmm(
			&mm_type,		// Specifies operator, transpose or not
			&nrows_A,		// Number of rows in matrix A (rows of C)
			&ncols_C,		// Number of cols in matrix C (rows of B)
			&ncols_A, 		// Number of internal cols/rows of A/B (summed over)
			&alpha,			// scalar ALPHA
			matdescra,		// Specifies type of matrix, *char
			val,			// nonzero elements of matrix A
			colIndex,		// row indices
			rowPointer,		// begin of col pointer
			rowPointer + 1,	// end of col pointer
			val_B,			// pointer for values of matrix B
			&nrows_B,		// leading dimension of B (different if this is a submatrix problem)
			&beta,			// scalar BETA
			val_C,			// pointer for values of matrix C
			&nrows_C		// leading dimension of C (different if this is a submatrix problem)
			);
	#else
	*/
		for (int r = 0; r < nrows_A; ++r){
			for (int c = 0; c < ncols_B; ++c){

				double temp_sum = 0.0;

				int begin = rowPointer[r];
				int end = rowPointer[r+1];

				for (int i = begin; i < end; ++i){
					temp_sum += val[i]*val_B[c*nrows_B + colIndex[i]];
				}

				val_C[c*nrows_C + r] = alpha*temp_sum + beta*val_C[c*nrows_C + r];
			}
		}

	// #endif

}

// Sparse Matrix - Dense Matrix Multiplication
// For complex matrices
void SpMatrix::denseMatrixMultiply(DMatrix &C, DMatrix &B, std::complex<double> alpha, std::complex<double> beta) {

  // A = this matrix
	// C := alpha*A*B + beta*C,

	int ncols_A = ncols;
	int nrows_A = nrows;

	int ncols_B = B.getNumCols();
	int nrows_B = B.getNumRows();

	if (ncols_A != nrows_B){
		throw std::invalid_argument("Matrix Size Mismatch for A,B inner dimensions in C = A*B \n");
	}

	std::complex<double>* val_B;
	val_B = B.getCpxValPtr();

	std::complex<double>* val_C;
	int nval_C;

	int ncols_C;
	int nrows_C;

	if (C.getCpxValPtr() == NULL){

		ncols_C = ncols_B;
		nrows_C = nrows_A;

		C.setup(nrows_C,ncols_C,1);
		nval_C = nrows_C*ncols_C;
		val_C = C.allocCpxVal();

	} else {
		ncols_C = C.getNumCols();
		nrows_C= C.getNumRows();

		if (nrows_A != ncols_C || ncols_B != nrows_C){
			throw std::invalid_argument("Matrix Size Mismatch for C in C = A*B \n");
		}

		nval_C = ncols_C*nrows_C;
		val_C = C.getCpxValPtr();
	}

	/*

	// !! DO NOT USE mkl_zcsrmm HERE !!
	// If one uses 0-based indexing (i.e. matdescra[4] == 'C', i.e. C,C++) then the dense matrix MUST be in row major order
	// For column-major order dense matrices, one must use 1-based indexing (matdescra[4] == 'F', i.e. Fortran)
	//
	// We want to do lots of dense x dense matrix multiplication, but only one sparse x dense matrix multiplication
	// Unfortunately, I think this means we just ignore the MKL package for this method -Stephen


	#ifdef USE_MKL
		char mm_type = 'N';
		char matdescra[6] = {'G',' ',' ','C',' ',' '};

		// vec_out = alpha*A*vec_in + beta*vec_out if mv_type = 'N'
		mkl_zcsrmm(
			&mm_type,		// Specifies operator, transpose or not
			&nrows_A,		// Number of rows in matrix A (rows of C)
			&ncols_C,		// Number of cols in matrix C (rows of B)
			&ncols_A, 		// Number of internal cols/rows of A/B (summed over)
			&alpha,			// scalar ALPHA
			matdescra,		// Specifies type of matrix, *char
			val,			// nonzero elements of matrix A
			colIndex,		// row indices
			rowPointer,		// begin of col pointer
			rowPointer + 1,	// end of col pointer
			val_B,			// pointer for values of matrix B
			&nrows_B,		// leading dimension of B (different if this is a submatrix problem)
			&beta,			// scalar BETA
			val_C,			// pointer for values of matrix C
			&nrows_C		// leading dimension of C (different if this is a submatrix problem)
			);
	#else
	*/

		for (int r = 0; r < nrows_A; ++r){
			for (int c = 0; c < ncols_B; ++c){

				std::complex<double> temp_sum(0,0);

				int begin = rowPointer[r];

				int end = rowPointer[r+1];

				for (int i = begin; i < end; ++i){
					temp_sum += val_c[i]*val_B[c*nrows_B + colIndex[i]];
				}
				val_C[c*nrows_C + r] = alpha*temp_sum + beta*val_C[c*nrows_C + r];
			}
		}

	// #endif

}

void SpMatrix::denseConvert(DMatrix &Mat_in){

	// dense matrix is in Column Major order
	Mat_in.setup(nrows, ncols, type);

	if (type == 0){

		double* mat_val = Mat_in.allocRealVal();

		for (size_t r = 0; r < nrows; ++r){

			int begin = rowPointer[r];
			int end = rowPointer[r+1];

			for (size_t i = begin; i < end; ++i){
				mat_val[colIndex[i]*nrows + r] = val[i];
			}
		}

	} else if (type == 1){


		std::complex<double>* mat_val_c = Mat_in.allocCpxVal();

		for (size_t r = 0; r < nrows; ++r){
			size_t begin = rowPointer[r];
			size_t end = rowPointer[r+1];
			//printf("begin = %lu, end = %lu \n",begin,end);
			//printf("maxnnz = %d, nrows*ncols = %lu, nrows = %d, ncols = %d \n", maxnnz, (size_t)nrows*(size_t)ncols, nrows, ncols);
			for (size_t i = begin; i < end; ++i){
				//printf("i = %lu, colIndex[i] = %d, nrows = %d, r = %lu \n",i,colIndex[i],nrows,r);
				mat_val_c[((size_t)colIndex[i])*((size_t)nrows) + r] = val_c[i];
			}
		}
		//printf("done with dense val_c population \n");

	}

}

/*
void SpMatrix::denseConvert(Eigen::MatrixXd &H_in){
	// Column Major order
	for (int r = 0; r < nrows; ++r){

		int begin = rowPointer[r];
		int end = rowPointer[r+1];

		for (int i = begin; i < end; ++i){
			H_in(r,colIndex[i]) = val[i];
		}
	}

}

void SpMatrix::denseConvert(Eigen::MatrixXcd &H_in){

		for (int r = 0; r < nrows; ++r){

			int begin = rowPointer[r];
			int end = rowPointer[r+1];

			for (int i = begin; i < end; ++i){
				H_in(r,colIndex[i]) = val_c[i];
			}
		}

}
*/

void SpMatrix::eigenSolve(std::vector<double> &eigvals, DMatrix &eigvecs, char jobz, char diag_type, int il, int iu){

	DMatrix dense_mat;
	//printf("doing denseConvert() \n");
	denseConvert(dense_mat);
	//printf("entering dense_mat.eigenSolve() \n");
	dense_mat.eigenSolve(eigvals, eigvecs, jobz, diag_type, il, iu);
	//printf("done with dense_mat.eigenSolve() \n");

}

int SpMatrix::getType(){
	return type;
}

int SpMatrix::getNumRows(){
	return nrows;
}

int SpMatrix::getNumCols(){
	return ncols;
}

double SpMatrix::getFirstVal(){
	if (type == 0)
		return val[0];
	else
		return -999;
}
