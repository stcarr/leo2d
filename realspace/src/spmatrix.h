/* 
 * File:   locality.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:39 PM
 */

#ifndef SPMATRIX_H
#define SPMATRIX_H
class SpMatrix {

	protected:
		// number of rows and number of columns of *this
		//! Number of rows of <em>*this</em>.
		// number of rows of *this
		int nrows;
		//! Number of columns of <em>*this</em>.
		// number of columns of *this
		int ncols;

		// CSC format
		//! A storage for storing values of nonzero elements of <em>*this</em>.
		// a storage for storing values of nonzero elements of *this
		double* val;
		//! An array for column indices of <em>*this</em>.
		// an array for column indices of *this
		int *colIndex;
		//! An array for row pointers of <em>*this</em>.
		// an array for row pointers of *this
		int *rowPointer;
		//! Maximum number of nonzero elements allowed in the allocated memory.
		// maximum number of nonzero elements allowed in the allocated memory
		int maxnnz;
		
		double** ac;
		int** ka;
		int nz;
    
    public:
		// Empty Constructor
		SpMatrix();
		
		// Zero constructor
		SpMatrix(int,int);
		
		// Normal constructor
		SpMatrix(int,int,double*,int*,int*,int);
	
		// Destructor
		~SpMatrix();
		
		// Matrix-vector Multiplication
		void vectorMultiply(double*, double*, double, double);

};

#endif /* SPMATRIX_H */
