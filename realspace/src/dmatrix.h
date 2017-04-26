/* 
 * File:   dmatrix.h
 * Author: Stephen Carr
 *
 * Created on April 12, 2017, 3:38 PM
 */
#include <complex>
#include <Eigen/Dense>


#define MKL_Complex16 std::complex<double>

#ifndef DMATRIX_H
#define DMATRIX_H
class DMatrix {

	protected:
		// type of matrix
		// 0: real
		// 1: complex
		int type;
		
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
		std::complex<double> *val_c;
		
		//! Maximum number of elements allowed in the allocated memory.
		// maximum number of elements allowed in the allocated memory
		int nval;
    
    public:
		// Empty Constructor
		DMatrix();
		
		// Copy Constructor
		DMatrix(const DMatrix&);
		
		// Zero constructor
		DMatrix(int,int);
		
		// Normal real constructor
		DMatrix(int,int,double*);
		
		// Normal complex constructor
		DMatrix(int,int,std::complex<double>*);
		
		double* allocRealVal();

		std::complex<double>* allocCpxVal();
		
		// before construction setup
		
		void setup(int,int);
		
		// after construction setup
		void setup(int,int,double*);
		
		// and for complex
		void setup(int,int,std::complex<double>*);
	
		// Destructor
		~DMatrix();
		
		// Value accessing
		double* getValPtr();
		void getValCopy(double*) const;
		void setVal(double* ptr);
		void squareAllVals();
		
		// Matrix-vector Multiplication
		void vectorMultiply(double*, double*, double, double);
		void vectorMultiply(std::complex<double>*, std::complex<double>*, std::complex<double>, std::complex<double>);
		
		// Matrix-matrix Multiplication
		void matrixMultiply(DMatrix&, DMatrix&, double, double);
		void matrixMultiply(DMatrix&, DMatrix&, double, double, char, char);
		
		//
		void eleMatrixMultiply(DMatrix&, DMatrix&, double, double);
		
		// Diagnostic functions
		int getNumRows() const;
		int getNumCols() const;
		double getFirstVal();

};

#endif /* DMATRIX_H */
