/*
 * File:   dmatrix.h
 * Author: Stephen Carr
 *
 * Created on April 12, 2017, 3:38 PM
 */
#include <complex>
#include <vector>

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
		size_t nval;

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
		void setup(int,int,int);

		// after construction setup
		void setup(int,int,double*);

		// and for complex
		void setup(int,int,std::complex<double>*);

		// Destructor
		~DMatrix();

		// Debug call
		void debugPrint();

		// Value accessing
		double* getValPtr();
		std::complex<double>* getCpxValPtr();

		void getValCopy(double*) const;
		void getValCopy(std::complex<double>*) const;
		void getCPPValCopy(std::vector< std::vector<              double  > >& val_copy) const;
		void getCPPValCopy(std::vector< std::vector< std::complex<double> > >& val_copy) const;


		void setVal(double*);
		void setVal(std::complex<double>* ptr);

		void squareAllVals();

		// Matrix-vector Multiplication
		void vectorMultiply(double*, double*, double, double);
		void vectorMultiply(std::complex<double>*, std::complex<double>*, std::complex<double>, std::complex<double>);

		// Matrix-matrix Multiplication
		void matrixMultiply(DMatrix&, DMatrix&, double, double);
		void matrixMultiply(DMatrix&, DMatrix&, double, double, char, char);
		void matrixMultiply(DMatrix&, DMatrix&, std::complex<double>, std::complex<double>);
		void matrixMultiply(DMatrix&, DMatrix&, std::complex<double>, std::complex<double>, char, char);

		//
		void eleMatrixMultiply(DMatrix&, DMatrix&, double, double, char A_type, char B_type);

		// Eigensolvers
		void eigenSolve(std::vector<double> &eigvals, DMatrix &eigvecs, char jobz, char diag_type, int il, int iu);

		// Diagnostic functions
		int getNumRows() const;
		int getNumCols() const;
		int getType() const;
		double getFirstVal();

};

#endif /* DMATRIX_H */
