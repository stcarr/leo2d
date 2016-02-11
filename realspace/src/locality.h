
/* 
 * File:   locality.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:39 PM
 */

#ifndef LOCALITY_H
#define LOCALITY_H
#include "hstruct.h"

class Locality {
    private:
		int size;
		int rank;
		int print_rank;
		std::vector<Sdata> sdata;
		std::vector<double> heights;
		std::vector<double> angles;
		int max_index;
		int max_inter_pairs;
		int max_intra_pairs;
		int num_eigs;
		int root;
		int** index_to_grid;
		int** inter_pairs;
		double** intra_pairs;
		void rootMatrixSolve();
		void workerMatrixSolve();
		static const int WORKTAG = 0;
		static const int STOPTAG = 1;
        
    public:
        Locality(std::vector<Sdata>,std::vector<double>,std::vector<double>);
        Locality(const Locality& orig);
        ~Locality();
		void setup();
		void initMPI(int, char**);
		void constructGeom();
		void constructMatrix();
		void plot();
		void save();
		void finMPI();

};

#endif /* LOCALITY_H */