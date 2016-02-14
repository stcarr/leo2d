
/* 
 * File:   locality.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:39 PM
 */

#ifndef LOCALITY_H
#define LOCALITY_H
#include "hstruct.h"
#include "interlayer_coupling.h"

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
		int max_num_local_jobs;
		int root;
		void rootMatrixSolve(int*,double*,int*,int*,double*,int*);
		void workerMatrixSolve(int*,double*,int*,int*,double*,int*);
		static const int WORKTAG = 0;
		static const int STOPTAG = 1;
        
    public:
        Locality(std::vector<Sdata>,std::vector<double>,std::vector<double>);
        Locality(const Locality& orig);
        ~Locality();
		void setup();
		void initMPI(int, char**);
		void constructGeom();
		void constructMatrix(int*,double*,int*,int*,double*,int*);
		void plot();
		void save();
		void finMPI();

};

#endif /* LOCALITY_H */
