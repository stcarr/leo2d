
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
		int center_index;
		int max_inter_pairs;
		int max_intra_pairs;
		int root;
		int nShifts;
		int num_eigs;
		int num_samples;
		double interval_start;
		double interval_end;
		int solver_type;
		void rootEigenSolve(int*,double*,int*,int*,double*);
		void workerEigenSolve(int*,double*,int*,int*,double*);
		void rootChebSolve(int*,double*,int*,int*,double*);
		void workerChebSolve(int*,double*,int*,int*,double*);
		static const int WORKTAG = 0;
		static const int STOPTAG = 1;
        
    public:
        Locality(std::vector<Sdata>,std::vector<double>,std::vector<double>);
        Locality(const Locality& orig);
        ~Locality();
		void setup(int,int,int,double,double,int);
		void initMPI(int, char**);
		void constructGeom();
		void constructMatrix(int*,double*,int*,int*,double*);
		void plot();
		void save();
		void finMPI();

};

#endif /* LOCALITY_H */
