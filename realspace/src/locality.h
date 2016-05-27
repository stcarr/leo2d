
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
#include "spmatrix.h"
#include <time.h>
class Locality {
    private:
	
		// name for file outputs
		std::string job_name;
		
		// global MPI variables
		int size;
		int rank;
		int print_rank;
		
		// Geometery information
		std::vector<Sdata> sdata;
		std::vector<double> heights;
		std::vector<double> angles;
		int max_index;
		int center_index;
		int max_inter_pairs;
		int max_intra_pairs;
		int intra_searchsize;
		int inter_searchsize;
		
		// Timing information
		time_t constructStart;
		time_t constructEnd;
		std::vector<time_t> solverTimes;
		time_t solveStart;
		time_t solveEnd;

		// Solver settings
		int root;
		int nShifts;
		int num_eigs;
		int num_samples;
		double interval_start;
		double interval_end;
		double energy_rescale;
		double energy_shift;
		double cheb_width;
		int poly_order;
		int solver_type;
		
		// Private solver methods
		void rootChebSolve(int*,double*,int*,int*,double*);
		void workerChebSolve(int*,double*,int*,int*,double*);
		
		// MPI Communication flags
		static const int WORKTAG = 1;
		static const int STOPTAG = 0;
        
    public:
		// Main constructor
        Locality(std::vector<Sdata>,std::vector<double>,std::vector<double>);
		
		// Unimplemented copy and delete methods
        Locality(const Locality& orig);
        ~Locality();
		
		// Set the solver information
		void setup(std::string,int,int,int,double,double,double,double,double,int,int,int,int);
		
		// Starts MPI
		void initMPI(int, char**);
		
		// Root determines matrix structure and broadcasts
		void constructGeom();
		
		// Construct and solve the tight binding problem. Also saves output files
		void constructMatrix(int*,double*,int*,int*,double*);
		
		// Not implemented
		void plot();
		
		// Prints out timing information
		void save();
		
		// Ends MPI and finishes job
		void finMPI();

};

#endif /* LOCALITY_H */
