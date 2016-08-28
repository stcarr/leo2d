
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
#include "loc_params.h"

#include <time.h>
#include <math.h>  
class Locality {
    private:
	
		// name for file outputs
		std::string job_name;
		
		// global MPI variables
		int size;
		int rank;
		int print_rank;
		
		// Geometry information
		std::vector<Sdata> sdata;
		std::vector<double> heights;
		std::vector<double> angles;
		int max_index;
		int* center_index;
		int max_inter_pairs;
		int max_intra_pairs;
		int intra_searchsize;
		int inter_searchsize;
		
		int magOn;
		double B;
		int elecOn;
		double E;
		
		int num_target_sheets;
		std::vector<int> target_sheets;
		double vacancy_chance;
		
		// Timing information
		time_t constructStart;
		time_t constructEnd;
		std::vector<time_t> solverTimes;
		time_t solveStart;
		time_t solveEnd;

		// Solver settings
		int root;
		int nShifts;
		double energy_rescale;
		double energy_shift;
		int poly_order;
		int solver_type;
		
		// Private solver methods
		void sendRootWork(int, int, int,std::vector< std::vector<double> >, std::vector< std::vector<int> >);
		std::vector<std::vector<int> > v_work;
		
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
		void setup(Loc_params);
		
		// Starts MPI
		void initMPI(int, char**);
		
		// Root determines matrix structure and broadcasts
		void constructGeom();
		
		// Construct and solve the tight binding problem. Also saves output files
		void constructMatrix(int*,double*,int*,int*,double*);
		
		// Calculates Peierls phase between two hopping sites;
		double peierlsPhase(double, double, double, double, double);
		
		// Calculates the on-site energy in the presence of a gated electric field
		double onSiteE(double, double, double);
		
		// Prints out timing information
		void save();
		
		// Ends MPI and finishes job
		void finMPI();

};

#endif /* LOCALITY_H */
