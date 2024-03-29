
/*
 * File:   locality.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:39 PM
 */

#ifndef LOCALITY_SERIAL_H
#define LOCALITY_SERIAL_H

#include "geom/hstruct.h"
#include "geom/strain.h"
#include "momentum/momentum_coupling.h"
#include "matrix/spmatrix.h"
#include "matrix/dmatrix.h"
#include "params/job_params.h"
#include "mlmc/mlmc_handler.h"

#include "fftw3.h"

#include <time.h>
#include <math.h>
class Locality_serial {
    private:

  		// name for file outputs
  		std::string job_name;

  		// global MPI variables
  		int size;
  		//int rank;
  		int print_rank;

  		// Geometry information
  		std::vector<Sdata> sdata;
  		std::vector<double> heights;
  		std::vector<double> angles;
  		int max_index;
  		int* center_index;
  		int max_inter_pairs;
  		int max_intra_pairs;

  		/*
  		int intra_searchsize;
  		int inter_searchsize;
  		*/

  		Job_params opts;

      // saved matrix for matrix-only subroutine
      DMatrix saved_matrix;

      StrainCalc strainInfo;

  		/*
  		int magOn;
  		double B;
  		int elecOn;
  		double E;

  		int num_target_sheets;
  		std::vector<int> target_sheets;
  		double vacancy_chance;
  		*/

  		// Timing information
  		time_t constructStart;
  		time_t constructEnd;
  		std::vector<time_t> solverTimes;
  		time_t solveStart;
  		time_t solveEnd;

  		// Solver settings
  		int root;

      // Momentum FFT coupling data
      Momentum_coupling fftw_inter;

  		/*
  		int nShifts;
  		double energy_rescale;
  		double energy_shift;
  		int poly_order;
  		int solver_type;
  		*/

  		// Private solver methods
  		void sendRootWork(Job_params, int);
  		std::vector<std::vector<int> > v_work;
  		std::vector<std::vector<int> > target_indices;

  		void recursiveShiftCalc(std::vector<Job_params>&, std::vector< std::vector<double> >, int, int, int, int, int, std::vector<int>, std::vector< std::vector<int> >);

  		void rootChebSolve(int* index_to_grid, double* index_to_pos,
						int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs,
						int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > intra_sc_vecs,
						std::vector< std::vector<int> > v_work, std::vector< std::vector<int> > target_indices,
            std::vector< std::vector<double> > shift_configs);


		Job_params workerChebSolve(int* index_to_grid, double* index_to_pos,
							int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs,
							int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > intra_sc_vecs,
							std::vector< std::vector<double> > shift_configs, Job_params jobIn);

  		void getVacanciesFromFile(std::vector<std::vector<int> >&, std::vector<std::vector<int> >&, Job_params);


		double getLocalTheta(int k_i, int* index_to_grid, double* index_to_pos, int max);
		int isNearestNeighbor(int i1, int j1, int o1, int s1, int i2, int j2, int o2, int s2);


	    std::vector< std::vector<double> > getReciprocal(std::vector< std::vector<double> >);
  		std::vector< std::vector<double> > getReciprocal(int);
  		double crossProd(std::vector<double> x, std::vector<double> y, int dim);
  		// void writeBufferToFile(double*, int, std::string);

  		// MPI Communication flags
  		static const int WORKTAG = 1;
  		static const int STOPTAG = 0;

    public:
  		// Main constructor
      Locality_serial(std::vector<Sdata>,std::vector<double>,std::vector<double>);

  		// Unimplemented copy and delete methods
      Locality_serial(const Locality_serial& orig);
      ~Locality_serial();

  		// Set the solver information
  		void setup(Job_params);

  		// Root determines matrix structure and broadcasts
  		void constructGeom();

  		// Construct and solve the tight binding problem. Also saves output files
  		void constructMatrix(int* index_to_grid,double* index_to_pos, int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs,
					int* intra_pairs, double* intra_pairs_t, std::vector< std::vector<int> > intra_sc_vecs,
					std::vector< std::vector<double> > shift_configs, std::vector< std::vector<int> > v_work,
					std::vector< std::vector<int> > target_indices);

  		// Updates the index_to_pos array for a specific job's orbital positions (i.e. with shift and strain)
  		void setConfigPositions(double*, double*, int*, std::vector< std::vector<double> >&, std::vector< std::vector<double> >&,  std::vector< std::vector< std::vector<double> > >&, Job_params);

      // Returns the real-space dispalcement given a certain shift configuration of an atom in sheet s
      std::vector<double> getConfigDisp(std::vector<double> config_in, int s);

  		// Creates and returns real SpMatrix objects and an array of target vectors for Electron-Electron Correlation
  		void generateRealH(SpMatrix &H, SpMatrix &dxH, SpMatrix &dyH, double* alpha_0_x_arr, double* alpha_0_y_arr, Job_params jobIn, int* index_to_grid, double* i2pos,
			int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs, int* intra_pairs, double* intra_pairs_t,
			std::vector< std::vector<int> > intra_sc_vecs, std::vector< std::vector< std::vector<double> > > strain, std::vector<int> current_index_reduction, int local_max_index);

  		// Creates and returns complex SpMatrix objects and an array of target vectors for Electron-Electron Correlation
  		void generateCpxH(SpMatrix &H, SpMatrix &dxH, SpMatrix &dyH, SpMatrix &cd_plus, SpMatrix &cd_minus, SpMatrix &dH_0_minus, SpMatrix &dH_0_plus,
			double* alpha_0_x_arr, double* alpha_0_y_arr, Job_params jobIn, int* index_to_grid, double* i2pos,
			int* inter_pairs, std::vector< std::vector<int> > inter_sc_vecs, int* intra_pairs, double* intra_pairs_t,
			std::vector< std::vector<int> > intra_sc_vecs, std::vector< std::vector< std::vector<double> > > strain, std::vector<int> current_index_reduction, int local_max_index);

  		// Creates and returns SpMatrix object representing H for a specific Momentum-space job
  		void generateMomH(SpMatrix&, Job_params, int*, double*, int*, int*, double*, std::vector<int>, int);

  		// Computes Local DOS using Chebyshev KPM methods
  		void computeDosKPM(std::vector< std::vector<double> >&,SpMatrix&, Job_params,std::vector<int>,int);

  		// Computes Local Electron-Electron Correlation using 2D Chebyshev KPM methods
  		void computeCondKPM(double*, SpMatrix&, SpMatrix&, Job_params, std::vector<int>, int, double*);

  		// Direct solver using Eigen package
  		void computeEigen(std::vector<double>&, DMatrix&, DMatrix &kpm_dos, DMatrix&, DMatrix&, DMatrix&, SpMatrix&, SpMatrix&, SpMatrix&, Job_params, std::vector<int>, int);
  		void computeEigenComplex(std::vector<double>&, DMatrix&, DMatrix &kpm_dos, DMatrix&, DMatrix&, DMatrix&, SpMatrix&, SpMatrix&, SpMatrix&, Job_params, std::vector<int>, int);

		// Computes the material's response to matIn saved via <Cheb_1|matOut|Cheb_2>
		// Use param_tools conductivityTransform method to transform to <E_1|matOut|E_2>
		  void computeMatrixResponse(Job_params jobIn, std::vector<double> eigvecs, DMatrix& eigvals, SpMatrix& matIn, DMatrix& matOut);

  		// Calculates Peierls phase between two hopping sites;
  		double peierlsPhase(double, double, double, double, double);

  		// Calculates the on-site energy in the presence of a gated electric field
  		double onSiteE(double, double, double, double);

      // returns saved matrix
      DMatrix getSavedMatrix();

  		// Prints out timing information
  		void save();

  		// Prints out job-averaged timing information
  		void printTiming(std::vector<Job_params>);


};

#endif /* LOCALITY_H */
