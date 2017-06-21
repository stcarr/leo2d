/*
 * File:   mpi_job_params.h
 * Author: Stephen
 *
 * Created on September 15, 2016, 2:52 PM
 */

#ifndef MPI_JOB_PARAMS_H
#define MPI_JOB_PARAMS_H

#include "loc_params.h"
#include "job_params.h"
#include <string>
#include <vector>
#include <fstream>

class Mpi_job_params {
    private:

		int jobID;
		int max_jobs;

		double energy_rescale;
		double energy_shift;
		double vacancy_chance;

		int solver_type;
		int observable_type;
		int solver_space;
		int diagonalize;
		int d_weights;
		int d_vecs;
		int d_cond;

		int mlmc;
		int mlmc_clusterID;
		int mlmc_level;
		int mlmc_num_clusters;
		int mlmc_cluster_size;

		int num_target_sheets;
		int* target_sheets;
		int poly_order;

		int magOn;
		int elecOn;
		double B;
		double E;

		int num_sheets;
		double* shifts;

		int num_targets;
		int* target_list;

		int num_vacancies;
		int* vacancy_list;

	public:

    Mpi_job_params();
    Mpi_job_params(const Mpi_job_params& orig);
    ~Mpi_job_params();

		void loadLocParams(Job_params opts);

    void setParam(std::string, int);
		void setParam(std::string, double);
		void setParam(std::string, int*, int);
		void setParam(std::string, double*, int);
		void setParam(std::string, int*, int, int);
		void setParam(std::string, double*, int, int);

		int getInt(std::string) const;
		double getDouble(std::string) const;
		int* getIntVec(std::string) const;
		double* getDoubleVec(std::string) const;
		int* getIntMat(std::string) const;
		double* getDoubleMat(std::string) const;

		void printParams();
		void printCheb(std::ofstream&);

		void sendParams(int, int);
		void sendInt(int, int, int);
		void sendDouble(double, int, int);
		void sendIntVec(int*, int, int, int);
		void sendDoubleVec(double*, int, int, int);
		void sendIntMat(int*, int, int, int, int);
		void sendDoubleMat(double*, int, int, int, int);

		void recvParams(int);
		void recvInt(int, std::string);
		void recvDouble(int, std::string);
		void recvIntVec(int, std::string);
		void recvDoubleVec(int, std::string);
		void recvIntMat(int, std::string);
		void recvDoubleMat(int, std::string);
};

#endif /* MPI_JOB_PARAMS_H */
