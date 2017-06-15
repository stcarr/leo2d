/* 
 * File:   mpi_job_results.h
 * Author: Stephen
 * 
 * Created on May 24, 2017, 1:53 PM
 */

#ifndef MPI_JOB_RESULTS_H
#define MPI_JOB_RESULTS_H

#include "loc_params.h"
#include "mpi_job_params.h"

#include <string>
#include <vector>
#include <fstream>

class Mpi_job_results {
    private:
	
		// job variables (same as mpi_job_params)
		
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
		int mlmc_current_num_samples;
		
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
	
		// Result variables
		
		std::vector< std::vector<double> > cheb_coeffs;
		std::vector<double> eigenvalues;
		std::vector< std::vector<double> > eigenweights;
		std::vector< std::vector<double> > eigenvectors;
		std::vector< std::vector<double> > M_xx;
		std::vector< std::vector<double> > M_xy;
		std::vector< std::vector<double> > M_yy;
		
		std::vector<double> cpu_time;
		std::vector<std::string> cpu_time_type;
		
        Mpi_job_results();
        Mpi_job_results(const Mpi_job_results& orig);
        ~Mpi_job_results();
		
		void loadLocParams(Loc_params opts);
		void loadJobParams(Mpi_job_params opts);
		
		void saveTiming(double, std::string);
		
        void setParam(std::string, int);
		void setParam(std::string, double);
		void setParam(std::string, int*, int);
		void setParam(std::string, double*, int);
		void setParam(std::string, int*, int, int);
		void setParam(std::string, double*, int, int);
		void setParam(std::string, std::vector<double>);
		void setParam(std::string, std::vector< std::vector<double> >);
		void setParam(std::string, std::vector< std::string >);
		
		std::vector<double> getResultVec(std::string);
		std::vector< std::vector<double> > getResultMat(std::string);
		
		int getInt(std::string) const;
		double getDouble(std::string) const;
		int* getIntVec(std::string) const;
		double* getDoubleVec(std::string) const;
		int* getIntMat(std::string) const;
		double* getDoubleMat(std::string) const;
		
		void save(std::ofstream&);
		void printParams();
		void printHeader(std::ofstream&);
		
		void mlmc_load(std::string);
		void mlmc_save(std::string);
		void mlmc_average(Mpi_job_results);
		void mlmc_variance(Mpi_job_results);
		void mlmc_cluster_average(Mpi_job_results,Mpi_job_results);
		void mlmc_cluster_variance(Mpi_job_results,Mpi_job_results);
		
		void conductivtyTransform();
		
		void send(int, int);
		void sendInt(int, int, int);
		void sendDouble(double, int, int);
		void sendIntVec(int*, int, int, int);
		void sendDoubleVec(double*, int, int, int);
		void sendIntMat(int*, int, int, int, int);
		void sendDoubleMat(double*, int, int, int, int);
		
		void CPP_sendDoubleVec(std::vector<double>, int, int);
		void CPP_sendDoubleMat(std::vector< std::vector<double> >, int, int);		
		void CPP_sendStrVec(std::vector< std::string >, int, int);		
		
		int recv_spool();
		void recv(int);
		void recvInt(int, std::string);
		void recvDouble(int, std::string);
		void recvIntVec(int, std::string);
		void recvDoubleVec(int, std::string);
		void recvIntMat(int, std::string);
		void recvDoubleMat(int, std::string);
		
		void CPP_recvDoubleVec(int, std::string);
		void CPP_recvDoubleMat(int, std::string);
		void CPP_recvStrVec(int, std::string);

};

#endif /* MPI_JOB_RESULTS_H */