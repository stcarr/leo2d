/*
 * File:   mlmc_handler.h
 * Author: Stephen
 *
 * Created on May 24, 2017, 1:53 PM
 */

#ifndef MLMC_HANDLER_H
#define MLMC_HANDLER_H

#include "job_params.h"
#include "mpi_job_results.h"

#include <string>
#include <vector>
#include <fstream>

class Mlmc_handler {
    private:

		Job_params opts;

		int poly_order;
		double energy_rescale;
		double energy_shift;
		int d_cond;

		std::string out_root;
		std::string temp_root;
		std::string prefix;

		int mlmc_max_level;
		int mlmc_level;
		int mlmc_num_clusters;
		int mlmc_cluster_size;

		Mpi_job_results average;
		Mpi_job_results variance;

		Mpi_job_results cluster_average;
		Mpi_job_results cluster_variance;
		std::vector<Mpi_job_results> cluster_original;

	public:

		Mlmc_handler();
		Mlmc_handler(const Mlmc_handler&);
		~Mlmc_handler();

		void setup(Job_params);
		void process(Mpi_job_results);
		void save();

};

#endif /* MLMC_HANDLER_H */
