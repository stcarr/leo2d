/*
 * File:   mlmc_handler.h
 * Author: Stephen
 *
 * Created on May 24, 2017, 1:53 PM
 */

#ifndef MLMC_HANDLER_H
#define MLMC_HANDLER_H

#include "params/job_params.h"

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
    int d_kpm_dos;
    int kpm_trace;

		std::string out_root;
		std::string temp_root;
		std::string prefix;

		int mlmc_max_level;
		int mlmc_level;
		int mlmc_num_clusters;
		int mlmc_cluster_size;

		Job_params average;
		Job_params variance;

		Job_params cluster_average;
		Job_params cluster_variance;
		std::vector<Job_params> cluster_original;

		int k_sampling;
		int num_k;
		std::vector<std::vector<Job_params> > k_staging_results;
		std::vector<int> k_staging_jobID;
		std::vector<int> k_staging_k_count;

	public:

		Mlmc_handler();
		Mlmc_handler(const Mlmc_handler&);
		~Mlmc_handler();

		void setup(Job_params);
		void process(Job_params);
    void k_sampling_average(std::vector<Job_params> k_samples, Job_params& k_results);
		void save();

};

#endif /* MLMC_HANDLER_H */
