/*
 * File:   mlmc_handler.cpp
 * Author: Stephen
 *
 * Created on May 24, 2017, 1:53 PM
 */

#include "mlmc_handler.h"
#include "params/param_tools.h"

#include <stdio.h>
#include <sstream>
#include <stdlib.h>
#include <math.h>

Mlmc_handler::Mlmc_handler() {

}

Mlmc_handler::Mlmc_handler(const Mlmc_handler& orig){

}

Mlmc_handler::~Mlmc_handler(){

}

void Mlmc_handler::setup(Job_params opts_in){

	opts = opts_in;

	std::vector< std::vector<double> > blank_double_mat;

	average = Job_params(opts);
	average.setParam("mlmc_current_num_samples",1);
	average.setParam("M_xx",blank_double_mat);
	average.setParam("M_yy",blank_double_mat);
	average.setParam("M_xy",blank_double_mat);

	variance = Job_params(opts);
	variance.setParam("mlmc_current_num_samples",1);
	variance.setParam("M_xx",blank_double_mat);
	variance.setParam("M_yy",blank_double_mat);
	variance.setParam("M_xy",blank_double_mat);

	cluster_average = Job_params(opts);
	cluster_average.setParam("mlmc_current_num_samples",1);
	cluster_average.setParam("M_xx",blank_double_mat);
	cluster_average.setParam("M_yy",blank_double_mat);
	cluster_average.setParam("M_xy",blank_double_mat);

	cluster_variance = Job_params(opts);
	cluster_variance.setParam("mlmc_current_num_samples",1);
	cluster_variance.setParam("M_xx",blank_double_mat);
	cluster_variance.setParam("M_yy",blank_double_mat);
	cluster_variance.setParam("M_xy",blank_double_mat);

	poly_order = opts.getInt("poly_order");
	energy_rescale = opts.getDouble("energy_rescale");
	energy_shift = opts.getDouble("energy_shift");
	d_cond = opts.getInt("d_cond");

	k_sampling = opts.getInt("k_sampling");
	if (k_sampling == 1){
		num_k = opts.getInt("num_k1")*opts.getInt("num_k2");
	}

	out_root = opts.getString("mlmc_out_root");
	temp_root = opts.getString("mlmc_temp_root");
	prefix = opts.getString("job_name");

	mlmc_level = opts.getInt("mlmc_level");
	mlmc_max_level = opts.getInt("mlmc_max_level");

	if (mlmc_level < mlmc_max_level){

		mlmc_num_clusters = opts.getInt("mlmc_num_clusters");
		mlmc_cluster_size = opts.getInt("mlmc_cluster_size");

		cluster_original.resize(mlmc_num_clusters);

		for (int i = 0; i < mlmc_num_clusters; ++i){

			cluster_original[i] = Job_params(opts);

			std::string tar_file; // temp_root,prefix,mlmc_level+1,(i+1)
			tar_file = temp_root;
			tar_file.append("/");
			tar_file.append(prefix);
			tar_file.append("_L");
			std::ostringstream o1;
			o1 << (mlmc_level+1);
			tar_file.append(o1.str());
			tar_file.append("_J");
			std::ostringstream o2;
			o2 << (i+1);
			tar_file.append(o2.str());
			tar_file.append("_");

			Param_tools::mlmc_load(cluster_original[i],tar_file);
		}
	}

}

void Mlmc_handler::process(Job_params results){

	//printf("Entering mlmc_handler.process() \n");

	if (d_cond > 0){
		Param_tools::conductivityTransform(results);
	}

	int jobID = results.getInt("mlmcID");
	//printf("jobID = %d \n",jobID);
	int mlmc_clusterID = results.getInt("mlmc_clusterID") ;
	//printf("mlmc_clusterID = %d \n",mlmc_clusterID);

	if (k_sampling == 0){
		// First check if this is part of a cluster
		if (mlmc_clusterID == -1){

			// Not part of a cluster from previous level, so save as stats for this level
			Param_tools::mlmc_average(average, results);
			Param_tools::mlmc_variance(variance, results);

			//printf("mlmc_handler.process() completed E/V process \n");

			// Also save exact matrix for later cluster processing
			if (mlmc_level > 1){

				std::string tar_file; // temp_root,prefix,mlmc_level,jobID
				tar_file = temp_root;
				tar_file.append("/");
				tar_file.append(prefix);
				tar_file.append("_L");
				std::ostringstream o1;
				o1 << mlmc_level;
				tar_file.append(o1.str());
				tar_file.append("_J");
				std::ostringstream o2;
				o2 << jobID;
				tar_file.append(o2.str());
				tar_file.append("_");

				Param_tools::mlmc_save(results, tar_file);
			}

		} else {

			//printf("results.M_xx[0][0] = %lf \n",results.M_xx[0][0]);
			//printf("cluster_original[mlmc_clusterID].M_xx[0][0] = %lf \n",cluster_original[mlmc_clusterID-1].M_xx[0][0]);

			// Else save the cluster average

			Param_tools::mlmc_cluster_average(cluster_average, cluster_original[mlmc_clusterID-1], results);
			//printf("mlmc_handler.process() completed dE process \n");
			Param_tools::mlmc_cluster_variance(cluster_variance, cluster_original[mlmc_clusterID-1],results);
			//printf("mlmc_handler.process() completed dV process \n");
		}

	} else if (k_sampling == 1){

		// figure out how many jobIDs we have staged for k_sampling
		int k_stage_size = (int)k_staging_results.size();
		// keep track of where we are in the staging arrays
		int result_index = -1;

		// if we have none, we add this result immediately
		if (k_stage_size == 0){
			std::vector<Job_params> first_k_results;
			first_k_results.push_back(results);
			k_staging_results.push_back(first_k_results);
			k_staging_jobID.push_back(jobID);
			k_staging_k_count.push_back(1);
			result_index = 0;
		} else {

			// if not, we look if we have a matching jobID already staged
			for (int i = 0; i < k_stage_size; ++i){
				if (jobID == k_staging_jobID[i]){
					result_index = i;
				}
			}

			// if not, we make a new one
			if (result_index == -1){
				std::vector<Job_params> first_k_results;
				first_k_results.push_back(results);
				k_staging_results.push_back(first_k_results);
				k_staging_jobID.push_back(jobID);
				k_staging_k_count.push_back(1);
				result_index = k_stage_size;
			} else {
				// otherwise, we just add it it in and increment the count
				k_staging_results[result_index].push_back(results);
				k_staging_k_count[result_index] = k_staging_k_count[result_index] + 1;
			}

		}

		// check if we have all k samples for this jobID
		// if so, we create a new result which is the average of the saved ones and process it

		if(k_staging_k_count[result_index] == num_k){

			std::vector< std::vector<double> > M_xx;
			std::vector< std::vector<double> > M_xy;
			std::vector< std::vector<double> > M_yy;

			// kind of hacky, need to find a better way to do this!
			Job_params k_results(results);
			if (d_cond > 0){
				M_xx = k_staging_results[result_index][0].getDoubleMat("M_xx");
				if (d_cond > 1){
					M_xy = k_staging_results[result_index][0].getDoubleMat("M_xy");
					M_yy = k_staging_results[result_index][0].getDoubleMat("M_yy");
				}
			}

			for (int i = 1; i < num_k; ++i){

				for (int x = 0; x < (int)M_xx.size(); ++x){
					for (int y = 0; y < (int)M_xx[x].size(); ++y){
						M_xx[x][y] = M_xx[x][y] + (k_staging_results[result_index][i].getDoubleMat("M_xx")[x][y]);
					}
				}

				for (int x = 0; x < (int)M_xy.size(); ++x){
					for (int y = 0; y < (int)M_xy[x].size(); ++y){
						M_xy[x][y] = M_xy[x][y] + (k_staging_results[result_index][i].getDoubleMat("M_xx")[x][y]);
					}
				}

				for (int x = 0; x < (int)M_yy.size(); ++x){
					for (int y = 0; y < (int)M_yy[x].size(); ++y){
						M_yy[x][y] = M_yy[x][y] + (k_staging_results[result_index][i].getDoubleMat("M_xx")[x][y]);
					}
				}

			}

			// normalize
			///*
			for (int x = 0; x < (int)M_xx.size(); ++x){
				for (int y = 0; y < (int)M_xx[x].size(); ++y){
					M_xx[x][y] = M_xx[x][y]/(double)num_k;
				}
			}

			for (int x = 0; x < (int)M_xy.size(); ++x){
				for (int y = 0; y < (int)M_xy[x].size(); ++y){
					M_xy[x][y] = M_xy[x][y]/(double)num_k;
				}
			}

			for (int x = 0; x < (int)M_yy.size(); ++x){
				for (int y = 0; y < (int)M_yy[x].size(); ++y){
					M_yy[x][y] = M_yy[x][y]/(double)num_k;
				}
			}
			//*/
			if (d_cond > 0){
				k_results.setParam("M_xx", M_xx);
				if (d_cond > 1){
					k_results.setParam("M_xy", M_xy);
					k_results.setParam("M_yy", M_yy);
					}
			}
			// First check if this is part of a cluster
			if (mlmc_clusterID == -1){

				// Not part of a cluster from previous level, so save as stats for this level
				Param_tools::mlmc_average(average, k_results);
				Param_tools::mlmc_variance(variance, k_results);

				//printf("mlmc_handler.process() completed E/V process \n");

				// Also save exact matrix for later cluster processing
				if (mlmc_level > 1){

					std::string tar_file; // temp_root,prefix,mlmc_level,jobID
					tar_file = temp_root;
					tar_file.append("/");
					tar_file.append(prefix);
					tar_file.append("_L");
					std::ostringstream o1;
					o1 << mlmc_level;
					tar_file.append(o1.str());
					tar_file.append("_J");
					std::ostringstream o2;
					o2 << jobID;
					tar_file.append(o2.str());
					tar_file.append("_");

					Param_tools::mlmc_save(k_results, tar_file);
				}

			} else {

				//printf("results.M_xx[0][0] = %lf \n",results.M_xx[0][0]);
				//printf("cluster_original[mlmc_clusterID].M_xx[0][0] = %lf \n",cluster_original[mlmc_clusterID-1].M_xx[0][0]);

				// Else save the cluster average

				Param_tools::mlmc_cluster_average(cluster_average, cluster_original[mlmc_clusterID-1], k_results);
				//printf("mlmc_handler.process() completed dE process \n");
				Param_tools::mlmc_cluster_variance(cluster_variance, cluster_original[mlmc_clusterID-1], k_results);
				//printf("mlmc_handler.process() completed dV process \n");
			}

			// and we remove those results from the k_staging list
			k_staging_results.erase(k_staging_results.begin()+result_index);
			k_staging_jobID.erase(k_staging_jobID.begin()+result_index);
			k_staging_k_count.erase(k_staging_k_count.begin()+result_index);

		}

	}

}

void Mlmc_handler::save(){

	//printf("Entering mlmc_handler.save() \n");
	std::string tar_header_file; // out_root,prefix,mlmc_level

	// Print metadata about this mlmc level
	//
	//

	// Now save average and variance files

	std::ostringstream o_lvl;
	o_lvl << mlmc_level;

	std::string tar_average_file; // out_root,prefix,mlmc_level
	tar_average_file = out_root;
	tar_average_file.append("/");
	tar_average_file.append(prefix);
	tar_average_file.append("_L");
	tar_average_file.append(o_lvl.str());
	tar_average_file.append("_E_");

	std::string tar_variance_file; // out_root,prefix,mlmc_level
	tar_variance_file = out_root;
	tar_variance_file.append("/");
	tar_variance_file.append(prefix);
	tar_variance_file.append("_L");
	tar_variance_file.append(o_lvl.str());
	tar_variance_file.append("_V_");

	//printf("Saving %s and %s \n",tar_average_file.c_str(), tar_variance_file.c_str());

	Param_tools::mlmc_save(average, tar_average_file);
	Param_tools::mlmc_save(variance, tar_variance_file);

	std::string tar_E_file;
	tar_E_file = out_root;
	tar_E_file.append("/");
	tar_E_file.append(prefix);
	tar_E_file.append("_L");
	tar_E_file.append(o_lvl.str());
	tar_E_file.append("_ENERGIES.bin");

	double E[poly_order];
	for (int i = 0; i < poly_order; ++i){
		E[i] = energy_shift + energy_rescale*cos((i*1.0 + 0.5)*M_PI/poly_order);
	}

	std::ofstream file_E(tar_E_file.c_str(), std::ofstream::binary);
	char* buffer = (char*)E;
	file_E.write(buffer, poly_order*(sizeof(double)/sizeof(char)));


	// If we aren't on the top level also save cluster average and variance

	if (mlmc_level < mlmc_max_level){

		std::string tar_d_average_file; // out_root,prefix,mlmc_level
		tar_d_average_file = out_root;
		tar_d_average_file.append("/");
		tar_d_average_file.append(prefix);
		tar_d_average_file.append("_L");
		tar_d_average_file.append(o_lvl.str());
		tar_d_average_file.append("_DE_");

		std::string tar_d_variance_file; // out_root,prefix,mlmc_level
		tar_d_variance_file = out_root;
		tar_d_variance_file.append("/");
		tar_d_variance_file.append(prefix);
		tar_d_variance_file.append("_L");
		tar_d_variance_file.append(o_lvl.str());
		tar_d_variance_file.append("_DV_");

		//printf("Saving %s and %s \n",tar_d_average_file.c_str(),tar_d_variance_file.c_str());

		Param_tools::mlmc_save(cluster_average, tar_d_average_file);
		Param_tools::mlmc_save(cluster_variance, tar_d_variance_file);

	}

	//printf("Done with mlmc_handler.save() \n");

}
