/*
 * File:   mlmc_handler.cpp
 * Author: Stephen
 *
 * Created on May 24, 2017, 1:53 PM
 */

#include "mlmc_handler.h"

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
	average.loadLocParams(opts);
	variance.loadLocParams(opts);
	cluster_average.loadLocParams(opts);
	cluster_variance.loadLocParams(opts);

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

			cluster_original[i].loadLocParams(opts);

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

			cluster_original[i].mlmc_load(tar_file);
		}
	}

}

void Mlmc_handler::process(Mpi_job_results results){

	//printf("Entering mlmc_handler.process() \n");

	if (d_cond > 0){
		results.conductivityTransform();
	}

	int jobID = results.getInt("jobID");
	//printf("jobID = %d \n",jobID);
	int mlmc_clusterID = results.getInt("mlmc_clusterID") ;
	//printf("mlmc_clusterID = %d \n",mlmc_clusterID);

	if (k_sampling == 0){
		// First check if this is part of a cluster
		if (mlmc_clusterID == -1){

			// Not part of a cluster from previous level, so save as stats for this level
			average.mlmc_average(results);
			variance.mlmc_variance(results);

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

				results.mlmc_save(tar_file);
			}

		} else {

			//printf("results.M_xx[0][0] = %lf \n",results.M_xx[0][0]);
			//printf("cluster_original[mlmc_clusterID].M_xx[0][0] = %lf \n",cluster_original[mlmc_clusterID-1].M_xx[0][0]);

			// Else save the cluster average

			cluster_average.mlmc_cluster_average(cluster_original[mlmc_clusterID-1],results);
			//printf("mlmc_handler.process() completed dE process \n");
			cluster_variance.mlmc_cluster_variance(cluster_original[mlmc_clusterID-1],results);
			//printf("mlmc_handler.process() completed dV process \n");
		}
		
	} else if (k_sampling == 1){
	
		// figure out how many jobIDs we have staged for k_sampling
		int k_stage_size = (int)k_staging_results.size();
		// keep track of where we are in the staging arrays
		int result_index = -1;		
		
		// if we have none, we add this result immediately
		if (k_stage_size == 0){
			std::vector<Mpi_job_results> first_k_results;
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
				std::vector<Mpi_job_results> first_k_results;
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
		
			// kind of hacky, need to find a better way to do this!
			Mpi_job_results k_results(results);
			std::vector< std::vector<double> > M_xx = k_staging_results[result_index][0].M_xx;
			std::vector< std::vector<double> > M_xy = k_staging_results[result_index][0].M_xy;
			std::vector< std::vector<double> > M_yy = k_staging_results[result_index][0].M_yy;
			
			for (int i = 1; i < num_k; ++i){
			
				for (int x = 0; x < (int)M_xx.size(); ++x){
					for (int y = 0; y < (int)M_xx[x].size(); ++y){
						M_xx[x][y] = M_xx[x][y] + (k_staging_results[result_index][i].M_xx[x][y]);
					}
				}
				
				for (int x = 0; x < (int)M_xy.size(); ++x){
					for (int y = 0; y < (int)M_xy[x].size(); ++y){
						M_xy[x][y] = M_xy[x][y] + (k_staging_results[result_index][i].M_xy[x][y]);
					}
				}
				
				for (int x = 0; x < (int)M_yy.size(); ++x){
					for (int y = 0; y < (int)M_yy[x].size(); ++y){
						M_yy[x][y] = M_yy[x][y] + (k_staging_results[result_index][i].M_yy[x][y]);
					}
				}
				
			}
		
			// normalize
			
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

			k_results.M_xx = M_xx;
			k_results.M_xy = M_xy;
			k_results.M_yy = M_yy;

			// First check if this is part of a cluster
			if (mlmc_clusterID == -1){

				// Not part of a cluster from previous level, so save as stats for this level
				average.mlmc_average(k_results);
				variance.mlmc_variance(k_results);

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

					k_results.mlmc_save(tar_file);
				}

			} else {

				//printf("results.M_xx[0][0] = %lf \n",results.M_xx[0][0]);
				//printf("cluster_original[mlmc_clusterID].M_xx[0][0] = %lf \n",cluster_original[mlmc_clusterID-1].M_xx[0][0]);

				// Else save the cluster average

				cluster_average.mlmc_cluster_average(cluster_original[mlmc_clusterID-1],k_results);
				//printf("mlmc_handler.process() completed dE process \n");
				cluster_variance.mlmc_cluster_variance(cluster_original[mlmc_clusterID-1],k_results);
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

	average.mlmc_save(tar_average_file);
	variance.mlmc_save(tar_variance_file);

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

		cluster_average.mlmc_save(tar_d_average_file);
		cluster_variance.mlmc_save(tar_d_variance_file);

	}

	//printf("Done with mlmc_handler.save() \n");

}
