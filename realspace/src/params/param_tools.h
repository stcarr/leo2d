/*
 * File:   param_tools.h
 * Author: Stephen
 *
 * Created on August 21, 2017, 11:02 AM
 */
#ifndef PARAM_TOOLS_H
#define PARAM_TOOLS_H 
 
 
#include "job_params.h"

/**
 * A namespace for all the processing methods for the Job_params class
 */
namespace Param_tools{

	// Save data from job to outFile
	void save(Job_params job, std::ofstream& outFile);
	
	// Save header information for job to outFile
	void saveHeader(Job_params job, std::ofstream& outFile);

	// Print some settings for job to terminal
	// Maybe not use, already implemented in Job_params class pretty easily
	//void printParams(Job_params job);	
	
	// Used in locality to save the time it took to get from previous to current tag
	void saveTiming(Job_params& job, double t, std::string tag);

	// Transform from chebyshev basis to energy basis
	void densityTransform(Job_params& job);
	
	// Transform from chebyshev basis to energy basis
	void conductivityTransform(Job_params& job);

	// Transform any matrix from chebyshev basis to energy basis
	void matrixResponseTransform(Job_params& job, std::string tag);
	
	// load MLMC data from file_name into job
	void mlmc_load(Job_params& job, std::string file_name);

	// save MLMC data from job into file_name
	void mlmc_save(Job_params job, std::string file_name);

	// average the sample into the total
	void mlmc_average(Job_params& total, Job_params samp);

	// average the sample variance into the total
	void mlmc_variance(Job_params& total, Job_params samp);

	// average the change in the average between two MLMC levels (orig vs cluster) into the total
	void mlmc_cluster_average(Job_params& total, Job_params samp_orig, Job_params samp_cluster);

	// average the change in the variance between two MLMC levels (orig vs cluster) into the total
	void mlmc_cluster_variance(Job_params& total, Job_params samp_orig, Job_params samp_cluster);

}

#endif /* PARAM_TOOLS_H */

