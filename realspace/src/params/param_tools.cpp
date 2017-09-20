/*
 * File:   param_tools.cpp
 * Author: Stephen
 *
 * Created on August 21, 2017, 11:02 AM
 */

#include "param_tools.h"

#include <fftw3.h> // For the inverse discrete cosine transform

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void Param_tools::save(Job_params job, std::ofstream& outFile) {

	int jobID = job.getInt("jobID");
	int verbose_save = job.getInt("verbose_save");
	int diagonalize = job.getInt("diagonalize");
	int d_weights = job.getInt("d_weights");
	int d_vecs = job.getInt("d_vecs");
	int d_cond = job.getInt("d_cond");
	int observable_type = job.getInt("observable_type");

	int poly_order = job.getInt("poly_order");
	double energy_rescale = job.getDouble("energy_rescale");
	double energy_shift = job.getDouble("energy_shift");

	std::vector<int> target_list = job.getIntVec("target_list");
	int num_targets = (int)target_list.size();

	if (verbose_save != 0){
		saveHeader(job, outFile);
	}

	if (verbose_save != 0){

		if (diagonalize == 0){
			if (observable_type == 0){

				if (diagonalize == 0){

				std::vector< std::vector<double> > cheb_coeffs = job.getDoubleMat("cheb_coeffs");

					if (jobID == 1){
						//print E vals:
						double g[poly_order];
						double E[poly_order];
						outFile << "E: ";

						for (int i = 0; i < poly_order; ++i){
							// Jackson coefficients
							g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(poly_order);
							E[i] = energy_shift + energy_rescale*cos((i*1.0 + 0.5)*M_PI/poly_order);
							outFile << E[i] << " ";
						}
						outFile << "\n";
					}

					outFile << "T: \n";

					for(int t = 0; t < num_targets; ++t){
						outFile << target_list[t] << ": ";
						for(int j = 0; j < poly_order-1; ++j){
							outFile << cheb_coeffs[t][j] << ", ";
						}
						outFile << cheb_coeffs[t][poly_order-1] << "\n";
					}

					outFile << "\n";
				}

			} else if (observable_type == 1){

				/*

				outFile << "T: \n";

				for(int t = 0; t < num_targets; ++t){
					outFile << target_list[t] << ": \n";
					for (int r = 0; r < poly_order; ++r){
						for(int j = 0; j < poly_order-1; ++j){
							outFile << result_array[job][j + r*poly_order + t*poly_order*poly_order] << ", ";
						}
						outFile << result_array[job][poly_order-1 + r*poly_order + t*poly_order*poly_order] << "\n";
					}
				}

				outFile << "\n";
				*/

			}
		} else if (diagonalize == 1){

			std::vector<double> eigenvalues = job.getDoubleVec("eigenvalues");

			int local_max_index = eigenvalues.size();
			outFile << "EIGS: ";
			for(int j = 0; j < local_max_index - 1; ++j){
				outFile << eigenvalues[j] << ", ";
			}
			outFile << eigenvalues[local_max_index - 1] << "\n";

			outFile << "\n";

			// Control for output printing
			// Depends on if eigenvectors (d_vecs) and conductivity (d_cond) are turned on or not

			if (d_weights == 1){

				std::vector< std::vector<double> > eigenweights = job.getDoubleMat("eigenweights");

				outFile << "WEIGHTS: \n";

				for(int t = 0; t < num_targets; ++t){
					outFile << target_list[t] << ": ";
					for(int j = 0; j < local_max_index - 1; ++j){
						outFile << eigenweights[t][j] << ", ";
					}
					outFile << eigenweights[t][local_max_index - 1] << "\n";
				}

				outFile << "\n";

			}

			if (d_vecs == 1){

				std::vector< std::vector<double> > eigenvectors = job.getDoubleMat("eigenvectors");

				outFile << "VECS: \n";
				for(int j = 0; j < local_max_index; ++j){
					for (int m = 0; m < local_max_index - 1; ++m){
						outFile << eigenvectors[j][m] << ", ";
					}
					outFile << eigenvectors[j][local_max_index - 1] << "\n";
				}

				outFile << "\n";

			}

			if (d_cond > 0){

				std::vector< std::vector<double> > M_xx = job.getDoubleMat("M_xx");

				outFile << "M_XX: \n";
				for(int j = 0; j < poly_order; ++j){
					for (int m = 0; m < poly_order - 1; ++m){
						outFile << M_xx[j][m] << ", ";
					}
					outFile << M_xx[j][poly_order - 1] << "\n";
				}

				outFile << "\n";

				if (d_cond > 1){

					std::vector< std::vector<double> > M_yy = job.getDoubleMat("M_yy");

					outFile << "M_YY: \n";
					for(int j = 0; j < poly_order; ++j){
						for (int m = 0; m < poly_order - 1; ++m){
							outFile << M_yy[j][m] << ", ";
						}
						outFile << M_yy[j][poly_order - 1] << "\n";
					}

					outFile << "\n";

					std::vector< std::vector<double> > M_xy = job.getDoubleMat("M_xy");

					outFile << "M_XY: \n";
					for(int j = 0; j < poly_order; ++j){
						for (int m = 0; m < poly_order - 1; ++m){
							outFile << M_xy[j][m] << ", ";
						}
						outFile << M_xy[j][poly_order - 1] << "\n";
					}

					outFile << "\n";
				}
				// Use to write M_xx to a binary file, M_XX_J<jobID>.dat

				/*

				int jobID = jobArray[job].getInt("jobID");

				std::string cwd = get_current_dir_name();

				system("mkdir temp");


				std::string M_xx_filename;
				M_xx_filename.append(cwd);
				M_xx_filename.append( "/temp/");

				if (jobID > 99){
					M_xx_filename.append("M_XX_J");
				} else if (jobID > 9){
					M_xx_filename.append("M_XX_J0");
				} else {
					M_xx_filename.append("M_XX_J00");
				}

				std::ostringstream temp_ss;

				temp_ss << jobID;
				M_xx_filename.append(temp_ss.str());
				M_xx_filename.append(".bin");

				// Debugging output
				//printf(M_xx_filename.c_str());
				//printf("\n");

				writeBufferToFile(&result_array[job][1 + b*local_max_index + a*local_max_index*local_max_index], poly_order*poly_order, M_xx_filename);
				*/

			}
		}
	} else { // non-verbose save
		if (diagonalize == 0){
			if (observable_type == 0){

				std::vector< std::vector<double> > cheb_coeffs = job.getDoubleMat("cheb_coeffs");

				if (jobID == 1){
					//print E vals:
					double g[poly_order];
					double E[poly_order];

					for (int i = 0; i < poly_order; ++i){
						// Jackson coefficients
						g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(poly_order);
						E[i] = energy_shift + energy_rescale*cos((i*1.0 + 0.5)*M_PI/poly_order);
						outFile << E[i] << " ";
					}
					outFile << "\n";
				}

				for(int t = 0; t < num_targets; ++t){
					for(int j = 0; j < poly_order-1; ++j){
						outFile << cheb_coeffs[t][j] << " ";
					}
					outFile << cheb_coeffs[t][poly_order-1] << "\n";
				}
			}
		} else if (diagonalize == 1){

			std::vector<double> eigenvalues = job.getDoubleVec("eigenvalues");

			int local_max_index = eigenvalues.size();
			for(int j = 0; j < local_max_index - 1; ++j){
				outFile << eigenvalues[j] << ", ";
			}
			outFile << eigenvalues[local_max_index - 1] << "\n";

			// Control for output printing
			// Depends on if eigenvectors (d_vecs) and conductivity (d_cond) are turned on or not
			/*
			if (d_weights == 1){

				outFile << "WEIGHTS: \n";

				for(int t = 0; t < num_targets; ++t){
					outFile << target_list[t] << ": ";
					for(int j = 0; j < local_max_index - 1; ++j){
						outFile << eigenweights[t][j] << ", ";
					}
					outFile << eigenweights[t][local_max_index - 1] << "\n";
				}

				outFile << "\n";

			}

			if (d_vecs == 1){


				outFile << "VECS: \n";
				for(int j = 0; j < local_max_index; ++j){
					for (int m = 0; m < local_max_index - 1; ++m){
						outFile << eigenvectors[j][m] << ", ";
					}
					outFile << eigenvectors[j][local_max_index - 1] << "\n";
				}

				outFile << "\n";

			}

			if (d_cond > 0){

				outFile << "M_XX: \n";
				for(int j = 0; j < poly_order; ++j){
					for (int m = 0; m < poly_order - 1; ++m){
						outFile << M_xx[j][m] << ", ";
					}
					outFile << M_xx[j][poly_order - 1] << "\n";
				}

				outFile << "\n";

					if (d_cond > 1){
					outFile << "M_YY: \n";
					for(int j = 0; j < poly_order; ++j){
						for (int m = 0; m < poly_order - 1; ++m){
							outFile << M_yy[j][m] << ", ";
						}
						outFile << M_yy[j][poly_order - 1] << "\n";
					}

					outFile << "\n";

					outFile << "M_XY: \n";
					for(int j = 0; j < poly_order; ++j){
						for (int m = 0; m < poly_order - 1; ++m){
							outFile << M_xy[j][m] << ", ";
						}
						outFile << M_xy[j][poly_order - 1] << "\n";
					}

					outFile << "\n";
				}
			}
			*/

		}
	}

}

void Param_tools::saveHeader(Job_params job, std::ofstream& outFile){

	int jobID = job.getInt("jobID");
	int solver_type = job.getInt("solver_type");

	std::vector<int> target_list = job.getIntVec("target_list");
	int num_targets = (int)target_list.size();

	int elecOn = job.getInt("elecOn");
	int magOn = job.getInt("magOn");
	double E = job.getDouble("E");
	double B = job.getDouble("B");

	outFile << "JOBID = " << jobID << " \n";

	if (solver_type == 1|| solver_type == 2) {

		int num_sheets = job.getInt("num_sheets");
		std::vector< std::vector<double> > shifts = job.getDoubleMat("shifts");

		outFile << "SHEET: SHIFT_X, SHIFT_Y, SHIFT_Z \n";
		for (int s = 0; s < num_sheets; ++s){
			outFile << s+1 << "    : ";
			outFile << shifts[s][0] << ", ";
			outFile << shifts[s][1] << ", ";
			outFile << shifts[s][2] << " \n";
		}

		outFile << "NUM_TAR = " << num_targets << "\n";
		if (num_targets != 0){
			outFile << "TAR_LIST: ";
			for (int t = 0; t < num_targets-1; ++t){
				outFile << target_list[t] << ", ";
			}
			outFile << target_list[num_targets - 1] << " \n";
		} else {
			outFile << "NO_TAR \n";
		}

	}

	if (solver_type == 3) {

	int mlmc_clusterID = job.getInt("mlmc_clusterID");

	outFile << "CLUSTERID = " <<  mlmc_clusterID << " \n";

		// /*
		outFile << "NUM_TAR = " << num_targets << "\n";
		if (num_targets != 0){
			outFile << "TAR_LIST: ";
			for (int t = 0; t < num_targets-1; ++t){
				outFile << target_list[t] << ", ";
			}
			outFile << target_list[num_targets - 1] << " \n";
		} else {
			outFile << "NO_TAR \n";
		}
		//*

		std::vector<int> vacancy_list = job.getIntVec("vacancy_list");
		int num_vacancies = (int)vacancy_list.size();

		outFile << "NUM_VAC = " << num_vacancies << "\n";
		if (vacancy_list[0] != -1){
			outFile << "VAC_LIST: ";
			for (int v = 0; v < num_vacancies-1; ++v){
				outFile << vacancy_list[v] << ", ";
			}
			outFile << vacancy_list[num_vacancies - 1] << " \n";
		} else {
			outFile << "NO_VAC \n";
		}

	}

	if (solver_type == 4) {

		outFile << "NUM_TAR = " << num_targets << "\n";
		if (num_targets != 0){
			outFile << "TAR_LIST: ";
			for (int t = 0; t < num_targets-1; ++t){
				outFile << target_list[t] << ", ";
			}
			outFile << target_list[num_targets - 1] << " \n";
		} else {
			outFile << "NO_TAR \n";
		}

		std::vector<int> vacancy_list = job.getIntVec("vacancy_list");
		int num_vacancies = (int)vacancy_list.size();

		outFile << "NUM_VAC = " << num_vacancies << "\n";
		if (vacancy_list[0] != -1){
			outFile << "VAC_LIST: ";
			for (int v = 0; v < num_vacancies-1; ++v){
				outFile << vacancy_list[v] << ", ";
			}
			outFile << vacancy_list[num_vacancies - 1] << " \n";
		} else {
			outFile << "NO_VAC \n";
		}

	}

	if (magOn == 1 || elecOn == 1) {
		outFile << "MAG_ON  = " << magOn  << ", B = " << B << " \n";
		outFile << "ELEC_ON = " << elecOn << ", E = " << E << " \n";
	}

}

void Param_tools::saveTiming(Job_params& job, double t, std::string tag){

	std::vector<double> temp_cpu_time;
	std::vector<string> temp_cpu_time_type;

	std::vector<string> double_vec_param_tags = job.getParamTags("doubleVec");

	// loop over doubleVec to see if we have a cpu_time variable already
	for (int i = 0; i < double_vec_param_tags.size(); ++i){
		if (double_vec_param_tags[i] == "cpu_time"){
			temp_cpu_time = job.getDoubleVec("cpu_time");
		}
	}

	std::vector<string> string_vec_param_tags = job.getParamTags("stringVec");

	// loop over stringVec to see if we have a cpu_time_type variable already
	for (int i = 0; i < string_vec_param_tags.size(); ++i){
		if (string_vec_param_tags[i] == "cpu_time_type"){
			temp_cpu_time_type = job.getStringVec("cpu_time_type");
		}
	}

	temp_cpu_time.push_back(t);
	temp_cpu_time_type.push_back(tag);

	job.setParam("cpu_time",temp_cpu_time);
	job.setParam("cpu_time_type",temp_cpu_time_type);

}

void Param_tools::densityTransform(Job_params& job) {

	int poly_order = job.getInt("poly_order");

	double energy_shift = job.getDouble("energy_shift");
	double energy_rescale = job.getDouble("energy_rescale");

	double g[poly_order];
	double E[poly_order];
	for (int i = 0; i < poly_order; ++i){
		// Jackson coefficients
		g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(poly_order);
		E[i] = energy_shift + energy_rescale*cos((i*1.0 + 0.5)*M_PI/poly_order);
	}

	//g[0] = 2*g[0];

	std::vector< std::vector<double> > cheb_coeffs = job.getDoubleMat("cheb_coeffs");

	for (int t = 0; t < cheb_coeffs.size(); ++t) {

		int p = cheb_coeffs[t].size();

		double* in = new double[p];

		for (int i = 0; i < p; ++i){
			in[i] = g[i]*cheb_coeffs[t][i];
		}

		double* out = new double[p];

		fftw_plan fftplan;
		fftplan = fftw_plan_r2r_1d(p,in,out,FFTW_REDFT01,FFTW_MEASURE);

		fftw_execute(fftplan);

		for (int i = 0; i < p; ++i){
			cheb_coeffs[t][i] = out[i]/(M_PI*sqrt( energy_rescale*energy_rescale - (E[i] - energy_shift)*(E[i] - energy_shift) ) );
		}

	}

	job.setParam("cheb_coeffs",cheb_coeffs);

}

void Param_tools::conductivityTransform(Job_params& job){

	int poly_order = job.getInt("poly_order");

	double g[poly_order];
	for (int i = 0; i < poly_order; ++i){
		// Jackson coefficients
		g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(sqrt(2.0)*poly_order);
	}

	int d_cond = job.getInt("d_cond");

	std::vector< std::vector<double> > M_xx;
	std::vector< std::vector<double> > M_yy;
	std::vector< std::vector<double> > M_xy;

	int dim_max = 0;
	if (d_cond > 0){
		dim_max = 1;
		M_xx = job.getDoubleMat("M_xx");
		if (d_cond > 1){
			dim_max = 3;
			M_yy = job.getDoubleMat("M_yy");
			M_xy = job.getDoubleMat("M_xy");
		}
	}

	for (int dim = 0; dim < dim_max; ++dim){

		double* in = new double[poly_order*poly_order];

		for (int i = 0; i < poly_order; ++i){
			for (int j = 0; j < poly_order; ++j){
				if (dim == 0) {
					in[i*poly_order + j] = 2*g[i]*g[j]*M_xx[i][j];
				} else if (dim == 1) {
					in[i*poly_order + j] = 2*g[i]*g[j]*M_yy[i][j];
				} else if (dim == 2) {
					in[i*poly_order + j] = 2*g[i]*g[j]*M_xy[i][j];
				}
			}
		}

		double* out = new double[poly_order*poly_order];

		fftw_plan p;
		p = fftw_plan_r2r_2d(poly_order,poly_order,in,out,FFTW_REDFT01,FFTW_REDFT01,FFTW_MEASURE);

		fftw_execute(p);

		for (int i = 0; i < poly_order; ++i){
			for (int j = 0; j < poly_order; ++j){
				if (dim == 0) {
					M_xx[i][j] = out[i*poly_order + j];
				} else if (dim == 1) {
					M_yy[i][j] = out[i*poly_order + j];
				} else if (dim == 2) {
					M_xy[i][j] = out[i*poly_order + j];
				}
			}
			// We put the Energies in the last column of the measure
		}

		fftw_destroy_plan(p);

		delete in;
		delete out;
	}

	if (d_cond > 0){
		job.setParam("M_xx",M_xx);
		if (d_cond > 1){
			job.setParam("M_yy",M_yy);
			job.setParam("M_xy",M_xy);
		}
	}

}

void Param_tools::mlmc_load(Job_params& job, std::string file_name){

	int d_cond = job.getInt("d_cond");

	// Use to load M_ij from a binary file
	if (d_cond > 0) {

		std::string file_name_xx = file_name;

		file_name_xx.append("M_xx.bin");

		std::vector< std::vector<double> > M_xx;


		std::ifstream file_xx(file_name_xx.c_str(), std::ifstream::binary);
		file_xx.seekg(0, file_xx.end);
		int sq_length_xx = (int) (file_xx.tellg()*sizeof(char))/sizeof(double);
		int length_xx = sqrt(sq_length_xx);
		file_xx.seekg(0, file_xx.beg);
		M_xx.resize(length_xx);
		for(int i = 0; i < length_xx; ++i){
			M_xx[i].resize(length_xx);
			char buffer[length_xx*(sizeof(double)/sizeof(char))];
			file_xx.read(buffer, length_xx*(sizeof(double)/sizeof(char)));
			double* doub_buffer = (double*) buffer;
			std::vector<double> temp_vec(doub_buffer, doub_buffer + length_xx);
			M_xx[i] = temp_vec;
		}
		file_xx.close();

		job.setParam("M_xx",M_xx);

		if (d_cond > 1){

			std::string file_name_yy = file_name;
			std::string file_name_xy = file_name;

			file_name_yy.append("M_yy.bin");
			file_name_xy.append("M_xy.bin");

			std::vector< std::vector<double> > M_yy;
			std::vector< std::vector<double> > M_xy;


			std::ifstream file_yy(file_name_yy.c_str(), std::ifstream::binary);
			file_yy.seekg(0, file_yy.end);
			int sq_length_yy = (file_yy.tellg()*sizeof(char))/sizeof(double);
			int length_yy = sqrt(sq_length_yy);
			file_yy.seekg(0, file_yy.beg);
			M_yy.resize(length_yy);
			for(int i = 0; i < length_yy; ++i){
				M_yy[i].resize(length_yy);
				char buffer[length_yy*(sizeof(double)/sizeof(char))];
				file_yy.read(buffer, length_yy*(sizeof(double)/sizeof(char)));
				double* doub_buffer = (double*) buffer;
				std::vector<double> temp_vec(doub_buffer, doub_buffer + length_yy);
				M_yy[i] = temp_vec;
			}
			file_yy.close();

			std::ifstream file_xy(file_name_xy.c_str(), std::ifstream::binary);
			file_xy.seekg(0, file_xy.end);
			int sq_length_xy = (file_xy.tellg()*sizeof(char))/sizeof(double);
			int length_xy = sqrt(sq_length_xy);
			file_xy.seekg(0, file_xy.beg);
			M_xy.resize(length_xy);
			for(int i = 0; i < length_xy; ++i){
				M_xy[i].resize(length_xy);
				char buffer[length_xy*(sizeof(double)/sizeof(char))];
				file_xy.read(buffer, length_xy*(sizeof(double)/sizeof(char)));
				double* doub_buffer = (double*) buffer;
				std::vector<double> temp_vec(doub_buffer, doub_buffer + length_xy);
				M_xy[i] = temp_vec;
			}
			file_xy.close();

			job.setParam("M_yy",M_yy);
			job.setParam("M_xy",M_xy);

		}
	}
}

void Param_tools::mlmc_save(Job_params job, std::string file_name){

	int d_cond = job.getInt("d_cond");

	// Use to write M_ij to a binary file
	if (d_cond > 0) {
		std::string file_name_xx = file_name;

		file_name_xx.append("M_xx.bin");

		std::vector< std::vector<double> > M_xx = job.getDoubleMat("M_xx");

		std::ofstream file_xx(file_name_xx.c_str(), std::ofstream::binary);
		for(size_t i = 0; i < M_xx.size(); i++ ){
			if ( M_xx[i].size() > 0 ){
			   char* buffer = (char*)(&M_xx[i][0]);
			   file_xx.write(buffer, M_xx[i].size()*(sizeof(double)/sizeof(char)));
			}
		}
		file_xx.close();

		if (d_cond > 1){

			std::string file_name_yy = file_name;
			std::string file_name_xy = file_name;

			file_name_yy.append("M_yy.bin");
			file_name_xy.append("M_xy.bin");

			std::vector< std::vector<double> > M_yy = job.getDoubleMat("M_yy");
			std::vector< std::vector<double> > M_xy = job.getDoubleMat("M_xy");


			std::ofstream file_yy(file_name_yy.c_str(), std::ofstream::binary);
			for(size_t i = 0; i < M_yy.size(); i++ ){
				if ( M_yy[i].size() > 0 ){
				   char* buffer = (char*)(&M_yy[i][0]);
				   file_yy.write(buffer, M_yy[i].size()*(sizeof(double)/sizeof(char)));
				}
			}
			file_yy.close();

			std::ofstream file_xy(file_name_xy.c_str(), std::ofstream::binary);
			for(size_t i = 0; i < M_xy.size(); i++ ){
				if ( M_xy[i].size() > 0 ){
				   char* buffer = (char*)(&M_xy[i][0]);
				   file_xy.write(buffer, M_xy[i].size()*(sizeof(double)/sizeof(char)));
				}
			}
			file_xy.close();

		}
	}
}

void Param_tools::mlmc_average(Job_params& total, Job_params samp){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	int d_cond = total.getInt("d_cond");

	if (d_cond > 0){

		std::vector< std::vector<double> > M_xx = total.getDoubleMat("M_xx");
		std::vector< std::vector<double> > samp_M_xx = samp.getDoubleMat("M_xx");

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			M_xx = samp_M_xx;
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + samp_M_xx[i][j]*s_s;
				}
			}
		}

		total.setParam("M_xx",M_xx);

		if (d_cond > 1){

			std::vector< std::vector<double> > M_yy = total.getDoubleMat("M_yy");
			std::vector< std::vector<double> > samp_M_yy = samp.getDoubleMat("M_yy");

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				M_yy = samp_M_yy;
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t +samp_M_yy[i][j]*s_s;
					}
				}
			}

			std::vector< std::vector<double> > M_xy = total.getDoubleMat("M_xy");
			std::vector< std::vector<double> > samp_M_xy = samp.getDoubleMat("M_xy");

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				M_xy = samp_M_xy;
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + samp_M_xy[i][j]*s_s;
					}
				}
			}

			total.setParam("M_yy",M_yy);
			total.setParam("M_xy",M_xy);

		}
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::mlmc_variance(Job_params& total, Job_params samp){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	int d_cond = total.getInt("d_cond");

	if (d_cond > 0){

		std::vector< std::vector<double> > M_xx = total.getDoubleMat("M_xx");
		std::vector< std::vector<double> > samp_M_xx = samp.getDoubleMat("M_xx");

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			M_xx.resize((int)samp_M_xx.size());
			for (int i = 0; i < (int) M_xx.size(); ++i){
				M_xx[i].resize((int)samp_M_xx[i].size());
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = pow(samp_M_xx[i][j],2);
				}
			}
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + pow(samp_M_xx[i][j],2)*s_s;
				}
			}
		}

		total.setParam("M_xx",M_xx);

		if (d_cond > 1){

			std::vector< std::vector<double> > M_yy = total.getDoubleMat("M_yy");
			std::vector< std::vector<double> > samp_M_yy = samp.getDoubleMat("M_yy");

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				M_yy.resize((int)samp_M_yy.size());
				for (int i = 0; i < (int) M_yy.size(); ++i){
					M_yy[i].resize((int)samp_M_yy[i].size());
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = pow(samp_M_yy[i][j],2);
					}
				}
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + pow(samp_M_yy[i][j],2)*s_s;
					}
				}
			}

			std::vector< std::vector<double> > M_xy = total.getDoubleMat("M_xy");
			std::vector< std::vector<double> > samp_M_xy = samp.getDoubleMat("M_xy");

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				M_xy.resize((int)samp_M_xy.size());
				for (int i = 0; i < (int) M_xy.size(); ++i){
					M_xy[i].resize((int)samp_M_xy[i].size());
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = pow(samp_M_xy[i][j],2);
					}
				}
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + pow(samp_M_xy[i][j],2)*s_s;
					}
				}
			}

			total.setParam("M_yy",M_yy);
			total.setParam("M_xy",M_xy);

		}
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::mlmc_cluster_average(Job_params& total, Job_params samp_orig, Job_params samp_cluster){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	int d_cond = total.getInt("d_cond");

	if (d_cond > 0){

		std::vector< std::vector<double> > M_xx = total.getDoubleMat("M_xx");
		std::vector< std::vector<double> > samp_orig_M_xx = samp_orig.getDoubleMat("M_xx");
		std::vector< std::vector<double> > samp_cluster_M_xx =samp_cluster.getDoubleMat("M_xx");

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			int size1 = (int) samp_orig_M_xx.size();
			int size2 = (int) samp_orig_M_xx[0].size();
			M_xx.resize(size1);
			for (int i = 0; i < size1; ++i){
				M_xx[i].resize(size2);
				for (int j = 0; j < size2; ++j){
					M_xx[i][j] = samp_orig_M_xx[i][j] - samp_cluster_M_xx[i][j];
				}
			}
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + (samp_orig_M_xx[i][j] - samp_cluster_M_xx[i][j])*s_s;
				}
			}
		}

		total.setParam("M_xx",M_xx);

		if (d_cond > 1) {

			std::vector< std::vector<double> > M_yy = total.getDoubleMat("M_yy");
			std::vector< std::vector<double> > samp_orig_M_yy = samp_orig.getDoubleMat("M_yy");
			std::vector< std::vector<double> > samp_cluster_M_yy =samp_cluster.getDoubleMat("M_yy");

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig_M_yy.size();
				int size2 = (int) samp_orig_M_yy[0].size();
				M_yy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_yy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_yy[i][j] = samp_orig_M_yy[i][j] - samp_cluster_M_yy[i][j];
					}
				}
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + (samp_orig_M_yy[i][j] - samp_cluster_M_yy[i][j])*s_s;
					}
				}
			}

			std::vector< std::vector<double> > M_xy = total.getDoubleMat("M_xy");
			std::vector< std::vector<double> > samp_orig_M_xy = samp_orig.getDoubleMat("M_xy");
			std::vector< std::vector<double> > samp_cluster_M_xy =samp_cluster.getDoubleMat("M_xy");


			if (M_xy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig_M_xy.size();
				int size2 = (int) samp_orig_M_xy[0].size();
				M_xy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_xy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_xy[i][j] = samp_orig_M_xy[i][j] - samp_cluster_M_xy[i][j];
					}
				}
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + (samp_orig_M_xy[i][j] - samp_cluster_M_xy[i][j])*s_s;
					}
				}
			}

			total.setParam("M_yy",M_yy);
			total.setParam("M_xy",M_xy);

		}
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::mlmc_cluster_variance(Job_params& total, Job_params samp_orig, Job_params samp_cluster){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	int d_cond = total.getInt("d_cond");

	if (d_cond > 0){

		std::vector< std::vector<double> > M_xx = total.getDoubleMat("M_xx");
		std::vector< std::vector<double> > samp_orig_M_xx = samp_orig.getDoubleMat("M_xx");
		std::vector< std::vector<double> > samp_cluster_M_xx =samp_cluster.getDoubleMat("M_xx");

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			int size1 = (int) samp_orig_M_xx.size();
			int size2 = (int) samp_orig_M_xx[0].size();
			M_xx.resize(size1);
			for (int i = 0; i < size1; ++i){
				M_xx[i].resize(size2);
				for (int j = 0; j < size2; ++j){
					M_xx[i][j] = pow(samp_orig_M_xx[i][j] - samp_cluster_M_xx[i][j],2);
				}
			}
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + pow(samp_orig_M_xx[i][j] - samp_cluster_M_xx[i][j],2)*s_s;
				}
			}
		}

		total.setParam("M_xx",M_xx);

		if (d_cond > 1){

			std::vector< std::vector<double> > M_yy = total.getDoubleMat("M_yy");
			std::vector< std::vector<double> > samp_orig_M_yy = samp_orig.getDoubleMat("M_yy");
			std::vector< std::vector<double> > samp_cluster_M_yy =samp_cluster.getDoubleMat("M_yy");

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig_M_yy.size();
				int size2 = (int) samp_orig_M_yy[0].size();
				M_yy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_yy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_yy[i][j] = pow(samp_orig_M_yy[i][j] - samp_cluster_M_yy[i][j],2);
					}
				}
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + pow(samp_orig_M_yy[i][j] - samp_cluster_M_yy[i][j],2)*s_s;
					}
				}
			}

			std::vector< std::vector<double> > M_xy = total.getDoubleMat("M_xy");
			std::vector< std::vector<double> > samp_orig_M_xy = samp_orig.getDoubleMat("M_xy");
			std::vector< std::vector<double> > samp_cluster_M_xy =samp_cluster.getDoubleMat("M_xy");

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig_M_xy.size();
				int size2 = (int) samp_orig_M_xy[0].size();
				M_xy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_xy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_xy[i][j] = pow(samp_orig_M_xy[i][j] - samp_cluster_M_xy[i][j],2);
					}
				}
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + pow(samp_orig_M_xy[i][j] - samp_cluster_M_xy[i][j],2)*s_s;
					}
				}
			}

			total.setParam("M_yy",M_yy);
			total.setParam("M_xy",M_xy);

		}
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}