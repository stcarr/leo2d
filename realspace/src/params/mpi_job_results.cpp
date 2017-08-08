/*
 * File:   mpi_job_results.cpp
 * Author: Stephen
 *
 * Created on May 24, 2017, 1:53 PM
 */

#include "mpi_job_results.h"

#include <fftw3.h> // For the inverse discrete cosine transform

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


Mpi_job_results::Mpi_job_results() {

	jobID = -1;
	max_jobs = -1;

	shifts = new double[3];
	shifts[0] = 0;
	shifts[1] = 0;
	shifts[2] = 0;

	num_target_sheets = 1;
	target_sheets = new int[1];
	target_sheets[0] = 0;

	num_targets = 1;
	target_list = new int[1];
	target_list[0] = 0;

	num_vacancies = 1;
	vacancy_list = new int[1];
	vacancy_list[0] = -1;

	energy_rescale = 20;
	energy_shift = 0;
	vacancy_chance = 0;

	solver_type = 0;
	observable_type = 0;
	solver_space = 0;
	verbose_save = 1;
	diagonalize = 0;
	d_weights = 1;
	d_vecs = 0;
	d_cond = 0;

	mlmc = 0;
	mlmc_clusterID = -1;
	mlmc_level = 1;
	mlmc_num_clusters = 0;
	mlmc_cluster_size = 4;
	mlmc_current_num_samples = 1;

	poly_order = 20;

	magOn = 0;
	elecOn = 0;
	B = 0;
	E = 0;

}

Mpi_job_results::Mpi_job_results(const Mpi_job_results& orig){

		jobID = orig.getInt("jobID");
		max_jobs = orig.getInt("max_jobs");

		energy_rescale = orig.getDouble("energy_rescale");
		energy_shift = orig.getDouble("energy_shift");
		vacancy_chance = orig.getDouble("vacancy_chance");

		solver_type = orig.getInt("solver_type");
		observable_type = orig.getInt("observable_type");
		solver_space = orig.getInt("solver_space");
		diagonalize = orig.getInt("diagonalize");
		d_weights = orig.getInt("d_weights");
		d_vecs = orig.getInt("d_vecs");
		d_cond = orig.getInt("d_cond");

		mlmc = orig.getInt("mlmc");
		mlmc_clusterID = orig.getInt("mlmc_clusterID");
		mlmc_level = orig.getInt("mlmc_level");
		mlmc_num_clusters = orig.getInt("mlmc_num_clusters");
		mlmc_cluster_size = orig.getInt("mlmc_cluster_size");

		num_target_sheets = orig.getInt("num_target_sheets");
		poly_order = orig.getInt("poly_order");

		magOn = orig.getInt("magOn");
		elecOn = orig.getInt("elecOn");
		B = orig.getDouble("B");
		E = orig.getDouble("E");

		target_sheets = orig.getIntVec("target_sheets");

		num_sheets = orig.getInt("num_sheets");
		shifts = orig.getDoubleMat("shifts");

		num_targets = orig.getInt("num_targets");
		target_list = orig.getIntVec("target_list");

		num_vacancies = orig.getInt("num_vacancies");
		vacancy_list = orig.getIntVec("vacancy_list");

		cheb_coeffs = orig.cheb_coeffs;
		eigenvalues = orig.eigenvalues;
		eigenweights = orig.eigenweights;
		eigenvectors = orig.eigenvectors;
		M_xx = orig.M_xx;
		M_yy = orig.M_yy;
		M_xy = orig.M_xy;

		cpu_time = orig.cpu_time;
		cpu_time_type = orig.cpu_time_type;

}

Mpi_job_results::~Mpi_job_results(){

}

void Mpi_job_results::loadLocParams(Job_params opts){

	energy_rescale = opts.getDouble("energy_rescale");
	energy_shift = opts.getDouble("energy_shift");
	vacancy_chance = opts.getDouble("vacancy_chance");

	solver_type = opts.getInt("solver_type");
	observable_type = opts.getInt("observable_type");
	solver_space = opts.getInt("solver_space");
	diagonalize = opts.getInt("diagonalize");
	verbose_save = opts.getInt("verbose_save");
	d_weights = opts.getInt("d_weights");
	d_vecs = opts.getInt("d_vecs");
	d_cond = opts.getInt("d_cond");

	mlmc = opts.getInt("mlmc");
	mlmc_level = opts.getInt("mlmc_level");
	mlmc_num_clusters = opts.getInt("mlmc_num_clusters");
	mlmc_cluster_size = opts.getInt("mlmc_cluster_size");

	num_target_sheets = opts.getInt("num_target_sheets");
	target_sheets = new int[num_target_sheets];
	std::vector<int> opts_target_sheets = opts.getIntVec("target_sheets");

	for (int i = 0; i < num_target_sheets; ++i){
		target_sheets[i] = opts_target_sheets[i];
	}

	poly_order = opts.getInt("poly_order");

	magOn = opts.getInt("magOn");
	elecOn = opts.getInt("elecOn");
	B = opts.getDouble("B");
	E = opts.getDouble("E");

}

void Mpi_job_results::loadJobParams(Job_params orig){

		jobID = orig.getInt("jobID");
		max_jobs = orig.getInt("max_jobs");

		energy_rescale = orig.getDouble("energy_rescale");
		energy_shift = orig.getDouble("energy_shift");
		vacancy_chance = orig.getDouble("vacancy_chance");

		solver_type = orig.getInt("solver_type");
		observable_type = orig.getInt("observable_type");
		solver_space = orig.getInt("solver_space");
		diagonalize = orig.getInt("diagonalize");
		verbose_save = orig.getInt("verbose_save");
		d_weights = orig.getInt("d_weights");
		d_vecs = orig.getInt("d_vecs");
		d_cond = orig.getInt("d_cond");

		mlmc = orig.getInt("mlmc");
		mlmc_clusterID = orig.getInt("mlmc_clusterID");
		mlmc_level = orig.getInt("mlmc_level");
		mlmc_num_clusters = orig.getInt("mlmc_num_clusters");
		mlmc_cluster_size = orig.getInt("mlmc_cluster_size");

		num_target_sheets = orig.getInt("num_target_sheets");
		poly_order = orig.getInt("poly_order");

		magOn = orig.getInt("magOn");
		elecOn = orig.getInt("elecOn");
		B = orig.getDouble("B");
		E = orig.getDouble("E");

		std::vector<int> target_sheets_in = orig.getIntVec("target_sheets");
		target_sheets = new int[num_target_sheets];
		for (int i = 0; i < num_target_sheets; ++i){
			target_sheets[i] = target_sheets_in[i];
		}

		num_sheets = orig.getInt("num_sheets");
		std::vector< std::vector<double> > shifts_in = orig.getDoubleMat("shifts");

		shifts = new double[(int)(shifts_in.size()*3)];
		for (int i = 0; i < shifts_in.size(); ++i){
			for (int j = 0; j < 3; ++j){
				shifts[i*3 + j] = shifts_in[i][j];
			}
		}

		num_targets = orig.getInt("num_targets");
		target_list = new int[num_targets];
		std::vector<int> target_list_in = orig.getIntVec("target_list");
		for (int i = 0; i < target_list_in.size();++i){
			target_list[i] = target_list_in[i];
		}

		num_vacancies = orig.getInt("num_vacancies");
		if (num_vacancies > 0){
		 	vacancy_list = new int[num_vacancies];
			std::vector<int> vacancy_list_in = orig.getIntVec("vacancy_list");
			for (int i = 0; i < vacancy_list_in.size();++i){
				vacancy_list[i] = vacancy_list_in[i];
			}
		}

}

void Mpi_job_results::saveTiming(double t, std::string tag){

	cpu_time.push_back(t);
	cpu_time_type.push_back(tag);

}

void Mpi_job_results::setParam(std::string tag, int val){

	if (tag == "jobID")
		jobID = val;
	if (tag == "max_jobs")
		max_jobs = val;
	if (tag == "solver_type")
		solver_type = val;
	if (tag == "observable_type")
		observable_type = val;
	if (tag == "solver_space")
		solver_space = val;
	if (tag == "diagonalize")
		diagonalize = val;
	if (tag == "verbose_save")
		verbose_save = val;
	if (tag == "d_weights")
		d_weights = val;
	if (tag == "d_vecs")
		d_vecs = val;
	if (tag == "d_cond")
		d_cond = val;
	if (tag == "poly_order")
		poly_order = val;
	if (tag == "magOn")
		magOn = val;
	if (tag == "elecOn")
		elecOn = val;
	if (tag == "mlmc")
		mlmc = val;
	if (tag == "mlmc_clusterID")
		mlmc_clusterID = val;
	if (tag == "mlmc_level")
		mlmc_level = val;
	if (tag == "mlmc_num_clusters")
		mlmc_num_clusters = val;
	if (tag == "mlmc_cluster_size")
		mlmc_cluster_size = val;

}

void Mpi_job_results::setParam(std::string tag, double val){

	if (tag == "energy_rescale")
		energy_rescale = val;
	if (tag == "energy_shift")
		energy_shift = val;
	if (tag == "B")
		B = val;
	if (tag == "E")
		E = val;
	if (tag == "vacancy_chance")
		vacancy_chance = val;

}

void Mpi_job_results::setParam(std::string tag, int *val, int dim){

	// target_sheets is soon to be removed! Replace with more general target_list.

	if (tag == "target_sheets") {
		num_target_sheets = dim;
		target_sheets = new int[num_target_sheets];
		for (int i = 0; i < num_target_sheets; ++i){
			target_sheets[i] = val[i];
		}
	} else if (tag == "target_list") {
		num_targets = dim;
		target_list = new int[num_targets];
		for (int i = 0; i < num_targets; ++i){
			target_list[i] = val[i];
		}
	} else if (tag == "vacancy_list") {
		num_vacancies = dim;
		vacancy_list = new int[num_vacancies];
		for (int i = 0; i < num_vacancies; ++i){
			vacancy_list[i] = val[i];
		}
	} else {
		printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());
		}
}

void Mpi_job_results::setParam(std::string tag, double *val, int dim){

}

void Mpi_job_results::setParam(std::string tag, int *val, int dim1, int dim2){

}

void Mpi_job_results::setParam(std::string tag, double *val, int dim1, int dim2){

	if (tag == "shifts") {

		if (dim2 != 3){
			printf("WARNING: (From: Mpi_job_results) cannot set shift parameter with dim2 != 3! \n");
		} else {

			num_sheets = dim1;

			shifts = new double[num_sheets*3];
			for (int s = 0; s < num_sheets; ++s){
				for (int i = 0; i < 3; ++i) {
					shifts[s*3 + i] = val[s*3 + i];
				}
			}
		}
	} else {
		printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());
	}
}

void Mpi_job_results::setParam(std::string tag, std::vector<double> val){

	if (tag == "eigenvalues")
		eigenvalues = val;
	if (tag == "cpu_time")
		cpu_time = val;

}

void Mpi_job_results::setParam(std::string tag, std::vector< std::vector<double> > val){

	if (tag == "cheb_coeffs")
		cheb_coeffs = val;
	if (tag == "eigenweights")
		eigenweights = val;
	if (tag == "eigenvectors")
		eigenvectors = val;
	if (tag == "M_xx")
		M_xx = val;
	if (tag == "M_yy")
		M_yy = val;
	if (tag == "M_xy")
		M_xy = val;

}

void Mpi_job_results::setParam(std::string tag, std::vector<std::string> val){

	if (tag == "cpu_time_type")
		cpu_time_type = val;

}

std::vector<double> Mpi_job_results::getResultVec(std::string tag) {

	if (tag == "eigenvalues")
		return eigenvalues;

}

std::vector< std::vector<double> > Mpi_job_results::getResultMat(std::string tag){

	if (tag == "cheb_coeffs")
		return cheb_coeffs;
	if (tag == "eigenweights")
		return eigenweights;
	if (tag == "eigenvectors")
		return eigenvectors;
	if (tag == "M_xx")
		return M_xx;
	if (tag == "M_yy")
		return M_yy;
	if (tag == "M_xy")
		return M_xy;

}

int Mpi_job_results::getInt(std::string tag) const{

	if (tag == "jobID")
		return jobID;
	if (tag == "max_jobs")
		return max_jobs;
	if (tag == "solver_type")
		return solver_type;
	if (tag == "observable_type")
		return observable_type;
	if (tag == "solver_space")
		return solver_space;
	if (tag == "diagonalize")
		return diagonalize;
	if (tag == "verbose_save")
		return verbose_save;
	if (tag == "d_weights")
		return d_weights;
	if (tag == "d_vecs")
		return d_vecs;
	if (tag == "d_cond")
		return d_cond;
	if (tag == "num_target_sheets")
		return num_target_sheets;
	if (tag == "poly_order")
		return poly_order;
	if (tag == "magOn")
		return magOn;
	if (tag == "elecOn")
		return elecOn;
	if (tag == "num_sheets")
		return num_sheets;
	if (tag == "num_targets")
		return num_targets;
	if (tag == "num_vacancies")
		return num_vacancies;
	if (tag == "mlmc")
		return mlmc;
	if (tag == "mlmc_clusterID")
		return mlmc_clusterID;
	if (tag == "mlmc_level")
		return mlmc_level;
	if (tag == "mlmc_num_clusters")
		return mlmc_num_clusters;
	if (tag == "mlmc_cluster_size")
		return mlmc_cluster_size;

	printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());
}

double Mpi_job_results::getDouble(std::string tag) const{

	if (tag == "energy_rescale")
		return energy_rescale;
	if (tag == "energy_shift")
		return energy_shift;
	if (tag == "B")
		return B;
	if (tag == "E")
		return E;
	if (tag == "vacancy_chance")
		return vacancy_chance;

	printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());

}

int* Mpi_job_results::getIntVec(std::string tag) const{

	if (tag == "target_sheets")
		return target_sheets;
	if (tag == "target_list")
		return target_list;
	if (tag == "vacancy_list")
		return vacancy_list;

	printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());

}

double* Mpi_job_results::getDoubleVec(std::string tag) const{

	printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());
}


int* Mpi_job_results::getIntMat(std::string tag) const{

	printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());
}

double* Mpi_job_results::getDoubleMat(std::string tag) const{

	if (tag == "shifts")
		return shifts;

	printf("WARNING: Mpi_job_results variable <%s> not found. \n", tag.c_str());

}

void Mpi_job_results::save(std::ofstream& outFile) {

	if (verbose_save != 0){
		printHeader(outFile);
	}

	if (verbose_save != 0){

		if (diagonalize == 0){
			if (observable_type == 0){

				if (diagonalize == 0){

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
	} else {
		if (diagonalize == 0){
			if (observable_type == 0){

				if (diagonalize == 0){

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
			}
		} else if (diagonalize == 1){

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

void Mpi_job_results::printParams(){

		printf("JOBID = %d settings: \n", jobID);
		printf("solver_type = %d \n", solver_type);
		printf("observable_type = %d \n", observable_type);
		printf("solver_space = %d \n",solver_space);
		printf("num_target_sheets = %d \n", num_target_sheets);
		printf("poly_order = %d \n", poly_order);
		printf("magOn = %d \n", magOn);
		printf("elecOn = %d \n", elecOn);
		printf("energy_rescale = %lf \n", energy_rescale);
		printf("energy_shift = %lf \n", energy_shift);
		printf("B = %lf \n", B);
		printf("E = %lf \n", E);
		printf("vacancy_chance = %lf \n",vacancy_chance);
}

void Mpi_job_results::printHeader(std::ofstream& outFile){

	outFile << "JOBID = " << jobID << " \n";

	if (solver_type == 1|| solver_type == 2) {

		outFile << "SHEET: SHIFT_X, SHIFT_Y, SHIFT_Z \n";
		for (int s = 0; s < num_sheets; ++s){
			outFile << s+1 << "    : ";
			outFile << shifts[s*3 + 0] << ", ";
			outFile << shifts[s*3 + 1] << ", ";
			outFile << shifts[s*3 + 2] << " \n";
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

void Mpi_job_results::mlmc_load(std::string file_name){

	// Use to load M_ij from a binary file
	if (d_cond > 0) {

		std::string file_name_xx = file_name;

		file_name_xx.append("M_xx.bin");


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

		if (d_cond > 1){

			std::string file_name_yy = file_name;
			std::string file_name_xy = file_name;

			file_name_yy.append("M_yy.bin");
			file_name_xy.append("M_xy.bin");

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

		}
	}
}

void Mpi_job_results::mlmc_save(std::string file_name){


	// Use to write M_ij to a binary file
	if (d_cond > 0) {
		std::string file_name_xx = file_name;

		file_name_xx.append("M_xx.bin");

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

void Mpi_job_results::mlmc_average(Mpi_job_results samp){

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	if (d_cond > 0){

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			M_xx = samp.M_xx;
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + samp.M_xx[i][j]*s_s;
				}
			}
		}

		if (d_cond > 1){

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				M_yy = samp.M_yy;
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + samp.M_yy[i][j]*s_s;
					}
				}
			}

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				M_xy = samp.M_xy;
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + samp.M_xy[i][j]*s_s;
					}
				}
			}

		}
	}

	mlmc_current_num_samples++;
}

void Mpi_job_results::mlmc_variance(Mpi_job_results samp){

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	if (d_cond > 0){

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			M_xx = samp.M_xx;
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + pow(samp.M_xx[i][j],2)*s_s;
				}
			}
		}

		if (d_cond > 1){

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				M_yy = samp.M_yy;
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + pow(samp.M_yy[i][j],2)*s_s;
					}
				}
			}

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				M_xy = samp.M_xy;
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + pow(samp.M_xy[i][j],2)*s_s;
					}
				}
			}

		}
	}

	mlmc_current_num_samples++;
}

void Mpi_job_results::mlmc_cluster_average(Mpi_job_results samp_orig,Mpi_job_results samp_cluster){

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	if (d_cond > 0){

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			int size1 = (int) samp_orig.M_xx.size();
			int size2 = (int) samp_orig.M_xx[0].size();
			M_xx.resize(size1);
			for (int i = 0; i < size1; ++i){
				M_xx[i].resize(size2);
				for (int j = 0; j < size2; ++j){
					M_xx[i][j] = samp_orig.M_xx[i][j] - samp_cluster.M_xx[i][j];
				}
			}
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + (samp_orig.M_xx[i][j] - samp_cluster.M_xx[i][j])*s_s;
				}
			}
		}

		if (d_cond > 1) {

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig.M_yy.size();
				int size2 = (int) samp_orig.M_yy[0].size();
				M_yy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_yy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_yy[i][j] = samp_orig.M_yy[i][j] - samp_cluster.M_yy[i][j];
					}
				}
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + (samp_orig.M_yy[i][j] - samp_cluster.M_yy[i][j])*s_s;
					}
				}
			}

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig.M_xy.size();
				int size2 = (int) samp_orig.M_xy[0].size();
				M_xy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_xy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_xy[i][j] = samp_orig.M_xy[i][j] - samp_cluster.M_xy[i][j];
					}
				}
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + (samp_orig.M_xy[i][j] - samp_cluster.M_xy[i][j])*s_s;
					}
				}
			}

		}
	}

	mlmc_current_num_samples++;
}

void Mpi_job_results::mlmc_cluster_variance(Mpi_job_results samp_orig,Mpi_job_results samp_cluster){

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	if (d_cond > 0){

		if (M_xx.empty() || mlmc_current_num_samples == 1){
			int size1 = (int) samp_orig.M_xx.size();
			int size2 = (int) samp_orig.M_xx[0].size();
			M_xx.resize(size1);
			for (int i = 0; i < size1; ++i){
				M_xx[i].resize(size2);
				for (int j = 0; j < size2; ++j){
					M_xx[i][j] = pow(samp_orig.M_xx[i][j] - samp_cluster.M_xx[i][j],2);
				}
			}
		} else {
			for (int i = 0; i < (int) M_xx.size(); ++i){
				for (int j = 0; j < (int) M_xx[i].size(); ++j){
					M_xx[i][j] = M_xx[i][j]*s_t + pow(samp_orig.M_xx[i][j] - samp_cluster.M_xx[i][j],2)*s_s;
				}
			}
		}

	if (d_cond > 1){

			if (M_yy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig.M_yy.size();
				int size2 = (int) samp_orig.M_yy[0].size();
				M_yy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_yy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_yy[i][j] = pow(samp_orig.M_yy[i][j] - samp_cluster.M_yy[i][j],2);
					}
				}
			} else {
				for (int i = 0; i < (int) M_yy.size(); ++i){
					for (int j = 0; j < (int) M_yy[i].size(); ++j){
						M_yy[i][j] = M_yy[i][j]*s_t + pow(samp_orig.M_yy[i][j] - samp_cluster.M_yy[i][j],2)*s_s;
					}
				}
			}

			if (M_xy.empty() || mlmc_current_num_samples == 1){
				int size1 = (int) samp_orig.M_xy.size();
				int size2 = (int) samp_orig.M_xy[0].size();
				M_xy.resize(size1);
				for (int i = 0; i < size1; ++i){
					M_xy[i].resize(size2);
					for (int j = 0; j < size2; ++j){
						M_xy[i][j] = pow(samp_orig.M_xy[i][j] - samp_cluster.M_xy[i][j],2);
					}
				}
			} else {
				for (int i = 0; i < (int) M_xy.size(); ++i){
					for (int j = 0; j < (int) M_xy[i].size(); ++j){
						M_xy[i][j] = M_xy[i][j]*s_t + pow(samp_orig.M_xy[i][j] - samp_cluster.M_xy[i][j],2)*s_s;
					}
				}
			}

		}
	}

	mlmc_current_num_samples++;
}

void Mpi_job_results::densityTransform() {

	double g[poly_order];
	double E[poly_order];
	for (int i = 0; i < poly_order; ++i){
		// Jackson coefficients
		g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(poly_order);
		E[i] = energy_shift + energy_rescale*cos((i*1.0 + 0.5)*M_PI/poly_order);
	}

	//g[0] = 2*g[0];

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
}

void Mpi_job_results::conductivityTransform(){

	double g[poly_order];
	for (int i = 0; i < poly_order; ++i){
		// Jackson coefficients
		g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(sqrt(2.0)*poly_order);
	}

	int dim_max = 0;
	if (d_cond > 0){
		dim_max = 1;
		if (d_cond > 1){
			dim_max = 3;
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

}


void Mpi_job_results::send(int target, int tag){

		sendInt(0,target,tag);

		sendInt(solver_type,target,tag);

		if (solver_type != -1) {

			sendInt(observable_type,target,tag);
			sendInt(solver_space,target,tag);
			sendInt(diagonalize,target,tag);
			sendInt(verbose_save,target,tag);
			sendInt(d_weights,target,tag);
			sendInt(d_vecs,target,tag);
			sendInt(d_cond,target,tag);

			sendInt(mlmc,target,tag);
			sendInt(mlmc_clusterID,target,tag);
			sendInt(mlmc_level,target,tag);
			sendInt(mlmc_num_clusters,target,tag);
			sendInt(mlmc_cluster_size,target,tag);

			sendInt(jobID,target,tag);
			sendInt(max_jobs,target,tag);

			sendDouble(energy_rescale,target,tag);
			sendDouble(energy_shift,target,tag);
			sendDouble(vacancy_chance,target,tag);

			sendInt(num_target_sheets,target,tag);
			sendIntVec(target_sheets,num_target_sheets,target,tag);
			sendInt(poly_order,target,tag);

			sendInt(magOn,target,tag);
			sendInt(elecOn,target,tag);
			sendDouble(B,target,tag);
			sendDouble(E,target,tag);

			sendInt(num_sheets,target,tag);
			sendDoubleMat(shifts,num_sheets,3,target,tag);

			sendIntVec(target_list, num_targets, target, tag);
			sendIntVec(vacancy_list, num_vacancies, target, tag);

			CPP_sendDoubleMat(cheb_coeffs, target, tag);
			CPP_sendDoubleVec(eigenvalues, target, tag);
			CPP_sendDoubleMat(eigenweights, target, tag);
			CPP_sendDoubleMat(eigenvectors, target, tag);
			CPP_sendDoubleMat(M_xx, target, tag);
			CPP_sendDoubleMat(M_yy, target, tag);
			CPP_sendDoubleMat(M_xy, target, tag);

			CPP_sendDoubleVec(cpu_time, target, tag);
			CPP_sendStrVec(cpu_time_type, target, tag);
		}

}

int Mpi_job_results::recv_spool(){

	int trash;

	MPI::Status status;
	MPI::COMM_WORLD.Recv(				// get size of incoming work
				&trash,
				1,
				MPI::INT,
				MPI::ANY_SOURCE,
				MPI::ANY_TAG,
				status);				// keeps tag and source information'

	int from = status.Get_source();
	recv(from);
	return from;

}


void Mpi_job_results::recv(int from){

		recvInt(from,"solver_type");

		if (solver_type != -1) {

			recvInt(from,"observable_type");
			recvInt(from,"solver_space");
			recvInt(from,"diagonalize");
			recvInt(from,"verbose_save");
			recvInt(from,"d_weights");
			recvInt(from,"d_vecs");
			recvInt(from,"d_cond");

			recvInt(from,"mlmc");
			recvInt(from,"mlmc_clusterID");
			recvInt(from,"mlmc_level,target");
			recvInt(from,"mlmc_num_clusters");
			recvInt(from,"mlmc_cluster_size");

			recvInt(from,"jobID");
			recvInt(from,"max_jobs");

			recvDouble(from,"energy_rescale");
			recvDouble(from,"energy_shift");
			recvDouble(from,"vacancy_chance");

			recvInt(from,"num_target_sheets");
			recvIntVec(from,"target_sheets");
			recvInt(from,"poly_order");

			recvInt(from,"magOn");
			recvInt(from,"elecOn");
			recvDouble(from,"B");
			recvDouble(from,"E");

			recvInt(from,"num_sheets");
			recvDoubleMat(from,"shifts");

			recvIntVec(from,"target_list");
			recvIntVec(from,"vacancy_list");

			CPP_recvDoubleMat(from,"cheb_coeffs");
			CPP_recvDoubleVec(from,"eigenvalues");
			CPP_recvDoubleMat(from,"eigenweights");
			CPP_recvDoubleMat(from,"eigenvectors");
			CPP_recvDoubleMat(from,"M_xx");
			CPP_recvDoubleMat(from,"M_yy");
			CPP_recvDoubleMat(from,"M_xy");

			CPP_recvDoubleVec(from,"cpu_time");
			CPP_recvStrVec(from,"cpu_time_type");

		}
}

void Mpi_job_results::sendInt(int val, int target, int tag){

	MPI::COMM_WORLD.Send(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label
}

void Mpi_job_results::sendDouble(double val, int target, int tag){

	MPI::COMM_WORLD.Send(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::DOUBLE,		// type of buffer
				target,				// rank to receive
				tag);				// tag to label
}

void Mpi_job_results::sendIntVec(int* val, int dim, int target, int tag){

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label


	MPI::COMM_WORLD.Send(
				val,		 		// input buffer
				dim,				// size of buffer
				MPI::INT,			// type of buffer
				target,				// r*ank to receive
				tag);				// tag to label

}

void Mpi_job_results::sendDoubleVec(double* val, int dim, int target, int tag){

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label


	MPI::COMM_WORLD.Send(
				val,		 		// input buffer
				dim,				// size of buffer
				MPI::DOUBLE,		// type of buffer
				target,				// rank to receive
				tag);				// tag to label

}


void Mpi_job_results::sendIntMat(int* val, int dim1, int dim2, int target, int tag){

	MPI::COMM_WORLD.Send(
				&dim1, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	MPI::COMM_WORLD.Send(
				&dim2, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	MPI::COMM_WORLD.Send(
				val,		 		// input buffer
				dim1*dim2,			// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// var to label

}

void Mpi_job_results::sendDoubleMat(double* val, int dim1, int dim2, int target, int tag){

	MPI::COMM_WORLD.Send(
				&dim1, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	MPI::COMM_WORLD.Send(
				&dim2, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	MPI::COMM_WORLD.Send(
				val,		 		// input buffer
				dim1*dim2,			// size of buffer
				MPI::DOUBLE,		// type of buffer
				target,				// rank to receive
				tag);				// var to label

}

void Mpi_job_results::CPP_sendStrVec(std::vector<std::string> val, int target, int tag){

	int dim;

	dim = val.size();

	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	if (dim != 0){
		for (int i = 0; i < dim; ++i) {

			int size = (int) val[i].length();

			MPI::COMM_WORLD.Send(
						&size, 				// input buffer
						1,					// size of buffer
						MPI::INT,			// type of buffer
						target,				// rank to receive
						tag);				// tag to label

			MPI::COMM_WORLD.Send(
						val[i].c_str(),		// input buffer
						size,				// size of buffer
						MPI::CHAR,			// type of buffer
						target,				// rank to receive
						tag);				// var to label

		}
	}

}

void Mpi_job_results::CPP_sendDoubleVec(std::vector<double> val, int target, int tag){

	int dim;

	dim = val.size();

	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	if (dim != 0){

		MPI::COMM_WORLD.Send(
					&val[0],		 	// input buffer
					dim,				// size of buffer
					MPI::DOUBLE,		// type of buffer
					target,				// rank to receive
					tag);				// var to label
	}

}

void Mpi_job_results::CPP_sendDoubleMat(std::vector< std::vector<double> > val, int target, int tag){

	int dim1;
	int dim2;

	dim1 = val.size();

	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim1, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label

	if (dim1 != 0){

		dim2 = val[0].size();

		MPI::COMM_WORLD.Send(
					&dim2, 				// input buffer
					1,					// size of buffer
					MPI::INT,			// type of buffer
					target,				// rank to receive
					tag);				// tag to label

		double temp_val[dim1*dim2];

		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2; ++j){
				temp_val[i*dim2 + j] = val[i][j];
			}
		}

		MPI::COMM_WORLD.Send(
					&temp_val,		 	// input buffer
					dim1*dim2,			// size of buffer
					MPI::DOUBLE,		// type of buffer
					target,				// rank to receive
					tag);				// var to label
	}

}

void Mpi_job_results::recvInt(int from, std::string var){

	int val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	setParam(var, val);

}

void Mpi_job_results::recvDouble(int from, std::string var){

	double val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::DOUBLE,		// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information


	setParam(var, val);

}

void Mpi_job_results::recvIntVec(int from, std::string var){

	int dim;
	int* val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
			&dim, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	val = new int[dim];

	MPI::COMM_WORLD.Recv(
			val,		 		// input buffer
			dim,				// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			status.Get_tag(),	// should be the same tag?
			status);			// keep MPI status information

	setParam(var,val,dim);

}

void Mpi_job_results::recvDoubleVec(int from, std::string var){

	int dim;
	double* val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
			&dim, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	val = new double[dim];

	MPI::COMM_WORLD.Recv(
			val,		 		// input buffer
			dim,				// size of buffer
			MPI::DOUBLE,		// type of buffer
			from,				// must come from "from"
			status.Get_tag(),	// should be the same tag?
			status);			// keep MPI status information

	setParam(var,val,dim);

}

void Mpi_job_results::recvIntMat(int from, std::string var){

	int dim1;
	int dim2;
	int* val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
			&dim1, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	MPI::COMM_WORLD.Recv(
			&dim2, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	val = new int[dim1*dim2];

	MPI::COMM_WORLD.Recv(
			val,		 		// input buffer
			dim1*dim2,			// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			status.Get_tag(),	// should be the same tag?
			status);			// keep MPI status information

	setParam(var,val,dim1,dim2);

}

void Mpi_job_results::recvDoubleMat(int from, std::string var){

	int dim1;
	int dim2;
	double* val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
			&dim1, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	MPI::COMM_WORLD.Recv(
			&dim2, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	val = new double[dim1*dim2];

	MPI::COMM_WORLD.Recv(
			val,		 		// input buffer
			dim1*dim2,			// size of buffer
			MPI::DOUBLE,		// type of buffer
			from,				// must come from "from"
			status.Get_tag(),	// should be the same tag?
			status);			// keep MPI status information

	setParam(var,val,dim1,dim2);

}

void Mpi_job_results::CPP_recvDoubleVec(int from, std::string var){

	int dim;
	std::vector<double> val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
			&dim, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	if (dim != 0){

		val.resize(dim);

		MPI::COMM_WORLD.Recv(
				&val[0],		 	// input buffer
				dim,				// size of buffer
				MPI::DOUBLE,		// type of buffer
				from,				// must come from "from"
				status.Get_tag(),	// should be the same tag?
				status);			// keep MPI status information

		setParam(var,val);
	}

}

void Mpi_job_results::CPP_recvStrVec(int from, std::string var){

	int dim;

	std::vector<std::string> string_vec;

	MPI::Status status;

	MPI::COMM_WORLD.Recv(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from "from"
				MPI::ANY_TAG,		// any tag
				status);			// keep MPI status information

	string_vec.resize(dim);

	if (dim != 0){
		for (int i = 0; i < dim; ++i) {

			int size;

			MPI::COMM_WORLD.Recv(
						&size, 				// input buffer
						1,					// size of buffer
						MPI::INT,			// type of buffer
						from,				// must come from "from"
						MPI::ANY_TAG,		// any tag
						status);			// keep MPI status information

			char input[size];

			MPI::COMM_WORLD.Recv(
						input,				// input buffer
						size,				// size of buffer
						MPI::CHAR,			// type of buffer
						from,				// must come from "from"
						MPI::ANY_TAG,	    // any tag
						status);			// keep MPI status information

			std::string string_temp(input,size);
			string_vec[i] = string_temp;

		}
		setParam(var,string_vec);
	}

}

void Mpi_job_results::CPP_recvDoubleMat(int from, std::string var){

	int dim1;
	int dim2;
	std::vector< std::vector<double> > val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
			&dim1, 				// input buffer
			1,					// size of buffer
			MPI::INT,			// type of buffer
			from,				// must come from "from"
			MPI::ANY_TAG,		// any tag
			status);			// keep MPI status information

	if (dim1 != 0){

		val.resize(dim1);

		MPI::COMM_WORLD.Recv(
				&dim2, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from "from"
				MPI::ANY_TAG,		// any tag
				status);			// keep MPI status information

		for (int i = 0; i < dim1; ++i){
			val[i].resize(dim2);
		}

		double temp_val[dim1*dim2];

		MPI::COMM_WORLD.Recv(
				temp_val,		 	// input buffer
				dim1*dim2,			// size of buffer
				MPI::DOUBLE,		// type of buffer
				from,				// must come from "from"
				status.Get_tag(),	// should be the same tag?
				status);			// keep MPI status information

		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2; ++j){
				val[i][j] = temp_val[i*dim2 + j];
				//printf("val[%d][%d] = %lf \n",i,j,val[i][j]);
			}
		}

		setParam(var,val);

	}

}
