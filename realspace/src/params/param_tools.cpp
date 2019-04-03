/*
 * File:   param_tools.cpp
 * Author: Stephen
 *
 * Created on August 21, 2017, 11:02 AM
 */

#include "param_tools.h"
#include "tools/numbers.h"

#include <fftw3.h> // For the inverse discrete cosine transform

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdexcept>

void Param_tools::save(Job_params job, std::ofstream& outFile) {

	outFile << fixed;

	int jobID = job.getInt("jobID");
	int verbose_save = job.getInt("verbose_save");
	int diagonalize = job.getInt("diagonalize");
	int d_kpm_dos = job.getInt("d_kpm_dos");
	int d_weights = job.getInt("d_weights");
	int d_vecs = job.getInt("d_vecs");
	int d_cond = job.getInt("d_cond");
	int observable_type = job.getInt("observable_type");
	int chiral_on = job.getInt("chiral_on");
	int ballistic_transport = job.getInt("ballistic_transport");

	int poly_order = job.getInt("poly_order");
	double energy_rescale = job.getDouble("energy_rescale");
	double energy_shift = job.getDouble("energy_shift");

	std::vector<int> target_list = job.getIntVec("target_list");
	int num_targets = (int)target_list.size();

	if (verbose_save != 0){
		saveHeader(job, outFile);
	}

	if (verbose_save == 1){

		if (diagonalize == 0){
			if (observable_type == 0){

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
		if (ballistic_transport == 1){

			std::vector<double> psi0_r = job.getDoubleVec("psi0_r");
		  std::vector<double> psi0_c = job.getDoubleVec("psi0_c");
			std::vector< std::vector<double> > psi_r = job.getDoubleMat("psi_r");
			std::vector< std::vector<double> > psi_c = job.getDoubleMat("psi_c");

			int N = psi0_r.size();
			int poly_max = psi_r.size();

			for (int k = 0; k < N; ++k){
				outFile << psi0_r[k] << " " << psi0_c[k] << " ";
				for (int p = 0; p < poly_max; ++p){
					outFile << psi_r[p][k] << " " << psi_c[p][k] << " ";
				}
				outFile << "\n";
			}

		} else if (diagonalize == 0){

			// KPM DOS
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

				// KPM Conductivity
			} else if (observable_type == 1){

				std::vector< std::vector<double> > kpm_M_xx = job.getDoubleMat("kpm_M_xx");

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
					// seperate E values from the conductivity matrices
					outFile << "\n";
					outFile << "\n";
				}

				for(int t = 0; t < num_targets; ++t){
					for (int i = 0; i < poly_order; ++i){
						for(int j = 0; j < poly_order; ++j){
							outFile << kpm_M_xx[t][j + i*poly_order] << " ";
						}
						outFile << "\n";
					}
					outFile << "\n";
				}

			}

		// now for diagonalization runs
		} else if (diagonalize == 1){

			if (d_kpm_dos == 1){

				if (jobID == 1){
					//print E vals:
					double E[poly_order];

					for (int i = 0; i < poly_order; ++i){
						E[i] = energy_shift + energy_rescale*cos((i*1.0 + 0.5)*M_PI/poly_order);
						outFile << E[i] << " ";
					}
					outFile << "\n";
				}

				std::vector<double> kpm_dos = job.getDoubleVec("kpm_dos");
				for (int i = 0; i < poly_order-1; ++i){
					outFile << kpm_dos[i] << ", ";
				}
				outFile << kpm_dos[poly_order - 1] << "\n";
			}

			if (chiral_on == 0 && d_kpm_dos == 0){
				if (d_cond  == 0) {
					std::vector<double> eigenvalues = job.getDoubleVec("eigenvalues");

					int local_max_index = eigenvalues.size();
					for(int j = 0; j < local_max_index - 1; ++j){
						outFile << eigenvalues[j] << ", ";
					}
					outFile << eigenvalues[local_max_index - 1] << "\n";
					if (d_weights == 1){

						std::vector< std::vector<double> > eigenweights = job.getDoubleMat("eigenweights");

						for(int t = 0; t < num_targets; ++t){
							for(int j = 0; j < local_max_index - 1; ++j){
								outFile << eigenweights[t][j] << ", ";
							}
							outFile << eigenweights[t][local_max_index - 1] << "\n";
						}

						outFile << "\n";

					}

					if (d_vecs == 1){

							int cpx_mat = job.getInt("complex_matrix");
							if (cpx_mat == 0){
								std::vector< std::vector<double> > eigenvectors = job.getDoubleMat("eigenvectors");
								for(int j = 0; j < local_max_index; ++j){
									for (int m = 0; m < local_max_index; ++m){
										outFile << j+1 << " " << m+1 << " " << eigenvectors[j][m] << " \n";
									}
								}
							} else if (cpx_mat == 1){
								std::vector< std::vector<double> > eigenvectors_r = job.getDoubleMat("eigenvectors_r");
								std::vector< std::vector<double> > eigenvectors_c = job.getDoubleMat("eigenvectors_c");

								for(int j = 0; j < local_max_index; ++j){
									for (int m = 0; m < local_max_index; ++m){
										outFile << j+1 << " " << m+1 << " " << eigenvectors_r[j][m] << " " << eigenvectors_c[j][m] << " \n";
									}
								}
							}

					}

				} else if (d_cond > 0){

					std::vector< std::vector<double> > M_xx = job.getDoubleMat("M_xx");

					for(int j = 0; j < poly_order; ++j){
						for (int m = 0; m < poly_order - 1; ++m){
							outFile << M_xx[j][m] << ", ";
						}
						outFile << M_xx[j][poly_order - 1] << "\n";
					}

					outFile << "\n";
				}

			} else if (chiral_on == 1){

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

				for (int matrix_count = 0; matrix_count < 6; ++matrix_count){

					std::string matrix_name;
					switch (matrix_count){
						case 0:
							matrix_name = "J_B_x";
						case 1:
							matrix_name = "J_B_y";
						case 2:
							matrix_name = "J_T_x";
						case 3:
							matrix_name = "J_T_y";
						case 4:
							matrix_name = "J_TB_x";
						case 5:
							matrix_name = "J_TB_y";
					}

					std::vector< std::vector< std::complex<double> > > tar_matrix = job.getCpxDoubleMat(matrix_name);

					int local_max_index = tar_matrix.size();

					for(int i = 0; i < local_max_index; ++i){
						for(int j = 0; j < local_max_index - 1; ++j){
							outFile << tar_matrix[i][j].real() << ", ";
						}
						outFile << tar_matrix[i][local_max_index - 1].real() << "\n";
					}
					outFile << "\n";

					for(int i = 0; i < local_max_index; ++i){
						for(int j = 0; j < local_max_index - 1; ++j){
							outFile << tar_matrix[i][j].imag() << ", ";
						}
						outFile << tar_matrix[i][local_max_index - 1].imag() << "\n";
					}
					outFile << "\n";

				}

			}

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

void Param_tools::saveKpts(Job_params job, std::ofstream& outFile) {

	outFile << fixed;
	outFile.precision(10);


	int jobID = job.getInt("jobID");
	int k_sampling = job.getInt("k_sampling");
	if (k_sampling == 1){
		std::vector<double> k_vec = job.getDoubleVec("k_vec");

		if (k_vec[0] >= 0){
			outFile << " ";
		}

		outFile << k_vec[0] << "		";
		if (k_vec[1] >= 0){
			outFile << " ";
		}

		outFile << k_vec[1] << " \n";
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

// Compute area of the brillioun zone for a giving real-space unitcell
double Param_tools::computeReciprocalArea(vector< vector<double> > uc){
	// b1 x b2 = (2 pi)^2/((a1 dot Ra2)*(a2 dot Ra1)) * (a1 x a2)

	// 90 degree CCW rotation matrix R acting on uc
	std::vector< std::vector<double> > R_uc = uc;
	R_uc[0][0] = -uc[0][1];
	R_uc[0][1] =  uc[0][0];
	R_uc[1][0] = -uc[1][1];
	R_uc[1][1] =  uc[1][0];

	double cross = uc[0][0]*uc[1][1] - uc[0][1]*uc[1][0];
	double denom_1 = uc[0][0]*R_uc[1][0] + uc[0][1]*R_uc[1][1];
	double denom_2 = uc[1][0]*R_uc[0][0] + uc[1][1]*R_uc[0][1];

	double area = (2*numbers::PI)*(2*numbers::PI)*cross/(denom_1 * denom_2);

	if (area < 0){
		area = -area;
	}

	//printf("computeReciprocalArea: uc = [%lf, %lf; %lf, %lf], area = %lf \n",uc[0][0],uc[0][1],uc[1][0],uc[1][1],area);

	return area;
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

void Param_tools::densityTransformTrace(Job_params& job, double scale_fac) {

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

	int n_s = cheb_coeffs.size(); // number of trace samples
	int n_E = cheb_coeffs[0].size(); // poly_order <-> Energy sampling

	std::vector<double> kpm_trace_moments;
	kpm_trace_moments.resize(n_E);

	for (int i_E = 0; i_E < n_E; ++i_E){

		kpm_trace_moments[i_E] = 0.0;

		// normalize by number of trace samples
		for (int i_s = 0; i_s < n_s; ++i_s){
			kpm_trace_moments[i_E] = kpm_trace_moments[i_E] + cheb_coeffs[i_s][i_E]/( (double)n_s  );
		}

		//printf("u[%d] = %le \n",i_E,kpm_trace_moments[i_E]);

	}

	int p = kpm_trace_moments.size();

	double* in = new double[p];

	for (int i = 0; i < p; ++i){
		in[i] = g[i]*kpm_trace_moments[i];
	}

	double* out = new double[p];

	fftw_plan fftplan;
	fftplan = fftw_plan_r2r_1d(p,in,out,FFTW_REDFT01,FFTW_MEASURE);

	fftw_execute(fftplan);

	std::vector<double> kpm_trace_dos;
	kpm_trace_dos.resize(p);
	for (int i = 0; i < p; ++i){
		kpm_trace_dos[i] = scale_fac*out[i]/(M_PI*sqrt( energy_rescale*energy_rescale - (E[i] - energy_shift)*(E[i] - energy_shift) ) );
	}

	job.setParam("kpm_trace_dos",kpm_trace_dos);

}

void Param_tools::conductivityTransform(Job_params& job){

	int diagonalize = job.getInt("diagonalize");
	int observable_type = job.getInt("observable_type");
	int chiral_on = job.getInt("chiral_on");
	int d_cond = job.getInt("d_cond");

	if (diagonalize == 0 && observable_type == 1){

		int num_targets = job.getInt("num_targets");
		int poly_order = job.getInt("poly_order");
		std::vector< std::vector<double> > kpm_M_xx_here = job.getDoubleMat("kpm_M_xx");
		std::vector< std::vector<double> > new_kpm_M_xx;
		new_kpm_M_xx.resize(num_targets);

		for (int t = 0; t < num_targets; ++t){
			new_kpm_M_xx[t].resize(poly_order*poly_order);

			Job_params temp_job(job);
			std::vector< std::vector<double> > temp_kpm_M_xx;
			temp_kpm_M_xx.resize(poly_order);

			for (int p1 = 0; p1 < poly_order; ++p1){
				temp_kpm_M_xx[p1].resize(poly_order);
				for (int p2 = 0; p2 < poly_order; ++p2){
					temp_kpm_M_xx[p1][p2] = kpm_M_xx_here[t][p2 + p1*poly_order];
				}
			}

			temp_job.setParam("temp_kpm_M_xx",temp_kpm_M_xx);
			Param_tools::matrixResponseTransform(temp_job,"temp_kpm_M_xx");
			temp_kpm_M_xx = temp_job.getDoubleMat("temp_kpm_M_xx");

			for (int p1 = 0; p1 < poly_order; ++p1){
				for (int p2 = 0; p2 < poly_order; ++p2){
					new_kpm_M_xx[t][p2 + p1*poly_order] = temp_kpm_M_xx[p1][p2];
				}
			}

			job.setParam("kpm_M_xx",new_kpm_M_xx);

		} // end t loop

	} else {
		if (d_cond > 0){

			Param_tools::matrixResponseTransform(job,"M_xx");

			if (d_cond > 1){
				Param_tools::matrixResponseTransform(job,"M_yy");
				Param_tools::matrixResponseTransform(job,"M_xy");
			}
		}

		if (chiral_on == 1){

			Param_tools::matrixResponseTransform(job,"J_B_x");
			Param_tools::matrixResponseTransform(job,"J_B_y");
			Param_tools::matrixResponseTransform(job,"J_T_x");
			Param_tools::matrixResponseTransform(job,"J_T_y");
			Param_tools::matrixResponseTransform(job,"J_TB_x");
			Param_tools::matrixResponseTransform(job,"J_TB_y");

		}
	}

}

void Param_tools::matrixResponseTransform(Job_params& job, std::string tag){

	int poly_order = job.getInt("poly_order");

	std::vector< std::vector<double> > matrixIn;
	std::vector< std::vector<double> > matrixOut;

	std::vector< std::vector< std::complex<double> > > matrixIn_cpx;
	std::vector< std::vector< std::complex<double> > > matrixOut_cpx;

	int type = 0;

	try{
		matrixIn = job.getDoubleMat(tag);
	} catch(const std::invalid_argument& ia){
		matrixIn_cpx = job.getCpxDoubleMat(tag);
		type = 1;
	}

	double g[poly_order];
	for (int i = 0; i < poly_order; ++i){
		// Jackson coefficients
		g[i] = ((poly_order-i)*cos((M_PI*i)/poly_order) + sin((M_PI*i)/poly_order)/tan(M_PI/poly_order))/(sqrt(2.0)*poly_order);
	}

	if (type == 0){

		//double* in = new double[poly_order*poly_order];
		//double* out = new double[poly_order*poly_order];

		double *in, *out;
		in = (double*) fftw_malloc( sizeof(double)*poly_order*poly_order);
		out = (double*) fftw_malloc( sizeof(double)*poly_order*poly_order);

		// fftw_plan p;
		// p = fftw_plan_r2r_2d(poly_order,poly_order,in,out,FFTW_REDFT01,FFTW_REDFT01,FFTW_ESTIMATE);

		fftw_plan p;
		p = fftw_plan_r2r_2d(poly_order,poly_order,in,out,FFTW_REDFT01,FFTW_REDFT01,FFTW_MEASURE);

		if (!p){
			throw std::runtime_error("FFTW plan is NULL (in Param_tools::matrixResponseTransform)");
		}

		for (int i = 0; i < poly_order; ++i){
			for (int j = 0; j < poly_order; ++j){
				in[i*poly_order + j] = 2*g[i]*g[j]*matrixIn[i][j];
				out[i*poly_order + j] = 99.123;
				//if (j < 20 && i > 950)
					//printf("in[%d][%d] = %le  |  out[%d][%d] = %le \n",i,j,in[i*poly_order + j],i,j,out[i*poly_order + j]);
				}
		}

		fftw_execute(p);

		matrixOut.resize(poly_order);
		for (int i = 0; i < poly_order; ++i){
			matrixOut[i].resize(poly_order);
			for (int j = 0; j < poly_order; ++j){
				matrixOut[i][j] = out[i*poly_order + j];
				//if (j < 20 && i > 950)
					//printf("in[%d][%d] = %le  |  out[%d][%d] = %le \n",i,j,in[i*poly_order + j],i,j,out[i*poly_order + j]);
			}
		}

		job.setParam(tag,matrixOut);

		fftw_destroy_plan(p);
		fftw_free(in); fftw_free(out);

		//delete in;
		//delete out;

	} else if (type == 1){

		// real component transformation
		double* in_real = new double[poly_order*poly_order];
		double* out_real = new double[poly_order*poly_order];

		fftw_plan p_real;
		p_real = fftw_plan_r2r_2d(poly_order,poly_order,in_real,out_real,FFTW_REDFT01,FFTW_REDFT01,FFTW_MEASURE);

		for (int i = 0; i < poly_order; ++i){
			for (int j = 0; j < poly_order; ++j){
				in_real[i*poly_order + j] = 2*g[i]*g[j]*matrixIn_cpx[i][j].real();
			}
		}


		fftw_execute(p_real);

		matrixOut_cpx.resize(poly_order);
		for (int i = 0; i < poly_order; ++i){
			matrixOut_cpx[i].resize(poly_order);
			for (int j = 0; j < poly_order; ++j){
				matrixOut_cpx[i][j].real(out_real[i*poly_order + j]);
			}
		}

		fftw_destroy_plan(p_real);

		delete in_real;
		delete out_real;

		// imaginary component transformation
		double* in_imag = new double[poly_order*poly_order];
		double* out_imag = new double[poly_order*poly_order];

		fftw_plan p_imag;
		p_imag = fftw_plan_r2r_2d(poly_order,poly_order,in_imag,out_imag,FFTW_REDFT01,FFTW_REDFT01,FFTW_MEASURE);

		for (int i = 0; i < poly_order; ++i){
			for (int j = 0; j < poly_order; ++j){
				in_imag[i*poly_order + j] = 2*g[i]*g[j]*matrixIn_cpx[i][j].imag();
			}
		}

		fftw_execute(p_imag);

		for (int i = 0; i < poly_order; ++i){
			for (int j = 0; j < poly_order; ++j){
				matrixOut_cpx[i][j].imag(out_imag[i*poly_order + j]);
			}
		}

		fftw_destroy_plan(p_imag);

		delete in_imag;
		delete out_imag;

		job.setParam(tag,matrixOut_cpx);

	}

}


void Param_tools::mlmc_load(Job_params& job, std::string file_name){

	int d_cond = job.getInt("d_cond");
	int d_kpm_dos = job.getInt("d_kpm_dos");
	int kpm_trace = job.getInt("kpm_trace");

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

	// Use to load dos from a binary file
	if (d_kpm_dos == 1) {

		std::string file_name_dos = file_name;

		file_name_dos.append("dos.bin");

		std::vector<double> dos;


		std::ifstream file_dos(file_name_dos.c_str(), std::ifstream::binary);
		file_dos.seekg(0, file_dos.end);
		int length_dos = (int) (file_dos.tellg()*sizeof(char))/sizeof(double);
		file_dos.seekg(0, file_dos.beg);

		char buffer[length_dos*(sizeof(double)/sizeof(char))];
		file_dos.read(buffer, length_dos*(sizeof(double)/sizeof(char)));
		double* doub_buffer = (double*) buffer;
		std::vector<double> temp_vec(doub_buffer, doub_buffer + length_dos);
		dos = temp_vec;

		file_dos.close();

		job.setParam("kpm_dos",dos);
	}

	// Use to load dos from a binary file
	if (kpm_trace == 1) {

		std::string file_name_dos = file_name;

		file_name_dos.append("dos.bin");

		std::vector<double> dos;


		std::ifstream file_dos(file_name_dos.c_str(), std::ifstream::binary);
		file_dos.seekg(0, file_dos.end);
		int length_dos = (int) (file_dos.tellg()*sizeof(char))/sizeof(double);
		file_dos.seekg(0, file_dos.beg);

		char buffer[length_dos*(sizeof(double)/sizeof(char))];
		file_dos.read(buffer, length_dos*(sizeof(double)/sizeof(char)));
		double* doub_buffer = (double*) buffer;
		std::vector<double> temp_vec(doub_buffer, doub_buffer + length_dos);
		dos = temp_vec;

		file_dos.close();

		job.setParam("kpm_trace_dos",dos);
	}

}

void Param_tools::mlmc_save(Job_params job, std::string file_name){

	int d_cond = job.getInt("d_cond");
	int d_kpm_dos = job.getInt("d_kpm_dos");
	int kpm_trace = job.getInt("kpm_trace");

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

	if (d_kpm_dos == 1){

		std::string file_name_dos = file_name;
		file_name_dos.append("dos.bin");

		std::vector<double> dos = job.getDoubleVec("kpm_dos");

		std::ofstream file_dos(file_name_dos.c_str(), std::ofstream::binary);

		if ( dos.size() > 0 ){
			for (int p = 0; p < dos.size(); ++p){
				//printf("dos[%d] = %lf \n",p,dos[p]);
			}
			 char* buffer = (char*)(&dos[0]);
			 file_dos.write(buffer, dos.size()*(sizeof(double)/sizeof(char)));
		}

		file_dos.close();

	}

	if (kpm_trace == 1){

		std::string file_name_dos = file_name;
		file_name_dos.append("dos.bin");

		std::vector<double> dos = job.getDoubleVec("kpm_trace_dos");

		std::ofstream file_dos(file_name_dos.c_str(), std::ofstream::binary);

		if ( dos.size() > 0 ){
			for (int p = 0; p < dos.size(); ++p){
				//printf("dos[%d] = %lf \n",p,dos[p]);
			}
			 char* buffer = (char*)(&dos[0]);
			 file_dos.write(buffer, dos.size()*(sizeof(double)/sizeof(char)));
		}

		file_dos.close();

	}

}

void Param_tools::mlmc_average(Job_params& total, Job_params samp){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	//printf("mlmc_averaging: (%d) s_t = %lf, s_s = %lf\n",mlmc_current_num_samples,s_t,s_s);

	int d_cond = total.getInt("d_cond");
	int d_kpm_dos = total.getInt("d_kpm_dos");
	int kpm_trace = total.getInt("kpm_trace");

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

	if (d_kpm_dos == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_dos");
		std::vector<double> samp_dos = samp.getDoubleVec("kpm_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			dos = samp_dos;
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + samp_dos[i]*s_s;
				//printf("averaged_dos[%d] = %lf \n", i,dos[i]);
			}
		}

		total.setParam("kpm_dos",dos);
	}

	if (kpm_trace == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_trace_dos");
		std::vector<double> samp_dos = samp.getDoubleVec("kpm_trace_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			dos = samp_dos;
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + samp_dos[i]*s_s;
				//printf("averaged_dos[%d] = %lf \n", i,dos[i]);
			}
		}

		total.setParam("kpm_trace_dos",dos);

	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::mlmc_variance(Job_params& total, Job_params samp){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	int d_cond = total.getInt("d_cond");
	int d_kpm_dos = total.getInt("d_kpm_dos");
	int kpm_trace = total.getInt("kpm_trace");

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

	if (d_kpm_dos == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_dos");
		std::vector<double> samp_dos = samp.getDoubleVec("kpm_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			dos.resize((int)samp_dos.size());
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = pow(samp_dos[i],2);
			}
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + pow(samp_dos[i],2)*s_s;
			}
		}

		total.setParam("kpm_dos",dos);
	}

	if (kpm_trace == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_trace_dos");
		std::vector<double> samp_dos = samp.getDoubleVec("kpm_trace_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			dos.resize((int)samp_dos.size());
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = pow(samp_dos[i],2);
			}
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + pow(samp_dos[i],2)*s_s;
			}
		}

		total.setParam("kpm_trace_dos",dos);
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::mlmc_cluster_average(Job_params& total, Job_params samp_orig, Job_params samp_cluster){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp


	int d_cond = total.getInt("d_cond");
	int d_kpm_dos = total.getInt("d_kpm_dos");
	int kpm_trace = total.getInt("kpm_trace");

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

	if (d_kpm_dos == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_dos");
		std::vector<double> samp_orig_dos = samp_orig.getDoubleVec("kpm_dos");
		std::vector<double> samp_cluster_dos = samp_cluster.getDoubleVec("kpm_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			int size = (int) samp_orig_dos.size();
			dos.resize(size);
			for (int i = 0; i < size; ++i){
				dos[i] = samp_orig_dos[i] - samp_cluster_dos[i];
			}
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + (samp_orig_dos[i] - samp_cluster_dos[i])*s_s;
			}
		}

		total.setParam("kpm_dos",dos);
	}

	if (kpm_trace == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_trace_dos");
		std::vector<double> samp_orig_dos = samp_orig.getDoubleVec("kpm_trace_dos");
		std::vector<double> samp_cluster_dos =samp_cluster.getDoubleVec("kpm_trace_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			int size = (int) samp_orig_dos.size();
			dos.resize(size);
			for (int i = 0; i < size; ++i){
				dos[i] = samp_orig_dos[i] - samp_cluster_dos[i];
			}
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + (samp_orig_dos[i] - samp_cluster_dos[i])*s_s;
			}
		}

		total.setParam("kpm_trace_dos",dos);
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::mlmc_cluster_variance(Job_params& total, Job_params samp_orig, Job_params samp_cluster){

	int mlmc_current_num_samples = total.getInt("mlmc_current_num_samples");

	double s_t = (1.0*mlmc_current_num_samples - 1)/(1.0*mlmc_current_num_samples); // scale for total
	double s_s = 1.0/(1.0*mlmc_current_num_samples); // scale for samp

	int d_cond = total.getInt("d_cond");
	int d_kpm_dos = total.getInt("d_kpm_dos");
	int kpm_trace = total.getInt("kpm_trace");

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

	if (d_kpm_dos == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_dos");
		std::vector<double> samp_orig_dos = samp_orig.getDoubleVec("kpm_dos");
		std::vector<double> samp_cluster_dos =samp_cluster.getDoubleVec("kpm_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			int size = (int) samp_orig_dos.size();
			dos.resize(size);
			for (int i = 0; i < size; ++i){
				dos[i] = pow(samp_orig_dos[i] - samp_cluster_dos[i],2);
			}
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + pow(samp_orig_dos[i] - samp_cluster_dos[i],2)*s_s;
			}
		}

		total.setParam("kpm_dos",dos);
	}

	if (kpm_trace == 1){

		std::vector<double> dos = total.getDoubleVec("kpm_trace_dos");
		std::vector<double> samp_orig_dos = samp_orig.getDoubleVec("kpm_trace_dos");
		std::vector<double> samp_cluster_dos =samp_cluster.getDoubleVec("kpm_trace_dos");

		if (dos.empty() || mlmc_current_num_samples == 1){
			int size = (int) samp_orig_dos.size();
			dos.resize(size);
			for (int i = 0; i < size; ++i){
				dos[i] = pow(samp_orig_dos[i] - samp_cluster_dos[i],2);
			}
		} else {
			for (int i = 0; i < (int) dos.size(); ++i){
				dos[i] = dos[i]*s_t + pow(samp_orig_dos[i] - samp_cluster_dos[i],2)*s_s;
			}
		}

		total.setParam("kpm_trace_dos",dos);
	}

	total.setParam("mlmc_current_num_samples",mlmc_current_num_samples+1);
}

void Param_tools::wannier_raw_eig_save(std::vector<Job_params>& results, std::ofstream& outFile_vals, std::ofstream& outFile_vecs){

	outFile_vals << fixed;
	outFile_vecs << fixed;
	outFile_vals.precision(12);
	outFile_vecs.precision(12);

	int num_k1 = results[0].getInt("num_k1");
	int num_k2 = results[0].getInt("num_k2");
	int num_k = num_k1*num_k2;
	std::vector<double> eigenvalues = results[0].getDoubleVec("eigenvalues");
	int num_bands = eigenvalues.size();

	for (int k = 0; k < num_k; ++k){
		eigenvalues = results[k].getDoubleVec("eigenvalues");
		for (int b = 0; b < num_bands; ++b){
			outFile_vals << (eigenvalues[b] >= 0 ? " ":"") << eigenvalues[b] << "   ";
		}
		outFile_vals << "\n";
	}

	std::vector< std::vector<double> > vecs_r;
	std::vector< std::vector<double> > vecs_c;
	for (int k = 0; k < num_k; ++k){
		vecs_r = results[k].getDoubleMat("eigenvectors_r");
		vecs_c = results[k].getDoubleMat("eigenvectors_c");

		int num_orbs = vecs_r[0].size();

		for (int b = 0; b < num_bands; ++b){
			for (int o = 0; o < num_orbs; ++o){
			outFile_vecs << vecs_r[b][o] << (vecs_c[b][o] >= 0 ? "+":"") << vecs_c[b][o] << "j   ";
			}
			outFile_vecs << "\n";
		}
		outFile_vecs << "\n";
	}

}

void Param_tools::wannier_win_save(std::vector<Job_params>& results, std::ofstream& outFile){

	outFile << fixed;
	outFile.precision(12);

	int num_k1 = results[0].getInt("num_k1");
	int num_k2 = results[0].getInt("num_k2");
	int num_k = num_k1*num_k2;
	std::vector<double> eigenvalues = results[0].getDoubleVec("eigenvalues");
	int num_bands = eigenvalues.size();
	int num_wan = 4;
	std::vector< std::vector<double> > unitcell = results[0].getDoubleMat("unit_cell");

	outFile << "# *.win file generated by LEO2D \n";
	outFile << "mp_grid " << num_k1 << " " << num_k2 << " 1 \n";
	outFile << "num_bands " << num_bands << " \n";
	outFile << "num_wann " << num_wan << " \n";
	outFile << "\n";

	outFile << "begin unit_cell_cart \n";
	outFile << "Angstroms \n";
		for (int d = 0; d < 2; ++d){
			outFile << "   " << (unitcell[d][0]  >= 0 ? " ":"") << unitcell[d][0] << " " << (unitcell[d][1]  >= 0 ? " ":"") << unitcell[d][1] << " " << 0.0 << " \n";
		}
	outFile << "    " << 0.0 << "  " << 0.0 << "  " << 100.0 << " \n";
	outFile << "end unit_cell_cart \n";
	outFile << "\n";

	outFile << "begin kpoints \n";
	for (int k = 0; k < num_k ; ++k){
			std::vector<double> k_vec = results[k].getDoubleVec("k_vec_grid");
			outFile << "   " << (k_vec[0]  >= 0 ? " ":"") << k_vec[0] << " " << (k_vec[1]  >= 0 ? " ":"") << k_vec[1] << "  " << 0.0 << " \n";
	}
	outFile << "end kpoints \n";
	outFile << " \n";

	/*
	// generate list of neighbors
	int nntot = 8; // number of nearest neighbors for a point in the k-grid

	// generate k grid
	int k_idx = 0;
	std::vector< std::vector<int> > ij_to_k;
	std::vector< std::vector<int> > k_to_ij;
	ij_to_k.resize(num_k1);
	k_to_ij.resize(num_k);
	for (int i = 0; i < num_k1; ++i){
		ij_to_k[i].resize(num_k2);
		for (int j = 0; j < num_k2; ++j){
			ij_to_k[i][j] = k_idx;
			k_to_ij[k_idx].resize(2);
			k_to_ij[k_idx][0] = i;
			k_to_ij[k_idx][1] = j;
			k_idx++;
		}
	}


	outFile << "begin nnkpts \n";
	outFile << "   " << nntot << " \n";
	for (int k = 0; k < num_k; ++k){

		int i0 = k_to_ij[k][0];
		int j0 = k_to_ij[k][1];

		int iplus = i0 + 1;
		int iminus = i0 - 1;
		int jplus = j0 + 1;
		int jminus = j0 - 1;

		if (iplus >= num_k1)
			iplus = iplus - num_k1;
		if (iminus < 0)
			iminus = iminus + num_k1;
		if (jplus >= num_k2)
			jplus = jplus - num_k2;
		if (jminus < 0)
			jminus = jminus + num_k2;

		outFile << "   " << k+1 << "   " << ij_to_k[iminus][j0] << "   0   0   0   \n";
		outFile << "   " << k+1 << "   " << ij_to_k[iplus][j0]<< "   0   0   0   \n";
		outFile << "   " << k+1 << "   " << ij_to_k[i0][jminus] << "   0   0   0   \n";
		outFile << "   " << k+1 << "   " << ij_to_k[i0][jplus] << "   0   0   0   \n";
		outFile << "   " << k+1 << "   " << ij_to_k[iminus][jminus] << "   0   0   0   \n";
		outFile << "   " << k+1 << "   " << ij_to_k[iplus][jplus] << "   0   0   0   \n";
		outFile << "   " << k+1 << "   " << k+1 << "   0   0   1   \n";
		outFile << "   " << k+1 << "   " << k+1 << "   0   0  -1   \n";
	}
	outFile << "end nnkpts \n";
	outFile << " \n";
	*/


}

void Param_tools::wannier_eig_save(std::vector<Job_params>& results, std::ofstream& outFile){

	outFile << fixed;
	outFile.precision(12);

	int num_k1 = results[0].getInt("num_k1");
	int num_k2 = results[0].getInt("num_k2");
	int num_k = num_k1*num_k2;
	std::vector<double> eigenvalues = results[0].getDoubleVec("eigenvalues");
	int num_bands = eigenvalues.size();

	for (int k = 0; k < num_k; ++k){
		eigenvalues = results[k].getDoubleVec("eigenvalues");
		for (int b = 0; b < num_bands; ++b){
			outFile << b+1 << " " << k+1 << " " << (eigenvalues[b] >= 0 ? " ":"") << eigenvalues[b] << " \n";
		}
	}

}

void Param_tools::wannier_mmn_save(std::vector<Job_params>& results, std::ofstream& outFile){

	outFile << fixed;
	outFile.precision(12);
	int num_k1 = results[0].getInt("num_k1");
	int num_k2 = results[0].getInt("num_k2");
	int num_k = num_k1*num_k2;
	std::vector< std::vector<double> > vecs = results[0].getDoubleMat("eigenvectors_r");
	int num_bands = vecs.size();
	int num_orbs = vecs[0].size();

	int nntot = 6; // number of nearest neighbors for a point in the k-grid

	outFile << "# *.mmn file generated by LEO2D \n";							// +2 for z-wrapping
	outFile << "   " << num_bands << " " << num_k << " " << nntot /* +2 */ << " \n";
	// generate k grid
	int k_idx = 0;
	std::vector< std::vector<int> > ij_to_k;
	std::vector< std::vector<int> > k_to_ij;
	ij_to_k.resize(num_k1);
	k_to_ij.resize(num_k);
	for (int i = 0; i < num_k1; ++i){
		ij_to_k[i].resize(num_k2);
		for (int j = 0; j < num_k2; ++j){
			ij_to_k[i][j] = k_idx;
			k_to_ij[k_idx].resize(2);
			k_to_ij[k_idx][0] = i;
			k_to_ij[k_idx][1] = j;
			k_idx++;
		}
	}

	// generate list of neighbors
	std::vector< std::vector<int> > neighbors;
	std::vector< std::vector<int> > sc_dx;
	std::vector< std::vector<int> >  sc_dy;
	neighbors.resize(num_k);
	sc_dx.resize(num_k);
	sc_dy.resize(num_k);

	for (int k = 0; k < num_k; ++k){

		int i0 = k_to_ij[k][0];
		int j0 = k_to_ij[k][1];

		int iplus = i0 + 1;
		int iminus = i0 - 1;
		int jplus = j0 + 1;
		int jminus = j0 - 1;

		int sc_dx_p = 0;
		int sc_dx_m = 0;
		int sc_dy_p = 0;
		int sc_dy_m = 0;

		if (iplus >= num_k1){
			iplus = iplus - num_k1;
			sc_dx_p = 1;
		}

		if (iminus < 0){
			iminus = iminus + num_k1;
			sc_dx_m = -1;
		}

		if (jplus >= num_k2){
			jplus = jplus - num_k2;
			sc_dy_p = 1;
		}

		if (jminus < 0){
			jminus = jminus + num_k2;
			sc_dy_m = -1;
		}
		neighbors[k].resize(6);
		neighbors[k][0] = ij_to_k[iminus][j0];
		neighbors[k][1] = ij_to_k[iplus][j0];
		neighbors[k][2] = ij_to_k[i0][jminus];
		neighbors[k][3] = ij_to_k[i0][jplus];
		neighbors[k][4] = ij_to_k[iminus][jminus];
		neighbors[k][5] = ij_to_k[iplus][jplus];

		sc_dx[k].resize(6);
		sc_dx[k][0] = sc_dx_m;
		sc_dx[k][1] = sc_dx_p;
		sc_dx[k][2] = 0;
		sc_dx[k][3] = 0;
		sc_dx[k][4] = sc_dx_m;
		sc_dx[k][5] = sc_dx_p;

		sc_dy[k].resize(6);
		sc_dy[k][0] = 0;
		sc_dy[k][1] = 0;
		sc_dy[k][2] = sc_dy_m;
		sc_dy[k][3] = sc_dy_p;
		sc_dy[k][4] = sc_dy_m;
		sc_dy[k][5] = sc_dy_p;

	}

	for (int k = 0; k < num_k; ++k){
		std::vector< std::vector<double> > vecs_kr = results[k].getDoubleMat("eigenvectors_r");
		std::vector< std::vector<double> > vecs_kc = results[k].getDoubleMat("eigenvectors_c");
		printf("k = %d \n",k);
		for (int nn = 0; nn < nntot; ++nn){
			int new_k = neighbors[k][nn];
			std::vector< std::vector<double> > vecs_newkr = results[new_k].getDoubleMat("eigenvectors_r");
			std::vector< std::vector<double> > vecs_newkc = results[new_k].getDoubleMat("eigenvectors_c");
			outFile << "   " << k+1 << " " << new_k+1 << " " << sc_dx[k][nn] << " " << sc_dy[k][nn] << " 0 \n";

			for (int i = 0; i < num_bands; ++i){
				for (int j = 0; j < num_bands; ++j){
					std::complex<double> temp_sum (0.,0.);
					for (int o = 0; o < num_orbs; ++o){
						temp_sum = temp_sum + std::conj(
																		std::complex<double>(vecs_kr[i][o],vecs_kc[i][o])
																	)*std::complex<double>(vecs_newkr[j][o],vecs_newkc[j][o]);
					}
					outFile << "      " << (temp_sum.real()  >= 0 ? " ":"") << temp_sum.real() << " " << (temp_sum.imag()  >= 0 ? " ":"") << temp_sum.imag() << " \n";
				}
			}

		}

		// For doing Z wrapping for wannier90
		/*
		outFile << "   " << k+1 << "   " << k+1 << "   0   0   1   \n";
		for (int i = 0; i < num_bands; ++i){
			for (int j = 0; j < num_bands; ++j){
				outFile << "      " << 0. <<  " " << 0. << " \n";
			}
		}
		outFile << "   " << k+1 << "   " << k+1 << "   0   0  -1   \n";
		for (int i = 0; i < num_bands; ++i){
			for (int j = 0; j < num_bands; ++j){
				outFile << "      " << 0. <<  " " << 0. << " \n";
			}
		}
		*/

	}

}
