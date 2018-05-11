/*
 * File:   job_params.cpp
 * Author: Stephen
 *
 * Created on August 24, 2016, 2:46 PM
 */

#include "params/job_params.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <mpi.h>

// --------------------------------------
// Old constructor, no use of class Sdata
// --------------------------------------
Job_params::Job_params() {

	// Set some default type defs here to avoid problems with incorrect casting
	int i_one = 1;
	int i_zero = 0;
	double d_one = 1.0;
	double d_zero = 0.0;

	// A default line-cut setting
	std::vector< std::vector<double> > lc_points;
	lc_points.resize(2);
	lc_points[0].resize(2);
	lc_points[1].resize(2);
	lc_points[0][0] = 0;
	lc_points[0][1] = 0;
	lc_points[1][0] = 1.0;
	lc_points[1][1] = 1.0;

	// Default job name
	setParam("job_name","LEO2D_JOB");

	// MPI sweep (job array) settings
	setParam("solver_type",i_zero);
	setParam("observable_type",i_zero);
	setParam("solver_space",i_zero);
	setParam("nShifts",i_one);
	setParam("num_shift_sheets",i_one);
	setParam("num_lc_points",i_one*2);
	setParam("lc_points",lc_points);
	setParam("cond_poly_par",0);

	// Momentum-space settings
	setParam("mom_vf_only",i_zero);

	// Geometery settings
	setParam("num_sheets",i_one);
	setParam("boundary_condition",i_zero);
	setParam("sc_search_size", i_one);
	setParam("global_shifts_on",i_zero);

	// Strain settings
	setParam("strain_type",i_zero);
	setParam("gsfe_z_on",i_zero);

	// Debug settings
	setParam("matrix_save",i_zero);
	setParam("matrix_pos_save", i_zero);
	setParam("verbose_save",i_one);

	// Diagonalization settings
	setParam("diagonalize",i_zero);
	setParam("d_type",i_zero);
	setParam("chiral_on",i_zero);
	setParam("d_kpm_dos",i_zero);
	setParam("d_weights",i_zero);
	setParam("d_vecs",i_zero);
	setParam("d_cond",i_zero);

	// K-sampling settngs (for periodic BCs)
	setParam("k_sampling",i_zero);
	setParam("k_type",i_zero);

	// settings for the FFT for momentum-space interlayer coupling
	setParam("fft_from_file",i_zero);
	setParam("fft_n_x",40*i_one);
	setParam("fft_n_y",40*i_one);
	setParam("fft_L_x",30*i_one);
	setParam("fft_L_y",30*i_one);
	setParam("fft_length_x",2*i_one);
	setParam("fft_length_y",2*i_one);
	setParam("fft_file","interlayer_fft.dat");

	// Hamiltonian pairing settings
	setParam("intra_searchsize",5*i_one);
	setParam("inter_searchsize",5*i_one);

	// KPM settings
	setParam("energy_rescale",20.0*d_one);
	setParam("energy_shift",d_zero);
	setParam("num_target_sheets",i_one);
	setParam("num_targets",i_zero);
	setParam("poly_order",20*i_one);
	setParam("dos_transform",i_one);
	setParam("cond_transform",i_one);

	// E&B Field settings
	setParam("magOn",i_zero);
	setParam("elecOn",i_zero);
	setParam("B",d_zero);
	setParam("E",d_zero);

	// MLMC settings
	setParam("mlmc",i_zero);
	setParam("mlmc_max_level",i_one);
	setParam("mlmc_level",i_zero);
	setParam("mlmc_num_clusters",i_zero);
	setParam("mlmc_cluster_size",4*i_one);
	setParam("mlmc_out_root","out");
	setParam("mlmc_temp_root","temp");
	setParam("mlmc_jobID",-1*i_one);
	setParam("mlmc_clusterID",-1*i_one);

	// Wannier settings
	setParam("wan_save",i_zero);
	setParam("wan_num_bands",-1*i_one);

	// Ballistic Transport
	setParam("ballistic_transport",i_zero);
	setParam("ballistic_sigma",10.0);
	setParam("ballistic_time_step",1.0);
	setParam("ballistic_max_steps",10*i_one);

	// Vacancy settings
	setParam("vacancy_file","vacancies.dat");
	setParam("num_vacancies",i_zero);
	setParam("vacancy_chance",d_zero);

}

Job_params::Job_params(const Job_params& orig){

	// loop over int
	std::vector<string> temp_tags = orig.getParamTags("int");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getInt(tag));
	}

	// loop over double
	temp_tags = orig.getParamTags("double");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getDouble(tag));
	}

	// loop over string
	temp_tags = orig.getParamTags("string");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getString(tag));
	}

	// loop over intVec
	temp_tags = orig.getParamTags("intVec");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getIntVec(tag));
	}

	// loop over doubleVec
	temp_tags = orig.getParamTags("doubleVec");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getDoubleVec(tag));
	}

	// loop over stringVec
	temp_tags = orig.getParamTags("stringVec");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getStringVec(tag));
	}

	// loop over intMat
	temp_tags = orig.getParamTags("intMat");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getIntMat(tag));
	}

	// loop over doubleMat
	temp_tags = orig.getParamTags("doubleMat");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getDoubleMat(tag));
	}

	// loop over cpxDoubleMat
	temp_tags = orig.getParamTags("cpxDoubleMat");
	for (int i = 0; i < temp_tags.size(); ++i){
		string tag = temp_tags[i];
		setParam(tag, orig.getCpxDoubleMat(tag));
	}
}

Job_params::~Job_params(){

}

void Job_params::setParam(std::string tag, int val){

	// first check if this tag already exists
	for (int i = 0; i < int_param_tags.size(); ++i){
		if (int_param_tags[i].compare(tag) == 0) {
			int_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	int_params.push_back(val);
	int_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, double val){

	// first check if this tag already exists

	for (int i = 0; i < double_param_tags.size(); ++i){
		if (double_param_tags[i].compare(tag) == 0) {
			double_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	double_params.push_back(val);
	double_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, string val){

	// first check if this tag already exists

	for (int i = 0; i < string_param_tags.size(); ++i){
		if (string_param_tags[i].compare(tag) == 0) {
			string_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	string_params.push_back(val);
	string_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, vector<int> val){

	// first check if this tag already exists

	for (int i = 0; i < int_vec_param_tags.size(); ++i){
		if (int_vec_param_tags[i].compare(tag) == 0) {
			int_vec_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	int_vec_params.push_back(val);
	int_vec_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, vector<double> val){

	// first check if this tag already exists

	for (int i = 0; i < double_vec_param_tags.size(); ++i){
		if (double_vec_param_tags[i].compare(tag) == 0) {
			double_vec_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	double_vec_params.push_back(val);
	double_vec_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, vector<string> val){

	// first check if this tag already exists

	for (int i = 0; i < string_vec_param_tags.size(); ++i){
		if (string_vec_param_tags[i].compare(tag) == 0) {
			string_vec_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	string_vec_params.push_back(val);
	string_vec_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, vector< vector<int> > val){

	// first check if this tag already exists

	for (int i = 0; i < int_mat_param_tags.size(); ++i){
		if (int_mat_param_tags[i].compare(tag) == 0) {
			int_mat_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	int_mat_params.push_back(val);
	int_mat_param_tags.push_back(tag);

}

void Job_params::setParam(std::string tag, vector< vector<double> > val){

	// first check if this tag already exists

	for (int i = 0; i < double_mat_param_tags.size(); ++i){
		if (double_mat_param_tags[i].compare(tag) == 0) {
			double_mat_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	double_mat_params.push_back(val);
	double_mat_param_tags.push_back(tag);

}


void Job_params::setParam(std::string tag, vector< vector< complex<double> > > val){

	// first check if this tag already exists

	for (int i = 0; i < cpx_double_mat_param_tags.size(); ++i){
		if (cpx_double_mat_param_tags[i].compare(tag) == 0) {
			cpx_double_mat_params[i] = val;
			return;
		}
	}

	// if not, add it to the list of params

	cpx_double_mat_params.push_back(val);
	cpx_double_mat_param_tags.push_back(tag);

}

int Job_params::getInt(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < int_param_tags.size(); ++i){
		if (int_param_tags[i].compare(tag) == 0) {
			return int_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getInt tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

double Job_params::getDouble(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < double_param_tags.size(); ++i){
		if (double_param_tags[i].compare(tag) == 0) {
			return double_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getDouble tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

string Job_params::getString(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < string_param_tags.size(); ++i){
		if (string_param_tags[i].compare(tag) == 0) {
			return string_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getString tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

vector<int> Job_params::getIntVec(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < int_vec_param_tags.size(); ++i){
		if (int_vec_param_tags[i].compare(tag) == 0) {
			return int_vec_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getIntVec tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

vector<double> Job_params::getDoubleVec(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < double_vec_param_tags.size(); ++i){
		if (double_vec_param_tags[i].compare(tag) == 0) {
			return double_vec_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getDoubleVec tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

vector<string> Job_params::getStringVec(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < string_vec_param_tags.size(); ++i){
		if (string_vec_param_tags[i].compare(tag) == 0) {
			return string_vec_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getStringVec tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

vector< vector<int> > Job_params::getIntMat(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < int_mat_param_tags.size(); ++i){
		if (int_mat_param_tags[i].compare(tag) == 0) {
			return int_mat_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getIntMat tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
}

vector< vector<double> > Job_params::getDoubleMat(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < double_mat_param_tags.size(); ++i){
		if (double_mat_param_tags[i].compare(tag) == 0) {
			return double_mat_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getDoubleMat tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
}

vector< vector< complex<double> > > Job_params::getCpxDoubleMat(std::string tag) const{

	// check if param exists, then return
	for (int i = 0; i < cpx_double_mat_param_tags.size(); ++i){
		if (cpx_double_mat_param_tags[i].compare(tag) == 0) {
			return cpx_double_mat_params[i];
		}
	}

	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getCpxDoubleMat tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
}

vector<string> Job_params::getParamTags(std::string name) const{

	if (name.compare("int") == 0)
		return int_param_tags;
	else if (name.compare("double") == 0)
		return double_param_tags;
	else if (name.compare("string") == 0)
		return string_param_tags;
	else if (name.compare("intVec") == 0)
		return int_vec_param_tags;
	else if (name.compare("doubleVec") == 0)
		return double_vec_param_tags;
	else if (name.compare("stringVec") == 0)
		return string_vec_param_tags;
	else if (name.compare("intMat") == 0)
		return int_mat_param_tags;
	else if (name.compare("doubleMat") == 0)
		return double_mat_param_tags;
	else if (name.compare("cpxDoubleMat") == 0)
		return cpx_double_mat_param_tags;

	string error_str = "!!LEO2D Critical Error!!: getParamTags target '";
	error_str = error_str + name + "' not recognized by job_params object\n";
	throw invalid_argument(error_str);
}

void Job_params::printParams(){

	printf(" T Input parameters for Locality object: \n");
	for(int i = 0; i < int_param_tags.size(); ++i){
		printf(" | %s = %d\n",int_param_tags[i].c_str(),int_params[i]);
	}

}

int Job_params::recvSpool(){

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
	recvParams(from);
	return from;

}

void Job_params::sendParams(int target){

	// loop over int
	int int_size = int_params.size();

	MPI::COMM_WORLD.Send(
				&int_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < int_size; ++i){
		sendInt(int_param_tags[i],int_params[i], target);
	}

	// loop over double
	int double_size = double_params.size();

	MPI::COMM_WORLD.Send(
				&double_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < double_size; ++i){
		sendDouble(double_param_tags[i],double_params[i], target);
	}

	// loop over string
	int string_size = string_params.size();

	MPI::COMM_WORLD.Send(
				&string_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < string_size; ++i){
		sendString(string_param_tags[i],string_params[i], target);
	}

	// loop over intVec
	int int_vec_size = int_vec_params.size();

	MPI::COMM_WORLD.Send(
				&int_vec_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < int_vec_size; ++i){
		sendIntVec(int_vec_param_tags[i],int_vec_params[i], target);
	}

	// loop over doubleVec
	int double_vec_size = double_vec_params.size();

	MPI::COMM_WORLD.Send(
				&double_vec_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < double_vec_size; ++i){
		sendDoubleVec(double_vec_param_tags[i],double_vec_params[i], target);
	}

	// loop over stringVec
	int string_vec_size = string_vec_params.size();

	MPI::COMM_WORLD.Send(
				&string_vec_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < string_vec_size; ++i){
		sendStringVec(string_vec_param_tags[i],string_vec_params[i], target);
	}

	// loop over intMat
	int int_mat_size = int_mat_params.size();

	MPI::COMM_WORLD.Send(
				&int_mat_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < int_mat_size; ++i){
		sendIntMat(int_mat_param_tags[i],int_mat_params[i], target);
	}

	// loop over doubleMat
	int double_mat_size = double_mat_params.size();

	MPI::COMM_WORLD.Send(
				&double_mat_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < double_mat_size; ++i){
		sendDoubleMat(double_mat_param_tags[i],double_mat_params[i], target);
	}

	// loop over cpxDoubleMat
	int cpx_double_mat_size = cpx_double_mat_params.size();

	MPI::COMM_WORLD.Send(
				&cpx_double_mat_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	for (int i = 0; i < cpx_double_mat_size; ++i){
		sendCpxDoubleMat(cpx_double_mat_param_tags[i],cpx_double_mat_params[i], target);
	}

}

void Job_params::recvParams(int from){

	MPI::Status status;

	// loop over int
	int int_size;

	MPI::COMM_WORLD.Recv(
				&int_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < int_size; ++i){
		recvInt(from);
	}

	// loop over double
	int double_size;

	MPI::COMM_WORLD.Recv(
				&double_size, 		// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < double_size; ++i){
		recvDouble(from);
	}

	// loop over string
	int string_size;

	MPI::COMM_WORLD.Recv(
				&string_size, 		// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < string_size; ++i){
		recvString(from);
	}

	// loop over intVec
	int int_vec_size ;

	MPI::COMM_WORLD.Recv(
				&int_vec_size, 			// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < int_vec_size; ++i){
		recvIntVec(from);
	}

	// loop over doubleVec
	int double_vec_size;

	MPI::COMM_WORLD.Recv(
				&double_vec_size, 	// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < double_vec_size; ++i){
		recvDoubleVec(from);
	}

	// loop over stringVec
	int string_vec_size;

	MPI::COMM_WORLD.Recv(
				&string_vec_size, 	// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < string_vec_size; ++i){
		recvStringVec(from);
	}

	// loop over intMat
	int int_mat_size;

	MPI::COMM_WORLD.Recv(
				&int_mat_size, 		// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < int_mat_size; ++i){
		recvIntMat(from);
	}

	// loop over doubleMat
	int double_mat_size;

	MPI::COMM_WORLD.Recv(
				&double_mat_size, 	// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < double_mat_size; ++i){
		recvDoubleMat(from);
	}

	// loop over cpxDoubleMat
	int cpx_double_mat_size;

	MPI::COMM_WORLD.Recv(
				&cpx_double_mat_size, 	// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	for (int i = 0; i < cpx_double_mat_size; ++i){
		recvCpxDoubleMat(from);
	}

}

void Job_params::sendTag(string tag, int target){

	int size = (int) tag.length();

	MPI::COMM_WORLD.Send(
				&size, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	MPI::COMM_WORLD.Send(
				tag.c_str(),		// input buffer
				size,				// size of buffer
				MPI::CHAR,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

}

string Job_params::recvTag(int from){

	int size;
	MPI::Status status;

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

	std::string string_out(input,size);
	return string_out;

}

void Job_params::sendInt(string tag, int val, int target){

	sendTag(tag, target);

	MPI::COMM_WORLD.Send(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

}

void Job_params::recvInt(int from){

	string tag = recvTag(from);

	int val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	setParam(tag, val);

}

void Job_params::sendDouble(string tag, double val, int target){

	sendTag(tag, target);

	MPI::COMM_WORLD.Send(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::DOUBLE,		// type of buffer
				target,				// rank to receive
				0);					// MPI label

}

void Job_params::recvDouble(int from){

	string tag = recvTag(from);

	double val;
	MPI::Status status;

	MPI::COMM_WORLD.Recv(
				&val, 				// input buffer
				1,					// size of buffer
				MPI::DOUBLE,		// type of buffer
				from,				// must come from root
				MPI::ANY_TAG,		// either WORKTAG or STOPTAG
				status);			// keep MPI status information

	setParam(tag, val);

}

void Job_params::sendString(string tag, string val, int target){

	sendTag(tag, target);

	int size = (int) val.length();

	MPI::COMM_WORLD.Send(
				&size, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	MPI::COMM_WORLD.Send(
				val.c_str(),		// input buffer
				size,				// size of buffer
				MPI::CHAR,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

}

void Job_params::recvString(int from){

	string tag = recvTag(from);

	int size;
	MPI::Status status;

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

	std::string val(input,size);

	setParam(tag, val);

}

void Job_params::sendIntVec(string tag, vector<int> val, int target){

	sendTag(tag, target);

	int dim;
	dim = val.size();
	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	if (dim != 0){

		MPI::COMM_WORLD.Send(
					&val[0],		 	// input buffer
					dim,				// size of buffer
					MPI::INT,			// type of buffer
					target,				// rank to receive
					0);					// MPI label
	}

}

void Job_params::recvIntVec(int from){

	string tag = recvTag(from);

	int dim;
	std::vector<int> val;
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
				MPI::INT,			// type of buffer
				from,				// must come from "from"
				status.Get_tag(),	// should be the same tag?
				status);			// keep MPI status information

		setParam(tag,val);

	}
}

void Job_params::sendDoubleVec(string tag, vector<double> val, int target){

	sendTag(tag, target);

	int dim;
	dim = val.size();
	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	if (dim != 0){

		MPI::COMM_WORLD.Send(
					&val[0],		 	// input buffer
					dim,				// size of buffer
					MPI::DOUBLE,		// type of buffer
					target,				// rank to receive
					0);					// MPI label
	}

}

void Job_params::recvDoubleVec(int from){

	string tag = recvTag(from);

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

		setParam(tag,val);

	}
}

void Job_params::sendStringVec(string tag, vector<string> val, int target){

	sendTag(tag, target);

	int dim;
	dim = val.size();
	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	if (dim != 0){
		for (int i = 0; i < dim; ++i){

			int size = (int) val[i].length();

			MPI::COMM_WORLD.Send(
						&size, 				// input buffer
						1,					// size of buffer
						MPI::INT,			// type of buffer
						target,				// rank to receive
						0);					// MPI label

			MPI::COMM_WORLD.Send(
						val[i].c_str(),		// input buffer
						size,				// size of buffer
						MPI::CHAR,			// type of buffer
						target,				// rank to receive
						0);					// MPI label

		}
	}

}

void Job_params::recvStringVec(int from){

	string tag = recvTag(from);

	int dim;
	std::vector<string> val;
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

		for (int i = 0; i < dim; ++i){

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

			std::string temp_string(input,size);
			val[i] = temp_string;

		}

		setParam(tag,val);

	}
}

void Job_params::sendIntMat(string tag, vector< vector<int> > val, int target){

	sendTag(tag, target);

	int dim1;
	int dim2;

	dim1 = val.size();

	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim1, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	if (dim1 != 0){

		int dim2_array[dim1];
		int total_dim = 0;

		for (int i = 0; i < dim1; ++i){
			dim2_array[i] = (int) val[i].size();
			total_dim += dim2_array[i];
		}

		MPI::COMM_WORLD.Send(
					&dim2_array, 		// input buffer
					dim1,				// size of buffer
					MPI::INT,			// type of buffer
					target,				// rank to receive
					0);					// MPI label

		int temp_val[total_dim];

		int counter = 0;
		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2_array[i]; ++j){
				temp_val[counter] = val[i][j];
				++counter;
			}
		}

		MPI::COMM_WORLD.Send(
					&temp_val,		 	// input buffer
					total_dim,			// size of buffer
					MPI::INT,			// type of buffer
					target,				// rank to receive
					0);					// MPI label
	}

}

void Job_params::recvIntMat(int from){

	string tag = recvTag(from);

	int dim1;

	std::vector< std::vector<int> > val;
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

		int dim2_array[dim1];

		MPI::COMM_WORLD.Recv(
					&dim2_array, 		// input buffer
					dim1,				// size of buffer
					MPI::INT,			// type of buffer
					from,				// must come from "from"
					MPI::ANY_TAG,		// any tag
					status);			// keep MPI status information

		int total_dim = 0;

		for (int i = 0; i < dim1; ++i){
			val[i].resize(dim2_array[i]);
			total_dim += dim2_array[i];
		}

		int temp_val[total_dim];

		MPI::COMM_WORLD.Recv(
				temp_val,		 	// input buffer
				total_dim,			// size of buffer
				MPI::INT,		// type of buffer
				from,				// must come from "from"
				status.Get_tag(),	// should be the same tag?
				status);			// keep MPI status information

		int counter = 0;
		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2_array[i]; ++j){
				val[i][j] = temp_val[counter];
				++counter;
			}
		}

		setParam(tag,val);

	}

}

void Job_params::sendDoubleMat(string tag, vector< vector<double> > val, int target){

	sendTag(tag, target);

	int dim1;
	int dim2;

	dim1 = val.size();

	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim1, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	if (dim1 != 0){

		int dim2_array[dim1];
		int total_dim = 0;

		for (int i = 0; i < dim1; ++i){
			dim2_array[i] = (int) val[i].size();
			total_dim += dim2_array[i];
		}

		MPI::COMM_WORLD.Send(
					&dim2_array, 		// input buffer
					dim1,				// size of buffer
					MPI::INT,			// type of buffer
					target,				// rank to receive
					0);					// MPI label

		double temp_val[total_dim];

		int counter = 0;
		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2_array[i]; ++j){
				temp_val[counter] = val[i][j];
				++counter;
			}
		}

		MPI::COMM_WORLD.Send(
					&temp_val,		 	// input buffer
					total_dim,			// size of buffer
					MPI::DOUBLE,		// type of buffer
					target,				// rank to receive
					0);					// MPI label
	}

}

void Job_params::recvDoubleMat(int from){

	string tag = recvTag(from);

	int dim1;

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

		int dim2_array[dim1];

		MPI::COMM_WORLD.Recv(
					&dim2_array, 		// input buffer
					dim1,				// size of buffer
					MPI::INT,			// type of buffer
					from,				// must come from "from"
					MPI::ANY_TAG,		// any tag
					status);			// keep MPI status information

		int total_dim = 0;

		for (int i = 0; i < dim1; ++i){
			val[i].resize(dim2_array[i]);
			total_dim += dim2_array[i];
		}

		double temp_val[total_dim];

		MPI::COMM_WORLD.Recv(
				temp_val,		 	// input buffer
				total_dim,			// size of buffer
				MPI::DOUBLE,		// type of buffer
				from,				// must come from "from"
				status.Get_tag(),	// should be the same tag?
				status);			// keep MPI status information

		int counter = 0;
		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2_array[i]; ++j){
				val[i][j] = temp_val[counter];
				++counter;
			}
		}

		setParam(tag,val);

	}

}

void Job_params::sendCpxDoubleMat(string tag, vector< vector< complex<double> > > val, int target){

	sendTag(tag, target);

	int dim1;
	int dim2;

	dim1 = val.size();

	MPI::Status status;

	MPI::COMM_WORLD.Send(
				&dim1, 				// input buffer
				1,					// size of buffer
				MPI::INT,			// type of buffer
				target,				// rank to receive
				0);					// MPI label

	if (dim1 != 0){

		int dim2_array[dim1];
		int total_dim = 0;

		for (int i = 0; i < dim1; ++i){
			dim2_array[i] = (int) val[i].size();
			total_dim += dim2_array[i];
		}

		MPI::COMM_WORLD.Send(
					&dim2_array, 		// input buffer
					dim1,				// size of buffer
					MPI::INT,			// type of buffer
					target,				// rank to receive
					0);					// MPI label

		double temp_val[2*total_dim];

		int counter = 0;
		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2_array[i]; ++j){
				temp_val[counter] = val[i][j].real();
				++counter;
				temp_val[counter] = val[i][j].imag();
				++counter;
			}
		}

		MPI::COMM_WORLD.Send(
					&temp_val,		 	// input buffer
					2*total_dim,		// size of buffer
					MPI::DOUBLE,		// type of buffer
					target,				// rank to receive
					0);					// MPI label
	}

}

void Job_params::recvCpxDoubleMat(int from){

	string tag = recvTag(from);

	int dim1;

	std::vector< std::vector< complex<double> > > val;
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

		int dim2_array[dim1];

		MPI::COMM_WORLD.Recv(
					&dim2_array, 		// input buffer
					dim1,				// size of buffer
					MPI::INT,			// type of buffer
					from,				// must come from "from"
					MPI::ANY_TAG,		// any tag
					status);			// keep MPI status information

		int total_dim = 0;

		for (int i = 0; i < dim1; ++i){
			val[i].resize(dim2_array[i]);
			total_dim += dim2_array[i];
		}

		double temp_val[2*total_dim];

		MPI::COMM_WORLD.Recv(
				temp_val,		 	// input buffer
				2*total_dim,		// size of buffer
				MPI::DOUBLE,		// type of buffer
				from,				// must come from "from"
				status.Get_tag(),	// should be the same tag?
				status);			// keep MPI status information

		int counter = 0;
		for (int i = 0; i < dim1; ++i){
			for (int j = 0; j < dim2_array[i]; ++j){
				val[i][j].real(temp_val[counter]);
				++counter;
				val[i][j].imag(temp_val[counter]);
				++counter;
			}
		}

		setParam(tag,val);

	}

}
