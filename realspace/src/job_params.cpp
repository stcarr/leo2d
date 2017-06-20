/* 
 * File:   job_params.cpp
 * Author: Stephen
 * 
 * Created on August 24, 2016, 2:46 PM
 */

#include "job_params.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <mpi.h>

// --------------------------------------
// Old constructor, no use of class Sdata
// --------------------------------------
Job_params::Job_params() {
	
}

Job_params::Job_params(Job_params& orig){

	// loop over int
	int_param_tags = orig.getParamTags("int");
	for (int i = 0; i < int_param_tags.size(); ++i){
		string tag = int_param_tags[i];
		setParam(tag, orig.getInt(tag));
	}
	
	// loop over double
	double_param_tags = orig.getParamTags("double");
	for (int i = 0; i < double_param_tags.size(); ++i){
		string tag = double_param_tags[i];
		setParam(tag, orig.getDouble(tag));
	}
	
	// loop over string
	string_param_tags = orig.getParamTags("string");
	for (int i = 0; i < string_param_tags.size(); ++i){
		string tag = string_param_tags[i];
		setParam(tag, orig.getString(tag));
	}
	
	// loop over intVec
	int_vec_param_tags = orig.getParamTags("intVec");
	for (int i = 0; i < int_vec_param_tags.size(); ++i){
		string tag = int_vec_param_tags[i];
		setParam(tag, orig.getIntVec(tag));
	}
	
	// loop over doubleVec
	double_vec_param_tags = orig.getParamTags("doubleVec");
	for (int i = 0; i < double_vec_param_tags.size(); ++i){
		string tag = double_vec_param_tags[i];
		setParam(tag, orig.getDoubleVec(tag));
	}
	
	// loop over stringVec
	string_vec_param_tags = orig.getParamTags("stringVec");
	for (int i = 0; i < string_vec_param_tags.size(); ++i){
		string tag = string_vec_param_tags[i];
		setParam(tag, orig.getStringVec(tag));
	}
	
	// loop over intMat
	int_mat_param_tags = orig.getParamTags("intMat");
	for (int i = 0; i < int_mat_param_tags.size(); ++i){
		string tag = int_mat_param_tags[i];
		setParam(tag, orig.getIntMat(tag));
	}
	
	// loop over doubleMat
	double_mat_param_tags = orig.getParamTags("doubleMat");
	for (int i = 0; i < double_mat_param_tags.size(); ++i){
		string tag = double_mat_param_tags[i];
		setParam(tag, orig.getDoubleMat(tag));
	}

}

Job_params::~Job_params(){

}

void Job_params::setParam(std::string tag, int val){

	// first check if this tag already exists
	
	for (int i = 0; i < int_param_tags.size(); ++i){
		if (int_param_tags[i].compare(tag)) {
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
		if (double_param_tags[i].compare(tag)) {
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
		if (string_param_tags[i].compare(tag)) {
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
		if (int_vec_param_tags[i].compare(tag)) {
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
		if (double_vec_param_tags[i].compare(tag)) {
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
		if (string_vec_param_tags[i].compare(tag)) {
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
		if (int_mat_param_tags[i].compare(tag)) {
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
		if (double_mat_param_tags[i].compare(tag)) {
			double_mat_params[i] = val;
			return;
		}
	}
	
	// if not, add it to the list of params
	
	double_mat_params.push_back(val);
	double_mat_param_tags.push_back(tag);

}

int Job_params::getInt(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < int_param_tags.size(); ++i){
		if (int_param_tags[i].compare(tag)) {
			return int_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getInt tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
	
}

double Job_params::getDouble(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < double_param_tags.size(); ++i){
		if (double_param_tags[i].compare(tag)) {
			return double_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getDouble tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
	
}

string Job_params::getString(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < string_param_tags.size(); ++i){
		if (string_param_tags[i].compare(tag)) {
			return string_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getString tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
	
}

vector<int> Job_params::getIntVec(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < int_vec_param_tags.size(); ++i){
		if (int_vec_param_tags[i].compare(tag)) {
			return int_vec_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getIntVec tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
	
}

vector<double> Job_params::getDoubleVec(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < double_vec_param_tags.size(); ++i){
		if (double_vec_param_tags[i].compare(tag)) {
			return double_vec_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getDoubleVec tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);

}

vector<string> Job_params::getStringVec(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < string_vec_param_tags.size(); ++i){
		if (string_vec_param_tags[i].compare(tag)) {
			return string_vec_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getStringVec tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
	
}

vector< vector<int> > Job_params::getIntMat(std::string tag) {

	// check if param exists, then return
	for (int i = 0; i < int_mat_param_tags.size(); ++i){
		if (int_mat_param_tags[i].compare(tag)) {
			return int_mat_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getIntMat tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
}

vector< vector<double> > Job_params::getDoubleMat(std::string tag){

	// check if param exists, then return
	for (int i = 0; i < double_mat_param_tags.size(); ++i){
		if (double_mat_param_tags[i].compare(tag)) {
			return double_mat_params[i];
		}
	}
	
	// if not, throw error
	string error_str = "!!LEO2D Critical Error!!: getDoubleMat tag '";
	error_str = error_str + tag + "' not found in job_params object\n";
	throw invalid_argument(error_str);
}

vector<string> Job_params::getParamTags(std::string name){

	if (name.compare("int"))
		return int_param_tags;
	else if (name.compare("double"))
		return double_param_tags;
	else if (name.compare("string"))
		return string_param_tags;
	else if (name.compare("intVec"))
		return int_vec_param_tags;
	else if (name.compare("doubleVec"))
		return double_vec_param_tags;
	else if (name.compare("stringVec"))
		return string_vec_param_tags;
	else if (name.compare("intMat"))
		return int_mat_param_tags;
	else if (name.compare("doubleMat"))
		return double_mat_param_tags;
	
	string error_str = "!!LEO2D Critical Error!!: getParamTags target '";
	error_str = error_str + name + "' not recognized by job_params object\n";
	throw invalid_argument(error_str);
}

void Job_params::printParams(){

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
