/* 
 * File:   mpi_job_params.h
 * Author: Stephen
 * 
 * Created on September 15, 2016, 2:52 PM
 */

#include "mpi_job_params.h"

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

Mpi_job_params::Mpi_job_params() {

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
	diagonalize = 0;
	d_weights = 1;
	d_vecs = 0;
	d_cond = 0;
	
	poly_order = 20;
	
	magOn = 0;
	elecOn = 0;
	B = 0;
	E = 0;
	
}

Mpi_job_params::Mpi_job_params(const Mpi_job_params& orig){

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
	
}

Mpi_job_params::~Mpi_job_params(){

}

void Mpi_job_params::loadLocParams(Loc_params opts){

	energy_rescale = opts.getDouble("energy_rescale");
	energy_shift = opts.getDouble("energy_shift");
	vacancy_chance = opts.getDouble("vacancy_chance");

	solver_type = opts.getInt("solver_type");
	observable_type = opts.getInt("observable_type");
	solver_space = opts.getInt("solver_space");
	diagonalize = opts.getInt("diagonalize");
	d_weights = opts.getInt("d_weights");
	d_vecs = opts.getInt("d_vecs");
	d_cond = opts.getInt("d_cond");	
	
	num_target_sheets = opts.getInt("num_target_sheets");
	target_sheets = new int[num_target_sheets];
	std::vector<int> opts_target_sheets = opts.getVecInt("target_sheets");
	
	for (int i = 0; i < num_target_sheets; ++i){
		target_sheets[i] = opts_target_sheets[i];
	}
	
	poly_order = opts.getInt("poly_order");
	
	magOn = opts.getInt("magOn");
	elecOn = opts.getInt("elecOn");
	B = opts.getDouble("B");
	E = opts.getDouble("E");

}

void Mpi_job_params::setParam(std::string tag, int val){

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
}

void Mpi_job_params::setParam(std::string tag, double val){

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

void Mpi_job_params::setParam(std::string tag, int *val, int dim){
	
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
		printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());
		}
}

void Mpi_job_params::setParam(std::string tag, double *val, int dim){

}

void Mpi_job_params::setParam(std::string tag, int *val, int dim1, int dim2){

}

void Mpi_job_params::setParam(std::string tag, double *val, int dim1, int dim2){

	if (tag == "shifts") {
		
		if (dim2 != 3){
			printf("WARNING: (From: Mpi_job_params) cannot set shift parameter with dim2 != 3! \n"); 
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
		printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());
	}
}

int Mpi_job_params::getInt(std::string tag) const{

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
		
	printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());
}

double Mpi_job_params::getDouble(std::string tag) const{

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
		
	printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());

}

int* Mpi_job_params::getIntVec(std::string tag) const{
	
	if (tag == "target_sheets")
		return target_sheets;
	if (tag == "target_list")
		return target_list;
	if (tag == "vacancy_list")
		return vacancy_list;
	
	printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());

}

double* Mpi_job_params::getDoubleVec(std::string tag) const{

	printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());
}


int* Mpi_job_params::getIntMat(std::string tag) const{

	printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());
}

double* Mpi_job_params::getDoubleMat(std::string tag) const{
	
	if (tag == "shifts")
		return shifts;
		
	printf("WARNING: Mpi_job_params variable <%s> not found. \n", tag.c_str());

}

void Mpi_job_params::printParams(){

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

void Mpi_job_params::printCheb(std::ofstream& outFile){

	outFile << "JOBID = " << jobID << " \n";
	
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
	
	outFile << "MAG_ON  = " << magOn  << ", B = " << B << " \n";
	outFile << "ELEC_ON = " << elecOn << ", E = " << E << " \n";
	
	
}

void Mpi_job_params::sendParams(int target, int tag){


		sendInt(solver_type,target,tag);
		
		if (solver_type != -1) {
		
			sendInt(observable_type,target,tag);
			sendInt(solver_space,target,tag);
			sendInt(diagonalize,target,tag);
			sendInt(d_weights,target,tag);
			sendInt(d_vecs,target,tag);
			sendInt(d_cond,target,tag);

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
			
		}
	
}


void Mpi_job_params::recvParams(int from){
		
		recvInt(from,"solver_type");
		
		if (solver_type != -1) {
		
			recvInt(from,"observable_type");
			recvInt(from,"solver_space");
			recvInt(from,"diagonalize");
			recvInt(from,"d_weights");
			recvInt(from,"d_vecs");
			recvInt(from,"d_cond");
		
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
		}
}

void Mpi_job_params::sendInt(int val, int target, int tag){

	MPI::COMM_WORLD.Send(	
				&val, 				// input buffer
				1,					// size of buffer 
				MPI::INT,			// type of buffer
				target,				// rank to receive
				tag);				// tag to label
}

void Mpi_job_params::sendDouble(double val, int target, int tag){

	MPI::COMM_WORLD.Send(	
				&val, 				// input buffer
				1,					// size of buffer 
				MPI::DOUBLE,		// type of buffer
				target,				// rank to receive
				tag);				// tag to label
}

void Mpi_job_params::sendIntVec(int* val, int dim, int target, int tag){

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

void Mpi_job_params::sendDoubleVec(double* val, int dim, int target, int tag){

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


void Mpi_job_params::sendIntMat(int* val, int dim1, int dim2, int target, int tag){

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

void Mpi_job_params::sendDoubleMat(double* val, int dim1, int dim2, int target, int tag){

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

void Mpi_job_params::recvInt(int from, std::string var){

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

void Mpi_job_params::recvDouble(int from, std::string var){

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

void Mpi_job_params::recvIntVec(int from, std::string var){

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

void Mpi_job_params::recvDoubleVec(int from, std::string var){

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

void Mpi_job_params::recvIntMat(int from, std::string var){

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

void Mpi_job_params::recvDoubleMat(int from, std::string var){

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
