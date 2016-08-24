/* 
 * File:   loc_params.cpp
 * Author: Stephen
 * 
 * Created on August 24, 2016, 2:46 PM
 */

#include "loc_params.h"

// --------------------------------------
// Old constructor, no use of class Sdata
// --------------------------------------
Loc_params::Loc_params() {
	
		nShifts = 1;
		solver_type = 0;
		intra_searchsize = 5;
		inter_searchsize = 5;
		num_target_sheets = 1;
		poly_order = 20;
		magOn = 0;
		elecOn = 0;
		
		energy_rescale = 20;
		energy_shift = 0;
		B = 0;
		E = 0;
		vacancy_chance = 0;
		
		job_name = "UNKNOWN_JOB";
		
		target_sheets.push_back(0);
}

Loc_params::Loc_params(const Loc_params& orig){

		nShifts = orig.getInt("nShifts");
		solver_type = orig.getInt("solver_type");
		intra_searchsize = orig.getInt("intra_searchsize");
		inter_searchsize = orig.getInt("inter_searchsize");
		num_target_sheets = orig.getInt("num_target_sheets");
		poly_order = orig.getInt("poly_order");
		magOn = orig.getInt("magOn");
		elecOn = orig.getInt("elecOn");
		
		energy_rescale = orig.getDouble("energy_rescale");
		energy_shift = orig.getDouble("energy_shift");
		B = orig.getDouble("B");
		E = orig.getDouble("E");
		vacancy_chance = orig.getDouble("vacancy_chance");
		
		job_name = orig.getString("job_name");
		
		target_sheets = orig.getVecInt("target_sheets");
}

Loc_params::~Loc_params(){

}

void Loc_params::setParam(std::string tag, int val){

	if (tag == "nShifts")
		nShifts = val;
	if (tag == "solver_type")
		solver_type = val;
	if (tag == "intra_searchsize")
		intra_searchsize = val;
	if (tag == "inter_searchsize")
		inter_searchsize = val;
	if (tag == "num_target_sheets")
		num_target_sheets = val;
	if (tag == "poly_order")
		poly_order = val;
	if (tag == "magOn")
		magOn = val;
	if (tag == "elecOn")
		elecOn = val;
		
}

void Loc_params::setParam(std::string tag, double val){

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

void Loc_params::setParam(std::string tag, std::string val){
	
	if (tag == "job_name")
		job_name = val;
	
}

void Loc_params::setParam(std::string tag, std::vector<int> val){

	if (tag == "target_sheets")
		target_sheets = val;
		
}

int Loc_params::getInt(std::string tag) const{

	if (tag == "nShifts")
		return nShifts;
	if (tag == "solver_type")
		return solver_type;
	if (tag == "intra_searchsize")
		return intra_searchsize;
	if (tag == "inter_searchsize")
		return inter_searchsize;
	if (tag == "num_target_sheets")
		return num_target_sheets;
	if (tag == "poly_order")
		return poly_order;
	if (tag == "magOn")
		return magOn;
	if (tag == "elecOn")
		return elecOn;

}

double Loc_params::getDouble(std::string tag) const{

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

}

std::string Loc_params::getString(std::string tag) const{
	
	if (tag == "job_name")
		return job_name;
	
}

std::vector<int> Loc_params::getVecInt(std::string tag) const{

	if (tag == "target_sheets")
		return target_sheets;

}

