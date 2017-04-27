/* 
 * File:   loc_params.cpp
 * Author: Stephen
 * 
 * Created on August 24, 2016, 2:46 PM
 */

#include "loc_params.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

// --------------------------------------
// Old constructor, no use of class Sdata
// --------------------------------------
Loc_params::Loc_params() {
	
		nShifts = 1;
		num_shift_sheets = 1;
		solver_type = 0;
		strain_type = 0;
		observable_type = 0;
		solver_space = 0;
		diagonalize = 0;
		d_weights = 1;
		d_vecs = 0;
		d_cond = 0;
		fft_from_file = 0;
		intra_searchsize = 5;
		inter_searchsize = 5;
		num_target_sheets = 1;
		poly_order = 20;
		magOn = 0;
		elecOn = 0;
		
		mlmc = 0;
		mlmc_max_level = 1;
		mlmc_level = 0;
		mlmc_num_clusters = 0;
		mlmc_cluster_size = 4;
		
		energy_rescale = 20.0;
		energy_shift = 0.0;
		B = 0.0;
		E = 0.0;
		vacancy_chance = 0.0;
		
		int* shift_sheets;
		
		job_name = "UNKNOWN_JOB";
		
		target_sheets.push_back(0);
}

Loc_params::Loc_params(const Loc_params& orig){

		nShifts = orig.getInt("nShifts");
		num_shift_sheets = orig.getInt("num_shift_sheets");
		solver_type = orig.getInt("solver_type");
		strain_type = orig.getInt("strain_type");
		observable_type = orig.getInt("observable_type");
		solver_space = orig.getInt("solver_space");
		diagonalize = orig.getInt("diagonalize");
		d_weights = orig.getInt("d_weights");
		d_vecs = orig.getInt("d_vecs");
		d_cond = orig.getInt("d_cond");
		fft_from_file = orig.getInt("fft_from_file");
		intra_searchsize = orig.getInt("intra_searchsize");
		inter_searchsize = orig.getInt("inter_searchsize");
		num_target_sheets = orig.getInt("num_target_sheets");
		poly_order = orig.getInt("poly_order");
		magOn = orig.getInt("magOn");
		elecOn = orig.getInt("elecOn");
		
		mlmc = orig.getInt("mlmc");
		mlmc_max_level = orig.getInt("mlmc_max_level");
		mlmc_level = orig.getInt("mlmc_level");
		mlmc_num_clusters = orig.getInt("mlmc_num_clusters");
		mlmc_cluster_size = orig.getInt("mlmc_cluster_size");
		
		energy_rescale = orig.getDouble("energy_rescale");
		energy_shift = orig.getDouble("energy_shift");
		B = orig.getDouble("B");
		E = orig.getDouble("E");
		vacancy_chance = orig.getDouble("vacancy_chance");
		
		int* temp_shift_sheets = orig.getIntVec("shift_sheets");
		shift_sheets = new int[num_shift_sheets];
		for (int s = 0; s < num_shift_sheets; ++s){
			shift_sheets[s] = temp_shift_sheets[s];
		}
		
		job_name = orig.getString("job_name");
		
		target_sheets = orig.getVecInt("target_sheets");
}

Loc_params::~Loc_params(){

}

void Loc_params::setParam(std::string tag, int val){

	if (tag == "nShifts")
		nShifts = val;
	if (tag == "num_shift_sheets")
		num_shift_sheets = val;
	if (tag == "solver_type")
		solver_type = val;
	if (tag == "strain_type")
		strain_type = val;
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
	if (tag == "fft_from_file")
		fft_from_file = val;
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
	if (tag == "mlmc")
		mlmc = val;
	if (tag == "mlmc_max_level")
		mlmc_max_level = val;		
	if (tag == "mlmc_level")
		mlmc_level = val;
	if (tag == "mlmc_num_clusters")
		mlmc_num_clusters = val;
	if (tag == "mlmc_cluster_size")
		mlmc_cluster_size = val;
		
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

void Loc_params::setParam(std::string tag, int* val){

	if (tag == "shift_sheets"){
		shift_sheets = new int[num_shift_sheets];
		for (int s = 0; s < num_shift_sheets; ++s){
			shift_sheets[s] = val[s];
		}
	}

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
	if (tag == "num_shift_sheets")
		return num_shift_sheets;
	if (tag == "solver_type")
		return solver_type;
	if (tag == "strain_type")
		return strain_type;
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
	if (tag == "fft_from_file")
		return fft_from_file;
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
	if (tag == "mlmc")
		return mlmc;
	if (tag == "mlmc_max_level")
		return mlmc_max_level ;		
	if (tag == "mlmc_level")
		return mlmc_level;
	if (tag == "mlmc_num_clusters")
		return mlmc_num_clusters;
	if (tag == "mlmc_cluster_size")
		return mlmc_cluster_size;

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

int* Loc_params::getIntVec(std::string tag) const{

	if (tag == "shift_sheets")
		return shift_sheets;

}

std::string Loc_params::getString(std::string tag) const{
	
	if (tag == "job_name")
		return job_name;
	
}

std::vector<int> Loc_params::getVecInt(std::string tag) const{

	if (tag == "target_sheets")
		return target_sheets;

}

void Loc_params::printParams(){

		printf(" T Input parameters for Locality object: \n");
		printf(" | nShifts = %d \n", nShifts);
		
		if (solver_type == 1)
			printf(" | solver_type = SQUARE \n");
		else if (solver_type == 2)
			printf(" | solver_type = LINECUT \n");
		else if (solver_type == 3)
			printf(" | solver_type = VACANCY DEFECTS ABOUT CENTER \n");
		else if (solver_type == 4)
			printf(" | solver_type = VACANDY DEFECTS FROM FILE \n");
		else if (solver_type == 5)
			printf(" | solver_type = STRAIN SAMPLING ABOUT CENTER \n");
			
		if (strain_type == 0)
			printf(" | strain_type = NONE \n");
		else if (solver_type == 1)
			printf(" | strain_type = REALSPACE FILE \n");
		else if (solver_type == 2)
			printf(" | strain_type = CONFIGURATION FILE \n");
	
	
		if (observable_type == 0)
			printf(" | observable_type = DOS \n");
		else if (observable_type == 1)
			printf(" | observable_type = CONDUCTIVITY \n");
		
		if (solver_space == 0)
			printf(" | solver_space = REALSPACE \n");
		else if (solver_space == 1)
			printf(" | solver_space = MOMENTUMSPACE \n");
		
		printf(" | diagonalize = %d \n", diagonalize);
		if (diagonalize == 1){
			printf(" | d_weights = %d \n",d_weights);
			printf(" | d_vecs = %d \n", d_vecs);
			printf(" | d_cond = %d \n", d_cond);
		}
		printf(" | poly_order = %d \n",poly_order);
		printf(" | intra_searchsize = %d \n", intra_searchsize);
		printf(" | inter_searchsize = %d \n", inter_searchsize);
		printf(" | num_target_sheets = %d \n", num_target_sheets);
		printf(" | poly_order = %d \n", poly_order);
		printf(" | magOn = %d \n", magOn);
		printf(" | elecOn = %d \n", elecOn);
		
		printf(" | energy_rescale = %lf \n", energy_rescale);
		printf(" | energy_shift = %lf \n", energy_shift);
		printf(" | B = %lf \n", B);
		printf(" | E = %lf \n", E);
		printf(" | vacancy_chance = %lf \n",vacancy_chance);
		
		std::cout << " | job_name = " << job_name << "\n";
		
		std::cout << " L target_sheets = [";
		std::cout << target_sheets[0] + 1;
		for (int s = 1; s < num_target_sheets; ++s){
			std::cout << ", " << target_sheets[s] + 1;
		}
		std::cout << "] \n";
}

