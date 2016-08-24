/* 
 * File:   loc_params.h
 * Author: Stephen
 * 
 * Created on August 24, 2016, 2:46 PM
 */

#ifndef LOC_PARAMS_H
#define LOC_PARAMS_H

#include <string>
#include <vector>

class Loc_params {
    private:
		
		int nShifts;
		int solver_type;
		int intra_searchsize;
		int inter_searchsize;
		int num_target_sheets;
		int poly_order;
		int magOn;
		int elecOn;
		
		double energy_rescale;
		double energy_shift;
		double B;
		double E;
		double vacancy_chance;
		
		std::string job_name;
		
		std::vector<int> target_sheets;
	
	public:
        Loc_params();
        Loc_params(const Loc_params& orig);
        ~Loc_params();
		
        void setParam(std::string, int);
		void setParam(std::string, double);
		void setParam(std::string, std::string);
		void setParam(std::string, std::vector<int>);
		
		int getInt(std::string) const;
		double getDouble(std::string) const;
		std::string getString(std::string) const;
		std::vector<int> getVecInt(std::string) const;
};

#endif /* LOC_PARAMS_H */