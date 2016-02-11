
/* 
 * File:   hstruct.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:43 PM
 */

#ifndef HSTRUCT_H
#define HSTRUCT_H

#include "sheet.h"

class Hstruct {
    private:
        int max_index;
        int max_sheets;
        std::vector<Sheet> sheets;
        std::vector<double> angles;
        std::vector<double> heights;
        std::vector<std::vector<double> > shifts;
		std::vector<std::vector<int> > index_array;
        void setIndex();
        
    public:
        Hstruct(std::vector<Sheet>,std::vector<double>,std::vector<double>);
        Hstruct(const Hstruct& orig);
        ~Hstruct();
        void setShift(int, std::vector<double>);
        double posAtomIndex(int, int);
        double posAtomGrid(int (&grid_index)[3], int, int);
        int findNearest(double (&pos)[3], int, int);
        double interHamiltonian();
        double totalHamiltonian();
        double localitySweep(int, double, double);
        double hUpdate();
        int getMaxIndex();
        std::vector<std::vector<int> > getInterPairs();
		void getIntraPairs(std::vector<int>&,std::vector<int>&,std::vector<double>&);
		std::vector<std::vector<int> > getIndexArray();
		int gridToIndex(int (&grid_index)[4]);
		void getIndexToPos(double*,int);

};

#endif /* HSTRUCT_H */
