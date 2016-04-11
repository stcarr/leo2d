/* 
 * File:   sheet.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:42 PM
 */

#ifndef SHEET_H
#define SHEET_H
#include <vector>
#include <string>
#include "sdata.h"


class Sheet {
    private:
        std::vector<std::vector<double> > a;
        std::vector<int> max_shape, min_shape;
        int max_index;
        std::vector<int> atom_types;
        std::vector<std::vector<double> > atom_pos;
        std::vector<std::vector<std::vector<int> > > grid_array;
        std::vector<std::vector<int> > index_array;
		bool ranSetup;
        void setIndex();
		void setInverse();
		double a_inverse[2][2];
        
    public:
        Sheet(std::vector<std::vector<double> >, std::vector<int>, std::vector<std::vector<double> >, std::vector<int>, std::vector<int>);
		Sheet(Sdata);
        Sheet(const Sheet& orig);
        ~Sheet();
        bool checkShape(double (&pos)[3]);
        double posAtomIndex(int, int);
        double posAtomGrid(int (&grid_index)[3], int);
        int indexToGrid(int, int);
        int gridToIndex(int (&grid_index)[3]);
        double intraHamiltonian();
        int getMaxIndex();
        double getUnit(int, int);
		int getShape(int, int);
        int getNumAtoms();
		double getInverse(int, int);

};


#endif /* SHEET_H */

