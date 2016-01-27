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


class Sheet {
    private:
        std::vector<std::vector<double> > a;
        std::vector<int> max_shape, min_shape;
        int max_index;
        std::vector<int> atom_types;
        std::vector<std::vector<double> > atom_pos;
        std::vector<std::vector<std::vector<int> > > grid_array;
        std::vector<std::vector<int> > index_array;
        void setIndex();
        
    public:
        Sheet(std::vector<std::vector<double> >, std::vector<int>, std::vector<std::vector<double> >, std::vector<int>, std::vector<int>);
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
        

};


#endif /* SHEET_H */

