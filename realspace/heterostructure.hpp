/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   heterostructure.hpp
 * Author: Stephen
 *
 * Created on January 13, 2016, 3:16 PM
 */

#ifndef HETEROSTRUCTURE_H
#define HETEROSTRUCTURE_H
#include "sheet.hpp"

class Heterostructure {
    private:
        int max_index;
        int max_sheets;
        std::vector<Sheet> sheets;
        std::vector<double> angles;
        std::vector<double> heights;
        std::vector<std::vector<double> > shifts;
        void setIndex();
        
    public:
        Heterostructure(std::vector<Sheet>,std::vector<double>,std::vector<double>);
        Heterostructure(const Heterostructure& orig);
        ~Heterostructure();
        void setShift(int, std::vector<double>);
        double posAtomIndex(int, int);
        double posAtomGrid(int (&grid_index)[3], int, int);
        int findNearest(double (&pos)[3], int, int);
        double interHamiltonian();
        double totalHamiltonian();
        double localitySweep(int, double, double);
        double hUpdate();
        int getMaxIndex();
        

};

#endif /* HETEROSTRUCTURE_H */

