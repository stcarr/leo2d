
/*
 * File:   hstruct.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:43 PM
 */

#ifndef HSTRUCT_H
#define HSTRUCT_H

#include "geom/sheet.h"
#include "momentum/momentum_coupling.h"
#include "params/job_params.h"

class Hstruct {
    private:
      int max_index;
      int max_sheets;
  		int solver_space;

      std::vector<Sheet> sheets;
      std::vector<double> angles;
      std::vector<double> heights;
      std::vector<std::vector<double> > shifts;
      std::vector<std::vector<int> > index_array;

      void setIndex();

    public:
      Hstruct(std::vector<Sheet>,std::vector<double>,std::vector<double>,int);
      Hstruct(const Hstruct& orig);
      ~Hstruct();

      void setShift(int, std::vector<double>);
      int indexToSheet(int);

      double posAtomIndex(int, int);

      int findNearest(double (&pos)[3], int, int);

      double interHamiltonian();
      double totalHamiltonian();

      double localitySweep(int, double, double);

      double hUpdate();
      int getMaxIndex();

      void orderPairs(std::vector< std::vector<int> >&, std::vector< std::vector<double> >&);

      void getShiftConfigs(std::vector<std::vector<double> >&, Job_params);

      void getInterPairs(std::vector<std::vector<int> >&,std::vector<std::vector<double> >&,Job_params);
  		void getIntraPairs(std::vector<int>&,std::vector<int>&,std::vector<double>&,std::vector<std::vector<double> >&,Job_params);
  		std::vector<std::vector<int> > getIndexArray();
  		int gridToIndex(int (&grid_index)[4]);
  		void getIndexToPos(double*,int);
      std::vector<std::vector<double> > getSupercellVecs();
  		double getUnitArea(int);

  		std::vector<std::vector<int> > getVacancyList(int,int);
  		std::vector<std::vector<int> > getTargetList(Job_params);

  		void makeInterFFTFile(int, int, int, int, int, int, double, double, std::string);


};

#endif /* HSTRUCT_H */
