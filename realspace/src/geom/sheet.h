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
#include <sstream>
#include <fstream>
#include <cstdlib>
#include "params/job_params.h"
#include "geom/sdata.h"
#include "materials/materials.h"



class Sheet {
    private:
      std::vector<std::vector<double> > a;
  	  std::vector<std::vector<double> > b;
      std::vector<int> max_shape, min_shape;
      double max_R;
      int max_index;
  	  Materials::Mat mat;
      int n_orbitals;
      std::vector<std::vector<double> > atom_pos;
      std::vector<std::vector<std::vector<int> > > grid_array;
      std::vector<std::vector<int> > index_array;
  	  std::vector<std::vector<double> > pos_array;
      std::vector<std::vector<double> > supercell_pair_array;
  	  bool ranSetup;
      void setIndex();
  	  void loadIndexRealspace();
  	  void loadIndexConfiguration();
  	  void setInverse();
  	  std::vector< std::vector<double> > a_inverse;
  	  int boundary_condition;
  	  int solver_space;
  	  int strain_type;
  	  std::string strain_file;
      std::vector< std::vector<double> > supercell;
      std::vector< std::vector<double> > supercell_inv;
      std::vector< std::vector<int> > supercell_stride;


    public:
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
  		double getOrbPos(int, int);
      int getNumAtoms();
  		Materials::Mat getMat();
  		void getIntraPairs(std::vector<int>&, std::vector<int>&, std::vector<double>&, std::vector< std::vector<int> >&, Job_params, int);
      void orderPairs(std::vector<int>&, std::vector<int>&, std::vector<double>&, std::vector< std::vector<int> >&);
      int findNearest(double (&pos)[3], int);

      void setReciprocal();
      void setSupercell(std::vector< std::vector<double> >);
  		double crossProd(std::vector<double>, std::vector<double>, int);
  		double getReciprocal(int, int);
  		double getInverse(int, int);

};


#endif /* SHEET_H */
