
/*
 * File:   hstruct.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:43 PM
 */

#ifndef STRAIN_H
#define STRAIN_H

#include "params/job_params.h"

#include <string>
#include <vector>


class StrainCalc {
    private:
      int num_sheets;
      std::vector<int> num_orb;
      std::vector< std::vector<int>> grid;
      std::vector< std::vector< std::vector< std::vector<double> > > > disp_x;
      std::vector< std::vector< std::vector< std::vector<double> > > > disp_y;
      std::vector< std::vector< std::vector< std::vector<double> > > > disp_z;

      // terms for Fourier expanded config. relaxation
      int max_k;
      std::vector< std::vector<double> > coeffs;

      // variables for Fourier basis strain
      std::vector<double> theta_list;
      std::vector< std::vector<double> > coeffs_x;
      std::vector< std::vector<double> > coeffs_y;
      std::vector< std::vector<double> > coeffs_z;

      Job_params opts;


    public:

      StrainCalc();
      StrainCalc(const StrainCalc& orig);
      ~StrainCalc();

      void loadConfigFile(std::string config_filename);
      void loadFourierConfigFile(std::string config_filename);
      void loadFourierConfigFile_interp(std::string thetas_filename, std::string x_filename, std::string y_filename, std::string z_filename);
	    void setOpts(Job_params opts_in);

      std::vector<double> fourierStrainDisp_config(std::vector<double> r, double* b1, double* b2, int sheet);
      std::vector<double> fourierStrainDisp(double* r, double* b1, double* b2);
      std::vector<double> fourierStrainDisp_sc(double* r, double* b1, double* b2, int s);

      std::vector< std::vector<double> > fourierStrainDisp_vectorized(std::vector< std::vector<double> > r, double* b1, double* b2);
      std::vector<double> fourierStrainDisp_old(std::vector<double> config_in, int sheet, int orb);
      std::vector<double> interpStrainDisp(std::vector<double> config_in, int sheet, int orb);

      double interp_4point(double x, double y, double v1, double v2, double v3, double v4);

  	  std::vector<double> supercellDisp(std::vector<double> pos_in, int sheet, int orb);
  	  std::vector< std::vector<double> > supercellStrain(std::vector<double> pos_in, int sheet, int orb);

  	  std::vector<double> realspaceDisp(std::vector<double> pos_in, int sheet, int orb);
  	  std::vector< std::vector<double> > realspaceStrain(std::vector<double> pos_in, int sheet, int orb);

      double gsfeHeight(std::vector<double> config_in);

};

#endif /* STRAIN_H */
