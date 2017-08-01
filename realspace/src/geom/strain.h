
/*
 * File:   hstruct.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:43 PM
 */

#ifndef STRAIN_H
#define STRAIN_H

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


    public:

      StrainCalc();
      StrainCalc(const StrainCalc& orig);
      ~StrainCalc();

      void loadConfigFile(std::string config_filename);

      std::vector<double> interpStrainDisp(std::vector<double> config_in, int sheet, int orb);

      double interp_4point(double x, double y, double v1, double v2, double v3, double v4);


};

#endif /* STRAIN_H */
