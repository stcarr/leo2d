/*
 * File:   spin_orbit.cpp
 * Author:  Stephen Carr
 *
 * Created on Oct 29, 2019, 3:41 PM
 */

#include <string>
#include <fstream>
#include <sstream>

#include "materials/spin_orbit.h"

using namespace std;

// get spin-orbit coupling pairs and terms
void SpinOrbit::generateSOC(vector< vector<int> >& soc_pairs, vector< complex<double> >& soc_terms,
   int* index_to_grid, int max_index, int local_max_index, vector<int> current_index_reduction,
   Job_params opts){

  double spin_val = opts.getDouble("spin_val");
  double soc_lambda_val = opts.getDouble("soc_lambda_val");

  for (int k0 = 0; k0 < max_index; ++k0){

    //int k_here = k0 - current_index_reduction[k0]; // cover case of vacancies...?
    int k_here = k0; // cover case of vacancies...?
    int grid_x_h = index_to_grid[4*k0 + 0]; // x grid of unit cell
    int grid_y_h = index_to_grid[4*k0 + 1]; // y grid of unit cell
    int o_h = index_to_grid[4*k0 + 2]; // tight-binding orbital of k0_here
    int s_h = index_to_grid[4*k0 + 3]; // sheet for k0_here

    int save_term = 0;

    int k_now = -1;
    int tar_o = -1;
    double lz_here = 0.0;

    //printf("k_here = %d, o_h = %d \n",k_here, o_h);

    if (o_h == 1){ // d_xz orbital

      k_now = k_here+1;
      tar_o = 2; // d_yz orbital
      lz_here = 1.0;
      save_term = 1;

    } else if (o_h == 2){ // d_yz orbital

          k_now = k_here-1;
          tar_o = 1; // d_xz orbital
          lz_here = -1.0;
          save_term = 1;

    } else if (o_h == 3) { // d_xy orbital

      k_now = k_here+1;
      tar_o = 4; // d_{x^2-y^2} orbital
      lz_here = 2.0;
      save_term = 1;

    } else if (o_h == 4) { // d_{x^2-y^2} orbital

      k_now = k_here-1;
      tar_o = 3; // d_xy orbital
      lz_here = -2.0;
      save_term = 1;

    }

    if (save_term == 1){

      int grid_x_now = index_to_grid[4*k_now + 0]; // x grid of unit cell
      int grid_y_now = index_to_grid[4*k_now + 1]; // y grid of unit cell
      int o_now = index_to_grid[4*k_now + 2]; // tight-binding orbital of k0_here
      int s_now = index_to_grid[4*k_now + 3]; // sheet for k0_here

      //printf("k_now = %d, grid_x_now = %d, grid_y_now = %d, o_now = %d, s_now = %d \n",k_now,grid_x_now, grid_y_now,o_now,s_now);

      if (grid_x_h != grid_x_now || grid_y_h != grid_y_now || s_h != s_now || o_now != tar_o){
        throw std::runtime_error("SpinOrbit::generateSOC() cannot find the correct onsite index! \n");
      }

      complex<double> term_here = complex<double>(0.0,1.0)*soc_lambda_val*spin_val*lz_here;

      //printf("soc_lambda_val = %lf, spin_val = %lf, lz_here = %lf \n", soc_lambda_val, spin_val, lz_here);
      vector<int> pair_here;
      pair_here.push_back(k_here); // (k_here - current_index_reduction[k_here]);
      pair_here.push_back(k_now); // (k_now - current_index_reduction[k_now]);
      soc_pairs.push_back(pair_here);
      soc_terms.push_back(term_here);
      //printf("saved SOC term [%d,%d]: %lf + %lfi \n",k_here, k_now, term_here.real(), term_here.imag());
    }


  }

}
