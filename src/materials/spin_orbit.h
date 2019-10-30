/*
 * File:   spin_orbit.h
 * Author:  Stephen Carr
 *
 * Created on Oct 29, 2019, 3:41 PM
 */

#ifndef spin_orbit_h
#define spin_orbit_h

#include <array>

#include "matrix/spmatrix.h"
#include "params/job_params.h"


/**
 * This module adds approximate spin-orbit terms based on passed parameters
 */

namespace SpinOrbit {

  // make pairs and complex values for spin-orbit correction
  void generateSOC(std::vector< std::vector<int> >& soc_pairs, std::vector< std::complex<double> >& soc_terms,
     int* index_to_grid, int max_index, int local_max_index, std::vector<int> current_index_reduction,
     Job_params opts);

}
#endif
