/*
 * File:   ballistic.h
 * Author: Stephen
 *
 * Created on May 09, 2018, 01:23 PM
 */
#ifndef BALLISTIC_H
#define BALLISTIC_H


#include "params/job_params.h"
#include "matrix/dmatrix.h"
#include "matrix/spmatrix.h"

/**
 * A namespace for all the processing methods for the ballistic class
 */
namespace Ballistic{

  void runBallisticTransport(Job_params& job, int* index_to_grid, double* i2pos, SpMatrix& H);

	// Create an initial wavepacket
  void generatePsi0(Job_params& job, int* index_to_grid, double* i2pos, int local_max_index);

  // Get an exponetial of H for time propogation
  DMatrix getExpH(Job_params& job, SpMatrix& H);

  // Propagate the wavepacket for the given number of iterations
  void propagate(Job_params& job, SpMatrix& H);
}

#endif /* BALLISTIC_H */
