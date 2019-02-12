/*
 * File:   ballistic.cpp
 * Author: Stephen
 *
 * Created on May 09, 2018, 01:23 PM
 */

#include "ballistic.h"
#include "tools/numbers.h"
#include "matrix/dmatrix.h"

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdexcept>


using namespace numbers;

void Ballistic::runBallisticTransport(Job_params& job, int* index_to_grid, double* i2pos, SpMatrix& H){

  int N = H.getNumRows();

  printf("generating Psi0 \n");
  generatePsi0(job, index_to_grid, i2pos, N);
  printf("Starting propagation \n");
  propagate(job, H);
  printf(" Done with ballistic transport \n");

}

// Matrix exponential function calculation
// Implementation by Po Liu using Intel's Math Kernel Library 11.0.4 in C++
// Algorithm based Arsigny's PhD thesis (http://www.ece.ubc.ca/~purang/content/abolmaesumi/research/interventional-ultrasound.)
// Modified for LEO2D

// "Scaling and Squaring" Method
DMatrix Ballistic::getExpH(Job_params& job, SpMatrix& H){

  // number of steps to approximate the small matrix's exponential
	int accuracy = 8;

	// Scaling factor
	int N = 5;

  // matrix size
  int s = H.getNumRows();

  double dt = job.getDouble("ballistic_time_step");

  // Scaling part
  printf("Starting scaling \n");
	//H_small = -1i*delta_t*H/(2^N);
  DMatrix H_dt;
  H.denseConvert(H_dt);
  std::complex<double> fact = std::complex<double>(0.0, -dt/pow(2.0,N));
  H_dt.scalarMultiply(fact);

  DMatrix H_exp;
  H_exp.setupAsIdentity(s,s, 1);

  DMatrix H_power;
  DMatrix H_power_prev;
  H_power_prev.setupAsIdentity(s,s, 1);


	// Exponentiate approximately
  printf("Starting H_exp \n");
  printf("step (out of %d): \n", accuracy);

	double factorial_i = 1.0;
	for(int i = 1; i < accuracy; i++) {
		factorial_i = factorial_i * i;
    printf("step %d \n",i);

    // Sets H_power = (H_dt)^i/(i!)
    H_dt.matrixMultiply(H_power, H_power_prev, std::complex<double>(1.0/factorial_i, 0.0), std::complex<double>(0.0, 0.0));
    H_exp.matrixAdd(H_power);
    H_power_prev = H_power;

	}

	// Squaring part
  printf("Starting H_sq \n");

  DMatrix H_sq;
  DMatrix H_sq_prev;
  H_exp.matrixMultiply(H_sq_prev, H_exp, std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0));

	for(int i = 1; i < N; i++) {
		// H_sq = H_sq_prev*H_sq_prev;
    H_sq_prev.matrixMultiply(H_sq, H_sq_prev, std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0));
    H_sq_prev = H_sq;
	}

  // Save results

  /*
  printf("Checking Unitary... \n");
  DMatrix unitary_check;
  H_sq.matrixMultiply(unitary_check, H_sq, std::complex<double>(1.0, 0.0),  std::complex<double>(0.0, 0.0), 'N', 'C');
  unitary_check.debugPrint();
  */

  printf("Returning the exp. H matrix \n");
  return H_sq;

}

// Create an initial wavepacket
void Ballistic::generatePsi0(Job_params& job, int* index_to_grid, double* i2pos, int local_max_index){

  	double sigma = job.getDouble("ballistic_sigma");

    // here we avoid the convenient std::complex because job_params does not support complex<doubles> (yet...)
    std::vector<double> psi0_r;
    std::vector<double> psi0_c;
    psi0_r.resize(local_max_index);
    psi0_c.resize(local_max_index);
    double psi0_mag = 0.0;

    for (int k_i = 0; k_i < local_max_index; ++k_i){
      // x,y,z position
      double x = i2pos[k_i*3 + 0];
      double y = i2pos[k_i*3 + 1];
      double z = i2pos[k_i*3 + 2];
      // [i,j] unitcell grid position, orbit index o, and sheet number s
      int i = index_to_grid[k_i*4 + 0];
      int j = index_to_grid[k_i*4 + 1];
      int o = index_to_grid[k_i*4 + 2];
      int s = index_to_grid[k_i*4 + 3];

      double psi0_r_here = pow(E,-(x*x + y*y)/(2*sigma*sigma));
      double psi0_c_here = 0.0;

      psi0_mag = psi0_mag + psi0_r_here*psi0_r_here + psi0_c_here*psi0_c_here;
      psi0_r[k_i] = psi0_r_here;
      psi0_c[k_i] = psi0_c_here;
    }

    double norm = sqrt(psi0_mag);

    for (int k_i = 0; k_i < local_max_index; ++k_i){
      psi0_r[k_i] = psi0_r[k_i]/norm;
      psi0_c[k_i] = psi0_c[k_i]/norm;

    }

    job.setParam("psi0_r",psi0_r);
    job.setParam("psi0_c",psi0_c);

}

// Propogate the wavepacket for the given number of iterations
void Ballistic::propagate(Job_params& job, SpMatrix& H){

  printf("Starting exponential generation \n");
  DMatrix H_dt_exp = getExpH(job, H);
  printf("Done with exponential \n");

  int T = job.getInt("ballistic_max_steps"); // max number of time-steps
  int N = H.getNumRows(); // number of orbitals

  std::complex<double>* psi = new std::complex<double>[N];
  std::complex<double>* psi_next = new std::complex<double>[N];

  std::vector<double> psi0_r = job.getDoubleVec("psi0_r");
  std::vector<double> psi0_c = job.getDoubleVec("psi0_c");

  for (int k = 0; k < N; ++k){
    psi[k] = std::complex<double>(psi0_r[k],psi0_c[k]);
    psi_next[k] = std::complex<double>(0.0);
  }

  std::complex<double> a = std::complex<double>(1.0,0.0);
  std::complex<double> b = std::complex<double>(0.0,0.0);

  std::vector< std::vector<double> > psi_r;
  std::vector< std::vector<double> > psi_c;
  psi_r.resize(T);
  psi_c.resize(T);

  for (int i = 0; i < T; ++i){
    psi_r[i].resize(N);
    psi_c[i].resize(N);

    H_dt_exp.vectorMultiply(psi, psi_next, a, b);

    double sq_norm = 0.0;
    for (int k1 = 0; k1 < N; ++k1){
      sq_norm += psi_next[k1].real()*psi_next[k1].real() + psi_next[k1].imag()*psi_next[k1].imag();
    }
    double norm = sqrt(sq_norm);
    for (int k = 0; k < N; ++k){
      psi_r[i][k] = psi_next[k].real()/norm;
      psi_c[i][k] = psi_next[k].imag()/norm;
      psi[k] = psi_next[k]/norm;
    }
  }

  job.setParam("psi_r",psi_r);
  job.setParam("psi_c",psi_c);
  printf("done with Ballistic Transport calculations!\n");

}
