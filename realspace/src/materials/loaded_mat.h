/*
 * File:   loaded_mat.h
 * Author:  Stephen Carr
 *
 * Created on May 21, 2018, 5:46 PM
 */

#ifndef loaded_mat_h
#define loaded_mat_h

#include <array>
#include <vector>
#include <mpi.h>

#include "tools/numbers.h"


/**
 * A namespace to group constants specific to the Graphene material.
 */
class LoadedMat{

        public:
          std::string name;
          std::array<std::array<double, 2>, 2> lattice;

          int num_orbs;
          std::vector< std::array<double, 3> > orb_pos;

          // intralyer terms, intralyer_terms[R_i][R_j][o1][o2]
          std::vector< std::vector< std::vector< std::vector<double> > > > intralayer_terms;


          // MPI Broadcast routines so we don't have to do multiple file reads
          void MPI_Bcast_root(int root){

                  MPI::COMM_WORLD.Bcast(&num_orbs, 1, MPI_INT, root);
                  int max_R = intralayer_terms.size();
                  MPI::COMM_WORLD.Bcast(&max_R, 1, MPI_INT, root);

                  MPI::Status status;

                  // Send name
                  int len = (int) name.size();
                  MPI::COMM_WORLD.Bcast(&len, 1, MPI_INT, root);
                  MPI::COMM_WORLD.Bcast(&name[0], len, MPI_CHAR, root);

                  // Send lattice
                  double temp_lat[4];
                  int counter = 0;
                  for (int a = 0; a < 2; ++a){
                          for (int d = 0; d < 2; ++d){
                                  temp_lat[counter] = lattice[a][d];
                                  counter++;
                          }
                  }
                  MPI::COMM_WORLD.Bcast(temp_lat, 4, MPI::DOUBLE, root);

                  // Send orbital positions
                  double temp_pos[3*num_orbs];
                  counter = 0;
                  for (int o = 0; o < num_orbs; ++o){
                          for (int d = 0; d < 3; ++d){
                                  temp_pos[counter] = orb_pos[o][d];
                                  counter++;
                          }
                  }
                  MPI::COMM_WORLD.Bcast(temp_pos, 3*num_orbs, MPI::DOUBLE, root);

                  // Send intralayer terms
                  int dim = max_R*max_R*num_orbs*num_orbs;
                  double temp_val[dim];

                  counter = 0;
                  for (int i = 0; i < max_R; ++i){
                          for (int j = 0; j < max_R; ++j){
                                  for (int k = 0; k < num_orbs; ++k){
                                          for (int l = 0; l < num_orbs; ++l){
                                                  temp_val[counter] = intralayer_terms[i][j][k][l];
                                                  ++counter;
                                          }
                                  }
                          }
                  }

                  MPI::COMM_WORLD.Bcast(temp_val, dim, MPI::DOUBLE, root);
          }

          void MPI_Bcast(int root){

                  MPI::COMM_WORLD.Bcast(&num_orbs, 1, MPI_INT, root);
                  int max_R;
                  MPI::COMM_WORLD.Bcast(&max_R, 1, MPI_INT, root);

                  MPI::Status status;

                  // Recieve name
                  int len;
                  MPI::COMM_WORLD.Bcast(&len, 1, MPI_INT, root);
                  char temp_name[len];
                  MPI::COMM_WORLD.Bcast(temp_name, len, MPI_CHAR, root);
                  name.assign(temp_name,len);

                  // Recieve lattice
                  double temp_lat[4];
                  MPI::COMM_WORLD.Bcast(temp_lat, 4, MPI::DOUBLE, root);


                  int counter = 0;
                  for (int a = 0; a < 2; ++a){
                          for (int d = 0; d < 2; ++d){
                                  lattice[a][d] = temp_lat[counter];
                                  counter++;
                          }
                  }

                  // Recieve orbital positions
                  double temp_pos[3*num_orbs];

                  MPI::COMM_WORLD.Bcast(temp_pos, 3*num_orbs, MPI::DOUBLE, root);

                  counter = 0;
                  orb_pos.resize(num_orbs);
                  for (int o = 0; o < num_orbs; ++o){
                          for (int d = 0; d < 3; ++d){
                                  orb_pos[o][d] = temp_pos[counter];
                                  counter++;
                          }
                  }

                  // Recieve intralayer terms
                  int dim = max_R*max_R*num_orbs*num_orbs;
                  double temp_val[dim];

                  counter = 0;
                  MPI::COMM_WORLD.Bcast(temp_val, dim, MPI::DOUBLE, root);

                  intralayer_terms.resize(max_R);
                  for (int i = 0; i < max_R; ++i){
                          intralayer_terms[i].resize(max_R);
                          for (int j = 0; j < max_R; ++j){
                                  intralayer_terms[i][j].resize(num_orbs);
                                  for (int k = 0; k < num_orbs; ++k){
                                          intralayer_terms[i][j][k].resize(num_orbs);
                                          for (int l = 0; l < num_orbs; ++l){
                                                  intralayer_terms[i][j][k][l] = temp_val[counter];
                                                  ++counter;
                                          }
                                  }
                          }
                  }
          }

    };
#endif
