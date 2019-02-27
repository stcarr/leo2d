/*
 * File:   mat_comm.cpp
 * Author:  Stephen Carr
 *
 * Created on Feb 27, 2019, 9:07 AM
 */

#include <mpi.h>
#include "mpi/mat_comm.h"


void MPI_Bcast_root_loadedMat(int root, LoadedMat &dat){

  MPI::COMM_WORLD.Bcast(&(dat.num_intra_data), 1, MPI_INT, root);
  MPI::COMM_WORLD.Bcast(&(dat.num_inter_data), 1, MPI_INT, root);

  for (int i = 0; i < dat.num_intra_data; ++i){
      MPI_Bcast_root_loadedIntra(root, dat.intra_data[i])
  }
  for (int i = 0; i < dat.num_inter_data; ++i){
      MPI_Bcast_root_loadedInter(root, dat.inter_data[i])
  }

}

void MPI_Bcast_loadedMat(int root, LoadedMat &dat){


  MPI::COMM_WORLD.Bcast(&(dat.num_intra_data), 1, MPI_INT, root);
  MPI::COMM_WORLD.Bcast(&(dat.num_inter_data), 1, MPI_INT, root);

  dat.intra_data.resize(dat.num_intra_data);
  dat.inter_data.resize(dat.num_inter_data);


  for (int i = 0; i < dat.num_intra_data; ++i){
      MPI_Bcast_loadedIntra(root, dat.intra_data[i])
  }
  for (int i = 0; i < dat.num_inter_data; ++i){
      MPI_Bcast_loadedInter(root, dat.inter_data[i])
  }

}

// MPI Broadcast routines so we don't have to do multiple file reads
void MPI_Bcast_root_loadedIntra(int root, LoadedIntraData &dat){

        MPI::COMM_WORLD.Bcast(&(dat.num_orbs), 1, MPI_INT, root);
        int max_R = dat.intralayer_terms.size();
        MPI::COMM_WORLD.Bcast(&(dat.max_R), 1, MPI_INT, root);

        MPI::Status status;

        // Send name
        int len = (int) dat.name.size();
        MPI::COMM_WORLD.Bcast(&len, 1, MPI_INT, root);
        MPI::COMM_WORLD.Bcast(&(dat.name[0]), len, MPI_CHAR, root);

        // Send lattice
        double temp_lat[4];
        int counter = 0;
        for (int a = 0; a < 2; ++a){
                for (int d = 0; d < 2; ++d){
                        temp_lat[counter] = dat.lattice[a][d];
                        counter++;
                }
        }
        MPI::COMM_WORLD.Bcast(temp_lat, 4, MPI::DOUBLE, root);

        // Send orbital positions
        double temp_pos[3*dat.num_orbs];
        counter = 0;
        for (int o = 0; o < dat.num_orbs; ++o){
                for (int d = 0; d < 3; ++d){
                        temp_pos[counter] = dat.orb_pos[o][d];
                        counter++;
                }
        }
        MPI::COMM_WORLD.Bcast(temp_pos, 3*dat.num_orbs, MPI::DOUBLE, root);

        // Send intralayer terms
        int dim = max_R*max_R*dat.num_orbs*dat.num_orbs;
        double temp_val[dim];

        counter = 0;
        for (int i = 0; i < max_R; ++i){
                for (int j = 0; j < max_R; ++j){
                        for (int k = 0; k < dat.num_orbs; ++k){
                                for (int l = 0; l < dat.num_orbs; ++l){
                                        temp_val[counter] = dat.intralayer_terms[i][j][k][l];
                                        ++counter;
                                }
                        }
                }
        }

        MPI::COMM_WORLD.Bcast(temp_val, dim, MPI::DOUBLE, root);
}

void MPI_Bcast(int root, LoadedIntraData &dat){
        MPI::COMM_WORLD.Bcast(&(dat.num_orbs), 1, MPI_INT, root);
        int max_R;
        MPI::COMM_WORLD.Bcast(&(dat.max_R), 1, MPI_INT, root);

        MPI::Status status;

        // Recieve name
        int len;
        MPI::COMM_WORLD.Bcast(&len, 1, MPI_INT, root);
        char temp_name[len];
        MPI::COMM_WORLD.Bcast(temp_name, len, MPI_CHAR, root);
        dat.name.assign(temp_name,len);

        // Recieve lattice
        double temp_lat[4];
        MPI::COMM_WORLD.Bcast(temp_lat, 4, MPI::DOUBLE, root);


        int counter = 0;
        for (int a = 0; a < 2; ++a){
                for (int d = 0; d < 2; ++d){
                        dat.lattice[a][d] = temp_lat[counter];
                        counter++;
                }
        }

        // Recieve orbital positions
        double temp_pos[3*dat.num_orbs];

        MPI::COMM_WORLD.Bcast(temp_pos, 3*dat.num_orbs, MPI::DOUBLE, root);

        counter = 0;
        dat.orb_pos.resize(dat.num_orbs);
        for (int o = 0; o < dat.num_orbs; ++o){
                for (int d = 0; d < 3; ++d){
                        dat.orb_pos[o][d] = temp_pos[counter];
                        counter++;
                }
        }

        // Recieve intralayer terms
        int dim = max_R*max_R*dat.num_orbs*dat.num_orbs;
        double temp_val[dim];

        counter = 0;
        MPI::COMM_WORLD.Bcast(temp_val, dim, MPI::DOUBLE, root);

        dat.intralayer_terms.resize(max_R);
        for (int i = 0; i < max_R; ++i){
                dat.intralayer_terms[i].resize(max_R);
                for (int j = 0; j < max_R; ++j){
                        dat.intralayer_terms[i][j].resize(dat.num_orbs);
                        for (int k = 0; k < dat.num_orbs; ++k){
                                dat.intralayer_terms[i][j][k].resize(dat.num_orbs);
                                for (int l = 0; l < dat.num_orbs; ++l){
                                        dat.intralayer_terms[i][j][k][l] = temp_val[counter];
                                        ++counter;
                                }
                        }
                }
        }
}

void MPI_Bcast_root_loadedInter(int root, LoadedInterData &dat){

}

void MPI_Bcast_loadedInter(int root, LoadedInterData &dat){

}
