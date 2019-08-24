/*
 * File:   mat_comm.h
 * Author:  Stephen Carr
 *
 * Created on Feb 27, 2019, 9:07 AM
 */

#ifndef mat_comm_h
#define mat_comm_h

#include "materials/loaded_mat.h"


void MPI_Bcast_root_loadedMat(int root, LoadedMat &dat);
void MPI_Bcast_loadedMat(int root, LoadedMat &dat);

void MPI_Bcast_root_loadedIntra(int root, LoadedIntraData &dat);
void MPI_Bcast_loadedIntra(int root, LoadedIntraData &dat);

void MPI_Bcast_root_loadedInter(int root, LoadedInterData &dat);
void MPI_Bcast_loadedInter(int root, LoadedInterData &dat);

#endif
