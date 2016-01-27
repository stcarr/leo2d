/* 
 * File:   hstruct.cpp
 * Author: Stephen Carr
 * 
 * Created on January 13, 2016, 3:16 PM
 */

#include "locality.h"
#include <math.h>
#include <mpi.h>

#include <stdio.h>


int main(int argc, char** argv) {
	Locality loc;
	loc.setup();
	loc.initMPI(argc, argv);
	loc.constructGeom();
	loc.finMPI();
}

Locality::Locality() {
    
}

Locality::Locality(const Locality& orig) {

}

Locality::~Locality() {

}

void Locality::setup() {

}

void Locality::initMPI(int argc, char** argv){

	MPI::Init(argc, argv);
	
	root = 0;
	
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();
}

void Locality::constructGeom(){
		
	if (rank == root){
		// construct sheet, hstruct objects
		//
		
		// max_index = hstruct.getMaxIndex();
		max_index = 1000;
	}
	
	
	MPI::COMM_WORLD.Bcast(&max_index, 1, MPI_INT, 0);
	int pairs_i[max_index];
	int pairs_j[max_index];
	
	if (rank == root){
	
		pairs = new int*[max_index];
		
		for(int k = 0; k < max_index; ++k){
			pairs[k] = new int[2];
			pairs[k][0] = 13;
			pairs[k][1] = 25;
		}
		
		for(int k = 0; k < max_index; ++k){
			pairs_i[k] = pairs[k][0];
			pairs_j[k] = pairs[k][1];
		}
		
	}
	
	
	MPI::COMM_WORLD.Bcast(pairs_i, max_index, MPI_INT, 0);
	MPI::COMM_WORLD.Bcast(pairs_j, max_index, MPI_INT, 0);
	
	printf("rank %d recieved pairs variables! \n", rank);

}

void Locality::constructMatrix(){ 

}

void Locality::solveMatrix(){ 

}

void Locality::getRaw(){ 

}

void Locality::getProcessed(){ 

}

void Locality::plot(){ 

}

void Locality::save(){ 

}

void Locality::finMPI(){ 

	MPI::Finalize();

}

