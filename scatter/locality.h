
/* 
 * File:   locality.h
 * Author: Stephen Carr
 *
 * Created on January 27, 2016, 2:39 PM
 */

#ifndef LOCALITY_H
#define LOCALITY_H
#include "hstruct.h"

class Locality {
    private:
		int size;
		int rank;
		int max_index;
		int root;
		int** pairs;
        
    public:
        Locality();
        Locality(const Locality& orig);
        ~Locality();
		void setup();
		void initMPI(int, char**);
		void constructGeom();
		void constructMatrix();
		void solveMatrix();
		void getRaw();
		void getProcessed();
		void plot();
		void save();
		void finMPI();

};

#endif /* LOCALITY_H */