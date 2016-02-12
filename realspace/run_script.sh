#!/bin/bash
#
#SBATCH -n 4 # Number of cores
#SBATCH --ntasks-per-node=64 # run 64 cores/node (2 nodes)
#SBATCH -J hstruct_test
#SBATCH -o hstruct_test_%j.out
#SBATCH -e hstruct_test_%j.err
#SBATCH -t 0-00:05 # Runtime limit
#SBATCH -p kaxiras # Partition to submit to
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem)

module load intel/15.0.0-fasrc01 mvapich2/2.0-fasrc03 slepc/3.5.4-fasrc02
module unload Anaconda/1.9.2-fasrc01
 
mpirun -n 4 ./main
