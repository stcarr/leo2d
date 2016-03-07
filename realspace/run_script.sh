#!/bin/bash
#
#SBATCH -n 4 # Number of cores
#SBATCH --ntasks-per-node=64 # run 64 cores/node (2 nodes)
#SBATCH -J hstruct_test
#SBATCH -o hstruct_test_%j.out
#SBATCH -e hstruct_test_%j.err
#SBATCH -t 0-00:30 # Runtime limit
#SBATCH -p general # Partition to submit to (lab group = -p kaxiras)
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem)

module load gcc/4.8.2-fasrc01 acml/5.3.1-fasrc01
module load mvapich2/2.0-fasrc03
 
mpiexec -n 2 ./main
