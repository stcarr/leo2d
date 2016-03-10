#!/bin/bash
#
#SBATCH -n 2 # Number of MPI ranks
#SBATCH --cpus-per-task=8 #CPU/Threads per rank
#SBATCH -J hstruct_test
#SBATCH -o hstruct_test_%j.out
#SBATCH -e hstruct_test_%j.err
#SBATCH -t 0-00:30 # Runtime limit
#SBATCH -p kaxiras # Partition to submit to (lab group = -p kaxiras)
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem)

module load intel/15.0.0-fasrc01
module load intel-mkl/11.0.0.079-fasrc02
module load mvapich2/2.0-fasrc03

mpiexec -n 2 ./main hstruct.in
