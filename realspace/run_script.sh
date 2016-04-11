#!/bin/bash
#
#SBATCH -n 9 # Number of MPI ranks
#SBATCH -c 4 #CPU/Threads per rank
#SBATCH -J hstruct_test
#SBATCH -o hstruct_test_%j.out
#SBATCH -e hstruct_test_%j.err
#SBATCH -t 0-02:00 # Runtime limit
#SBATCH -p kaxiras # Partition to submit to (lab group = -p kaxiras)
#SBATCH --mem-per-cpu=4000 # Memory per cpu in MB (see also --mem)

# Set OMP_NUM_THREADS to the same value as -c
# with a fallback in case it isn't set.
# SLURM_CPUS_PER_TASK is set to the value of -c, but only if -c is explicitly set
if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  omp_threads=$SLURM_CPUS_PER_TASK
else
  omp_threads=1
fi
export OMP_NUM_THREADS=$omp_threads


module load intel/15.0.0-fasrc01
module load intel-mkl/11.0.0.079-fasrc02
module load mvapich2/2.2b-fasrc02

echo $OMP_NUM_THREADS

srun -n $SLURM_NTASKS --mpi=pmi2 ./main hstruct.in
#mpirun -n $SLURM_NTASKS ./main hstruct.in
