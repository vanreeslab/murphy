#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=3925
#SBATCH --partition=batch,hmem

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

# run
#srun ./murphy --weak-scal --profile
mpirun -n ${SLURM_NPROCS} ./murphy --weak-scal --profile



