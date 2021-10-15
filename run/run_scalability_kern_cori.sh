#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=0:30:00
#SBATCH --ntasks-per-node=32
#SBATCH --constraint=haswell
#SBATCH --account=m3640

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

# run
echo "srun -n ${SLURM_NPROCS} -u ./murphy --weak-scal --profile"
srun -n ${SLURM_NPROCS} -u ./murphy --weak-scal --profile



