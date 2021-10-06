#!/bin/bash
#SBATCH --qos=debug
#SBATCH --time=0:29:00
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
#SBATCH --account=m3640

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

# run
srun ./murphy --weak-scal --profile



