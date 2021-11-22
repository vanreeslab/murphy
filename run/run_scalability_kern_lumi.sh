#!/bin/bash
#SBATCH --partition=debug
#SBATCH --time=0:30:00
#SBATCH --ntasks-per-node=128
#SBATCH --account=project_465000036

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

# run
srun -n 128 -u ./murphy --weak-scal --profile



