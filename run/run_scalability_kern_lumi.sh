#!/bin/bash
#SBATCH --partition=standard,debug
#SBATCH --time=00:20:00
#SBATCH --ntasks-per-node=128
#SBATCH --account=project_465000036
#SBATCH --exclusive

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "max rma calls: $MPICH_RMA_MAX_PENDING"
echo "------------------------"


# run
MPICH_RMA_MAX_PENDING=2048 srun -u ./murphy --weak-scal --profile
#srun -u ./murphy --weak-scal --profile
#srun -n 128 -u ./murphy --weak-scal --profile



