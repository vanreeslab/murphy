#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=weak-scalability
#SBATCH --time=01:00:00
#
##SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2625
# ##SBATCH --partition=batch,debug 
#
# #SBATCH --mail-user=thomas.gillis@uclouvain.be
# #SBATCH --mail-type=ALL

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy
SCRATCH=$GLOBALSCRATCH
RUN_DIR=murphy_weak_${SLURM_JOB_NUM_NODES}_${SLURM_JOB_ID}

# create the tmp directory
mkdir -p $SCRATCH/$RUN_DIR
mkdir -p $SCRATCH/$RUN_DIR/data
mkdir -p $SCRATCH/$RUN_DIR/prof

# copy what is needed
cp -r $HOME_MURPHY/murphy $SCRATCH/$RUN_DIR

# go there
cd $SCRATCH/$RUN_DIR 

# run that shit
#mpirun ./murphy --abc --profile --iter-max=500 --iter-dump=10 --dump-detail=1 --rtol=1e-2 --ctol=1e-4 --vr-sigma=0.025 --dom=1,1,$SLURM_JOB_NUM_NODES > log_$SLURM_JOB_ID.log
mpirun ./murphy --abc --profile --iter-max=100 --iter-dump=1000 --rtol=1e-2 --ctol=1e-4 --vr-sigma=0.01 --dom=1,1,$SLURM_JOB_NUM_NODES > log_$SLURM_JOB_ID.log
#mpirun -d valgrind --error-limit=no --suppressions=${EBROOTOPENMPI}/share/openmpi/openmpi-valgrind.supp ./murphy --ns --profile --iter-max=10 > $SLURM_JOB_ID.log
