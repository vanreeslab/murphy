#!/bin/bash
# Submission script for Engaging 
#SBATCH --job-name=murphy-abc
#SBATCH --time=12:00:00
#
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=5000mb
#SBATCH --partition=sched_mit_wvanrees
#
# #SBATCH --mail-user=thomas.gillis@uclouvain.be
# #SBATCH --mail-type=ALL

HOME_MURPHY=/home/tgillis/murphy
SCRATCH=/nobackup1/tgillis
RUN_DIR=murphy_abc_${SLURM_JOB_ID}
#RUN_DIR=murphy_weak/murphy_weak_${SLURM_JOB_NUM_NODES}_${SLURM_JOB_ID}

# create the tmp directory
mkdir -p $SCRATCH/$RUN_DIR
mkdir -p $SCRATCH/$RUN_DIR/data
mkdir -p $SCRATCH/$RUN_DIR/prof

# copy what is needed
cp -r $HOME_MURPHY/murphy $SCRATCH/$RUN_DIR

# go there
cd $SCRATCH/$RUN_DIR 


OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy --abc --ilevel=5 --iter-max=100000 --rtol=1e-4 --ctol=1e-5 --level-max=9 --iter-dump=10 > log_$SLURM_JOB_ID.log
#OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy --abc --ilevel=4 --iter-max=100000 --rtol=1e-2 --ctol=1e-5 --level-max=8 --iter-dump=50 > log_$SLURM_JOB_ID.log


