#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=weak-scalability
#SBATCH --time=08:00:00
#
##SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2625
##SBATCH --partition=batch,hmem
#SBATCH --mem-per-cpu=15625
#SBATCH --partition=hmem
#
# #SBATCH --mail-user=thomas.gillis@uclouvain.be
# #SBATCH --mail-type=ALL

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy
SCRATCH=$GLOBALSCRATCH
RUN_DIR=murphy_epsilon_${SLURM_JOB_ID}
#RUN_DIR=murphy_weak/murphy_weak_${SLURM_JOB_NUM_NODES}_${SLURM_JOB_ID}

# create the tmp directory
mkdir -p $SCRATCH/$RUN_DIR
mkdir -p $SCRATCH/$RUN_DIR/src
mkdir -p $SCRATCH/$RUN_DIR/data
mkdir -p $SCRATCH/$RUN_DIR/prof
mkdir -p $SCRATCH/$RUN_DIR/make_arch

# copy what is needed
cp -r $HOME_MURPHY/Makefile $SCRATCH/$RUN_DIR
cp -r $HOME_MURPHY/src $SCRATCH/$RUN_DIR/src
cp -r $HOME_MURPHY/make_arch/make.nic5 $SCRATCH/$RUN_DIR/make_arch


# go there
cd $SCRATCH/$RUN_DIR 

#==============================================================================
## wavelet 4.2
ARCH_FILE=make_arch/make.nic5 make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=4"
mv murphy murphy42
#mpirun ./murphy42 --sadv --profile --iter-max=5000 --iter-dump=10000 --rtol=1e-2 --ctol=1e-4 --grid-on-sol > log_42_$SLURM_JOB_ID.log
#mpirun ./murphy42 --sadv --profile --iter-max=5000 --iter-dump=10000 --rtol=1e-2 --ctol=1e-4 --grid-on-sol --no-weno > log_42_$SLURM_JOB_ID.log
mpirun ./murphy42 --eps --ilevel=7 > log_42_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.nic5 make destroy

#==============================================================================
## wavelet 4.4
ARCH_FILE=make_arch/make.nic5 make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=4 -DBLOCK_GS=6"
mv murphy murphy44
#mpirun ./murphy44 --sadv --profile --iter-max=5000 --iter-dump=10000 --rtol=1e-2 --ctol=1e-4 --grid-on-sol > log_44_$SLURM_JOB_ID.log
#mpirun ./murphy44 --sadv --profile --iter-max=5000 --iter-dump=10000 --rtol=1e-2 --ctol=1e-4 --grid-on-sol --no-weno > log_44_$SLURM_JOB_ID.log
mpirun ./murphy44 --eps --ilevel=7 > log_44_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.nic5 make destroy

#==============================================================================
## wavelet 4.0
ARCH_FILE=make_arch/make.nic5 make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=4"
mv murphy murphy40
#mpirun ./murphy40 --sadv --profile --iter-max=5000 --iter-dump=10000 --rtol=1e-2 --ctol=1e-4 --grid-on-sol > log_40_$SLURM_JOB_ID.log
#mpirun ./murphy40 --sadv --profile --iter-max=5000 --iter-dump=10000 --rtol=1e-2 --ctol=1e-4 --grid-on-sol --no-weno > log_40_$SLURM_JOB_ID.log
mpirun ./murphy40 --eps --ilevel=7 > log_40_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.nic5 make destroy




