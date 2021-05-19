#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=epsilon
#SBATCH --time=08:00:00
#
##SBATCH --ntasks=64
##SBATCH --ntasks-per-node=64
##SBATCH --mem-per-cpu=2625
##SBATCH --partition=batch,hmem
#SBATCH --mem-per-cpu=15625
#SBATCH --partition=batch,hmem
#
# #SBATCH --mail-user=thomas.gillis@uclouvain.be
# #SBATCH --mail-type=ALL

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy
SCRATCH=${GLOBALSCRATCH}
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
cp -r $HOME_MURPHY/make_arch/make.nic5_mpich $SCRATCH/$RUN_DIR/make_arch


# go there
cd $SCRATCH/$RUN_DIR 

MPIRUN=${HOME}/lib-mpich/mpich-3.4.1/bin/mpirun -env UCX_NET_DEVICES=mlx5_0:1

#==============================================================================
## wavelet 4.0
ARCH_FILE=make_arch/make.nic5_mpich make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy40
${MPIRUN} ./murphy40 --eps --ilevel=5 > log_40_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.nic5_mpich make destroy

#==============================================================================
## wavelet 4.2
ARCH_FILE=make_arch/make.nic5_mpich make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=8"
mv murphy murphy42
${MPIRUN} ./murphy42 --eps --ilevel=5 > log_42_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.nic5_mpich make destroy


#end of file
