#!/bin/bash
# Submission script for Engaging 
#SBATCH --job-name=murphy-epsilon
#SBATCH --time=12:00:00
#
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=5000mb
#SBATCH --partition=sched_mit_wvanrees
#

HOME_MURPHY=/home/tgillis/murphy
SCRATCH=/nobackup1/tgillis
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
cp -r $HOME_MURPHY/make_arch/make.engaging $SCRATCH/$RUN_DIR/make_arch


# go there
cd $SCRATCH/$RUN_DIR 

#==============================================================================
## wavelet 2.0
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=2 -DWAVELET_NT=0 -DBLOCK_GS=2"
mv murphy murphy20
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy20 --eps --ilevel=7 > log_20_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.engaging make destroy

#==============================================================================
## wavelet 2.2
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=2 -DWAVELET_NT=2 -DBLOCK_GS=4"
mv murphy murphy22
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy22 --eps --ilevel=7 > log_22_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.engaging make destroy

#==============================================================================
## wavelet 4.0
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy40
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy40 --eps --ilevel=7 > log_40_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.engaging make destroy

#==============================================================================
## wavelet 4.2
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=8"
mv murphy murphy42
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy42 --eps --ilevel=7 > log_42_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.engaging make destroy

#==============================================================================
## wavelet 6.0
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=6 -DWAVELET_NT=0 -DBLOCK_GS=10"
mv murphy murphy60
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy60 --eps --ilevel=7 > log_60_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.engaging make destroy

#==============================================================================
## wavelet 6.2
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=6 -DWAVELET_NT=2 -DBLOCK_GS=12"
mv murphy murphy62
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy62 --eps --ilevel=7 > log_62_$SLURM_JOB_ID.log
ARCH_FILE=make_arch/make.engaging make destroy




#end of file
