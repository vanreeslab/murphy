#!/bin/bash
# Submission script for Zenobe
#PBS -N murphy
#PBS -q large
#PBS -W group_list=examples
##PBS -l walltime=01:00:00
##PBS -l select=8:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1
#PBS -l walltime=10:00:00
#PBS -l select=30:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1

HOME_MURPHY=/home/acad/ucl-tfl/tgillis/murphy
SCRATCH=/SCRATCH/acad/examples/tgillis/
RUN_DIR=murphy_epsilon_${PBS_JOBID}

# create the tmp directory
mkdir -p $SCRATCH/$RUN_DIR
mkdir -p $SCRATCH/$RUN_DIR/src
mkdir -p $SCRATCH/$RUN_DIR/data
mkdir -p $SCRATCH/$RUN_DIR/prof
mkdir -p $SCRATCH/$RUN_DIR/make_arch

# copy what is needed
cp -r $HOME_MURPHY/Makefile $SCRATCH/$RUN_DIR
cp -r $HOME_MURPHY/src $SCRATCH/$RUN_DIR/src
cp -r $HOME_MURPHY/make_arch/make.zenobe $SCRATCH/$RUN_DIR/make_arch

# get the commit id
cd $HOME_MURPHY
GITCOMMIT=$(git describe --always --dirty)

# go there
cd $SCRATCH/$RUN_DIR

# get the machine file
#==============================================================================
COUNT_GS=2
for N in {2,4,6}
do
    for NT in {0,2}
    do
        echo "================================================================="
        echo "      WAVELET ${N}.${NT} -> COUNT = ${COUNT_GS}"
        #compilation 
        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.zenobe make clean
        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.zenobe make -j OPTS="-DWAVELET_N=${N} -DWAVELET_NT=${NT} -DBLOCK_GS=${COUNT_GS}"

        # move to the full name
        MURPHY_NAME=murphy_${N}_${NT}
        mv murphy ${MURPHY_NAME}

        # run that shit
        mpirun --version
        mpirun --mca osc pt2pt ./${MURPHY_NAME} --eps --ilevel=6 > log_${N}${NT}_${PBS_JOBID}.log

        # destroy the compilation infos
        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.zenobe make clean

        #increment the GP
        let "COUNT_GS+=2"
    done
done

#==============================================================================
#end of file
