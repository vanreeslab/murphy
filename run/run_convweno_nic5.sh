#!/bin/bash
#SBATCH --job-name=conv-weno
#SBATCH --time=02:00:00
#
##SBATCH --ntasks=64
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2625
#SBATCH --partition=batch,hmem
#

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy
SCRATCH=$GLOBALSCRATCH
RUN_DIR=murphy_convweno_${SLURM_JOB_ID}

# create the tmp directory
mkdir -p $SCRATCH/$RUN_DIR
mkdir -p $SCRATCH/$RUN_DIR/src
mkdir -p $SCRATCH/$RUN_DIR/data
mkdir -p $SCRATCH/$RUN_DIR/prof
mkdir -p $SCRATCH/$RUN_DIR/make_arch

# copy what is needed
cp -r $HOME_MURPHY/Makefile $SCRATCH/$RUN_DIR
cp -r $HOME_MURPHY/src $SCRATCH/$RUN_DIR/src
cp -r $HOME_MURPHY/make_arch/make.nic5_intel $SCRATCH/$RUN_DIR/make_arch

# get the commit id
cd $HOME_MURPHY
GITCOMMIT=$(git describe --always --dirty)
echo "commit = $GITCOMMIT"

# go there
cd $SCRATCH/$RUN_DIR 

#==============================================================================
COUNT_GS=4
for N in {2,4,6}
do
    for NT in {0,2}
    do
        echo "================================================================="
        echo "      WAVELET ${N}.${NT} -> COUNT = ${COUNT_GS}"
        #compilation 
        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.nic5_intel make -j OPTS="-DWAVELET_N=${N} -DWAVELET_NT=${NT} -DBLOCK_GS=${COUNT_GS}"
        # move to the full name
        MURPHY_NAME=murphy_${N}_${NT}
        cp murphy ${MURPHY_NAME}

        for RE in {-1,0.01,0.1,1,10,100}
        do
            echo "Reynolds = ${RE}"
            mpirun ./${MURPHY_NAME} --2lvl-weno --level-max=4 --reynolds=${RE} > log_w3_re${RE}_${N}${NT}_$SLURM_JOB_ID.log
            mpirun ./${MURPHY_NAME} --2lvl-weno --level-max=4 --fix-weno --reynolds=${RE} > log_c3_re${RE}_${N}${NT}_$SLURM_JOB_ID.log
        done

        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.nic5_intel make clean

        #increment the GP
        if ! { [ $N -eq 2 ] && [ $NT -eq 0 ]; }; then
            let "COUNT_GS+=2"
        fi
    done
done

#==============================================================================
#end of file
