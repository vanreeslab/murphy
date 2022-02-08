#!/bin/bash
# Submission script for Lemaitre3 
#SBATCH --job-name=epsilon
#SBATCH --time=06:00:00
#
##SBATCH --ntasks=64
##SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=2625
##SBATCH --partition=batch,hmem
##SBATCH --mem-per-cpu=15625
#SBATCH --partition=batch,hmem
#
# #SBATCH --mail-user=thomas.gillis@uclouvain.be
# #SBATCH --mail-type=ALL

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy
SCRATCH=$GLOBALSCRATCH
RUN_DIR=murphy_epsilon_${SLURM_JOB_ID}

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

# go there
cd $SCRATCH/$RUN_DIR 

# get the machine file
MACHINEFILE="nodes.$SLURM_JOBID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE
#==============================================================================
COUNT_GS=2
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
        mv murphy ${MURPHY_NAME}
        # run that shit
        mpirun ./${MURPHY_NAME} --eps --ilevel=5 > log_${N}${NT}_$SLURM_JOB_ID.log
        # destroy the compilation infos
        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.nic5_intel make destroy

        #increment the GP
        let "COUNT_GS+=2"
    done
done

#==============================================================================
#end of file
