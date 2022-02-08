#!/bin/bash
#SBATCH --job-name=conv-weno
#SBATCH --time=02:00:00
#
##SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
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

# go there
cd $SCRATCH/$RUN_DIR 

# get the machine file
MACHINEFILE="nodes.$SLURM_JOBID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE
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
        # weno 3
        cp murphy ${MURPHY_NAME}
        mpirun ./${MURPHY_NAME} --conv-weno --ilevel=4 --level-max=9 > log_${N}${NT}_$SLURM_JOB_ID.log
        # weno 5
        ## conservative 3
        #cp murphy ${MURPHY_NAME}_c3
        #mpirun --mca pml ucx --mca osc ucx --mca btl ^uct ./${MURPHY_NAME}_c3 --conv-weno --ilevel=4 --level-max=7 --weno=3 --fix-weno > log_${N}${NT}_c3_$SLURM_JOB_ID.log
        ## conservative 5
        #cp murphy ${MURPHY_NAME}_c5
        #mpirun --mca pml ucx --mca osc ucx --mca btl ^uct ./${MURPHY_NAME}_c5 --conv-weno --ilevel=4 --level-max=7 --weno=5 --fix-weno > log_${N}${NT}_c5_$SLURM_JOB_ID.log
        # destroy the compilation infos
        GIT_COMMIT=${GITCOMMIT} ARCH_FILE=make_arch/make.nic5_intel make destroy

        #increment the GP
        if ! { [ $N -eq 2 ] && [ $NT -eq 0 ]; }; then
            let "COUNT_GS+=2"
        fi
    done
done

#==============================================================================
#end of file
