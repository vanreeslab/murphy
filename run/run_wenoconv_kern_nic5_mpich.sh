#!/bin/bash
# Submission script for Engaging
#SBATCH --time=12:00:00
#
#SBATCH --ntasks-per-node=64
##SBATCH --mem-per-cpu=15625
#SBATCH --mem-per-cpu=4000
#SBATCH --partition=batch,hmem
#

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy

# create the tmp directory
mkdir -p ${RUN_DIR}

# go there
cd ${RUN_DIR} 
cp ${HOME_MURPHY}/make_arch ${RUN_DIR}/make_arch
cp ${HOME_MURPHY}/src ${RUN_DIR}/src

# get the machine file
MACHINEFILE="nodes.$SLURM_JOBID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE

# 4.0
M=4
NT=0
make destroy
ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mkdir -p w40
cd w40
k
/home/ucl/tfl/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun -f ${MACHINEFILE} ./${MNAME} --conv-weno > log_$SLURM_JOB_ID.log
/home/ucl/tfl/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun -f ${MACHINEFILE} ./${MNAME} --conv-weno --no-adapt > log_$SLURM_JOB_ID.log


echo "------------------------"

