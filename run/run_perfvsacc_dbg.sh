#!/bin/bash
# Submission script for Engaging
#SBATCH --time=12:00:00
##SBATCH --nodelist=node[1019-1022,1029-1030]
#
#SBATCH --ntasks-per-node=36
#SBATCH --mem-per-cpu=5000
#SBATCH --partition=sched_mit_wvanrees
#

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

HOME_MURPHY=/home/tgillis/murphy

cd ${HOME_MURPHY}
make destroy
ARCH_FILE=make_arch/make.engaging make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=0 -DBLOCK_GS=10"
mv murphy murphy60

# create the tmp directory
SCRATCH=/nobackup1/tgillis/perf_vs_acc
RUN_DIR=$SCRATCH/level-dbg2
mkdir -p ${RUN_DIR}
mkdir -p ${RUN_DIR}/data
mkdir -p ${RUN_DIR}/prof

# go there
cd ${RUN_DIR} 

cp ${HOME_MURPHY}/murphy60 .

#OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun ./murphy60 --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=3 --rtol=1e-2 --ctol=1e-4 --tfinal=0.1 --weno=5 > log_$SLURM_JOB_ID.log
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun valgrind --show-leak-kinds=all --leak-check=full --track-origins=yes --error-limit=no ./murphy60 --sadv --profile --iter-max=0 --iter-dump=10000 --level-min=0 --level-max=3 --ilevel=2 --rtol=1e-2 --ctol=1e-4 --tfinal=0.1 --weno=5 > log_$SLURM_JOB_ID.log



echo "------------------------"

