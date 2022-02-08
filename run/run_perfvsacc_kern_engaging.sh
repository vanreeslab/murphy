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

# create the tmp directory
mkdir -p ${RUN_DIR}
mkdir -p ${RUN_DIR}/data
mkdir -p ${RUN_DIR}/prof

# go there
cd ${RUN_DIR} 

cp ${HOME_MURPHY}/${MNAME} .

# get the machine file                                                          
MACHINEFILE="nodes.$SLURM_JOBID"                                                
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE 

#compile the 4.2 exec
#echo "mpirun ./murphy --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} > log_$SLURM_JOB_ID.log"
#mpirun valgrind ./murphy --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} > log_$SLURM_JOB_ID.log
#mpirun --mca osc pt2pt ./murphy --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} > log_$SLURM_JOB_ID.log
#OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun -f ${MACHINEFILE} valgrind --track-origins=yes --error-limit=no ./${MNAME} --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=0.1 --weno=5 > log_$SLURM_JOB_ID.log
OMP_NUM_THREADS=1 /home/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun -env UCX_NET_DEVICES=mlx5_0:1 ./${MNAME} --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=0 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=0.1 --weno=5 --cfl-max=0.75 > log_$SLURM_JOB_ID.log


echo "------------------------"

