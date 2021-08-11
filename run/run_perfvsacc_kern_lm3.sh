#!/bin/bash
# Submission script for Engaging
#SBATCH --time=6:00:00
#
##SBATCH --ntasks-per-node=24
##SBATCH --mem-per-cpu=15625
#SBATCH --mem-per-cpu=2625
#SBATCH --partition=batch,debug
##SBATCH --partition=batch
#

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo "------------------------"

HOME_MURPHY=/home/ucl/tfl/tgillis/murphy

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

#run that shit
#mpirun ./${MNAME} --sadv --profile --iter-max=50000 --iter-dump=50000 --level-min=2 --level-max=9 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=0.5 --weno=5 --fix-weno --cfl-max=1.0 > log_$SLURM_JOB_ID.log
mpirun --mca pml ucx --mca osc ucx --mca btl ^uct \
        ./${MNAME} --sadv --profile --iter-max=50000 --iter-dump=50000 --level-min=2 --level-max=9 \
                   --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=0.5 --weno=5 --cfl-max=1.0 > log_$SLURM_JOB_ID.log
#mpirun --mca pml ucx --mca btl ^uct --tag-output valgrind --leak-check=full ./${MNAME} --sadv --profile --iter-max=50000 --iter-dump=50000 --level-min=0 --level-max=4 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=0.1 --weno=5 --cfl-max=1.0 > log_$SLURM_JOB_ID.log


echo "------------------------"

