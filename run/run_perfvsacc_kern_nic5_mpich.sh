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
mkdir -p ${RUN_DIR}/data
mkdir -p ${RUN_DIR}/prof

# go there
cd ${RUN_DIR} 
cp ${HOME_MURPHY}/${MNAME} .

# get the machine file
MACHINEFILE="nodes.$SLURM_JOBID"
srun -l /bin/hostname | sort -n | awk '{print $2}' > $MACHINEFILE

#run that shit
#mpirun ./${MNAME} --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} > log_$SLURM_JOB_ID.log
/home/ucl/tfl/tgillis/lib-mpich/mpich-3.4.1/bin/mpirun -f ${MACHINEFILE} ./${MNAME} --sadv --profile --iter-max=5000 --iter-dump=10000 --level-min=2 --level-max=7 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=0.1 --weno=5 --cfl-max=1.0 > log_$SLURM_JOB_ID.log


echo "------------------------"

