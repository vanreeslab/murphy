#!/bin/bash
# Submission script for Zenobe 
#PBS -N murphy
#PBS -q large
##PBS -q main
#PBS -W group_list=examples
#PBS -l walltime=06:00:00 
#

exec > ${PBS_O_WORKDIR}/${PBS_JOBNAME}_${PBS_JOBID}.log 
echo "------------------ Work dir --------------------" 
cd ${PBS_O_WORKDIR} && echo ${PBS_O_WORKDIR} 
echo "------------------ Job Info --------------------" 
echo "jobid : $PBS_JOBID" 
echo "jobname : $PBS_JOBNAME" 
echo "job type : $PBS_ENVIRONMENT" 
echo "submit dir : $PBS_O_WORKDIR" 
echo "queue : $PBS_O_QUEUE" 
echo "user : $PBS_O_LOGNAME" 
echo "threads : $OMP_NUM_THREADS" 

module purge
module load OpenMPI/.4.1.2-GCC-8.3.0

HOME_MURPHY=/home/acad/ucl-tfl/tgillis/murphy

# create the tmp directory
mkdir -p ${RUN_DIR}
mkdir -p ${RUN_DIR}/data
mkdir -p ${RUN_DIR}/prof

# go there
cd ${RUN_DIR}
cp ${HOME_MURPHY}/${MNAME} .

# run that shit
mpirun --version
mpirun --mca osc pt2pt ${MNAME} --sadv --profile --iter-max=50000 --iter-dump=50000 --level-min=0 --level-max=12 --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=3.0 ${WENO} --cfl-max=0.25 > log_${PBS_JOBID}.log


qstat -f $PBS_JOBID
