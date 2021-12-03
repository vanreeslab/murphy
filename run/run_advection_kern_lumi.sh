#!/bin/bash
#SBATCH --partition=standard,debug
#SBATCH --time=4:00:00
##SBATCH --time=0:30:00
#SBATCH --ntasks-per-node=128
#SBATCH --account=project_465000036
#

echo "------------------------"
echo "welcome to the job: ${SLURM_JOB_NAME} -> id = ${SLURM_JOB_ID}"
echo ""
echo "RUN_DIR = ${RUN_DIR}"
echo "------------------------"

#HOME_MURPHY=/home/ucl/tfl/tgillis/murphy
HOME_MURPHY=/users/thgillis/murphy

# create the tmp directory
mkdir -p ${RUN_DIR}
mkdir -p ${RUN_DIR}/data
mkdir -p ${RUN_DIR}/prof

# go there
cd ${RUN_DIR} 
cp ${HOME_MURPHY}/${MNAME} .

#run that shit
#export MPICH_RMA_MAX_PENDING=2048
#export MPICH_RMA_MAX_PENDING=4096
#export MPICH_RMA_MAX_PENDING=8192
export MPICH_RMA_MAX_PENDING=32768
#export UCX_POSIX_USE_PROC_LINK=n
#export UCX_TLS=rc,ud,sm,self

ulimit -l

module list -t

#mpirun --mca mtl ofi -n ${SLURM_NPROCS} \
srun -u \
        ./${MNAME} --sadv --profile --iter-max=500000 --iter-dump=500000 --level-min=0 --level-max=12 \
                   --ilevel=${ILEVEL} --rtol=${RTOL} --ctol=${CTOL} ${NO_ADAPT} --tfinal=3.0 ${WENO} --cfl-max=0.25 > log_$SLURM_JOB_ID.log

echo "------------------------"

