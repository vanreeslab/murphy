#remove the old stuffs
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
export RUN_DIR=${GLOBALSCRATCH}/weak_nic5_${TAG}

# compile murphy
#make destroy
ARCH_FILE=make_arch/make.nic5_intel make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"

echo "-------------------------------------------------------------------------"
echo "RUN DIR = ${RUN_DIR}"
#create the dirs
mkdir -p ${RUN_DIR}/prof
mkdir -p ${RUN_DIR}/data

# move what is needed
cp murphy ${RUN_DIR}/
cp ./run/run_scalability_kern_nic5.sh ${RUN_DIR}/



cd ${RUN_DIR}

# let's go for the adaptive simulations
#submit the first
i=1
sbatch --nodes=${i} --job-name=weak_${i} ./run_scalability_kern_nic5.sh

# start the loop
i=4
#while [ $i -le 64 ]
while [ $i -le 10 ]
do
	# submit
	sbatch --nodes=${i} --job-name=weak_${i} ./run_scalability_kern_nic5.sh
	
	# increment the counter
	((i+=4))
done




