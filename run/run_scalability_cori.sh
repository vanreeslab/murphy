#remove the old stuffs
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
export RUN_DIR=${SCRATCH}/weak_cori_${TAG}
echo "RUN DIR = ${RUN_DIR}"

# compile murphy
#make destroy
ARCH_FILE=make_arch/make.cori make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"

#create the dirs
mkdir -p ${RUN_DIR}/prof
mkdir -p ${RUN_DIR}/data

# move what is needed
cp murphy ${RUN_DIR}/
cp ./run/run_scalability_kern_cori.sh ${RUN_DIR}/



cd ${RUN_DIR}

# let's go for the adaptive simulations
i=1
while [ $i -le 1 ]
do
	
	# submit
	echo "sbatch --nodes=${i} --job-name=weak_${i} .run_scalability_kern_cori.sh"
	sbatch --nodes=${i} --job-name=weak_${i} ./run_scalability_kern_cori.sh
	
	# increment the counter
	((i++))
done




