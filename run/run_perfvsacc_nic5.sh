# Let's gooooo

#remove the old stuffs
SCRATCH=$GLOBALSCRATCH/perf_vs_acc_nic5
rm -rf ${SCRATCH}

## produce the different exec
## 2.0
#make destroy
#ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=2 -DWAVELET_NT=0 -DBLOCK_GS=2"
#mv murphy murphy20
## 2.2
#make destroy
#ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=2 -DWAVELET_NT=2 -DBLOCK_GS=4"
#mv murphy murphy22
# 4.0
make destroy
ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy40
## 4.2
#make destroy
#ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=8"
#mv murphy murphy42
## 6.0
#make destroy
#ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=0 -DBLOCK_GS=10"
#mv murphy murphy60
## 6.2
#make destroy
#ARCH_FILE=make_arch/make.nic5_mpich make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=2 -DBLOCK_GS=12"
#mv murphy murphy62


NCPUS=384

# let's go for the adaptive simulations
i=5
while [ $i -le 7 ]
do
	
	# get the tolerances
	export RTOL=$(echo "1e-$i")
	export CTOL=$(echo "1e-$(($i+2))")
	
	# get the folder
	export RUN_DIR=$SCRATCH/level$i
	
	# do not not-adapt
	export NO_ADAPT=""
	export ILEVEL=3
	
	# submit
	#bash ./run/run_perfvsacc_kern_nic5_mpich.sh
	#export MNAME=murphy20
	#sbatch --ntasks=${NCPUS} --job-name=murphy20_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
	#export MNAME=murphy22
	#sbatch --ntasks=${NCPUS} --job-name=murphy22_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
	export MNAME=murphy40
	sbatch --ntasks=${NCPUS} --job-name=murphy40_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
	#export MNAME=murphy42
	#sbatch --ntasks=${NCPUS} --job-name=murphy42_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
	#export MNAME=murphy60
	#sbatch --ntasks=${NCPUS} --job-name=murphy60_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
	#export MNAME=murphy62
	#sbatch --ntasks=${NCPUS} --job-name=murphy62_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
	
	# increment the counter
	((i++))
	
done


#i=3
#while [ $i -le 7 ]
#do
#	
#	# get the tolerances
#	export RTOL=$(echo "10^(-($i+0))" | bc -l)
#	export CTOL=$(echo "10^(-($i+2))" | bc -l)
#	
#	# get the folder
#	export RUN_DIR=$SCRATCH/ref_level$i
#	
#	# do not not-adapt
#	export NO_ADAPT="--no-adapt"
#	export ILEVEL=$i
#	
#	# submit
#	#bash ./run/run_perfvsacc_kern_nic5_mpich.sh
#	export MNAME=murphy42
#	sbatch --ntasks=${NCPUS} --job-name=murphy_ref_${i} ./run/run_perfvsacc_kern_nic5_mpich.sh
#	
#	# increment the counter
#	((i++))
#	
#done





