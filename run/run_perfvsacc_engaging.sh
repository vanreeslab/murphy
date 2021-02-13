# Let's gooooo

#remove the old stuffs
SCRATCH=/nobackup1/tgillis/perf_vs_acc
rm -rf ${SCRATCH}

# produce the different exec
# 4.0
make destroy
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=4"
mv murphy murphy40
# 4.2
make destroy
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=4"
mv murphy murphy42
#
make destroy
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=4 -DWAVELET_NT=4 -DBLOCK_GS=6"
mv murphy murphy44
#
make destroy
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=6 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy60
#
make destroy
ARCH_FILE=make_arch/make.engaging make -j OPTS="-DWAVELET_N=6 -DWAVELET_NT=2 -DBLOCK_GS=6"
mv murphy murphy62


NCPUS=216

# let's go for the adaptive simulations
i=0
while [ $i -le 7 ]
do
	
	# get the tolerances
	export RTOL=$(echo "1e-$i")
	export CTOL=$(echo "1e-$(($i+2))")
	#export RTOL=$(echo "10^(-($i+0))" | bc -l | sed 's!\.0*$!!')
	#export CTOL=$(echo "10^(-($i+2))" | bc -l | sed 's!\.0*$!!')
	
	# get the folder
	export RUN_DIR=$SCRATCH/level$i
	
	# do not not-adapt
	export NO_ADAPT=""
	export ILEVEL=5
	
	# submit
	#bash ./run/run_perfvsacc_kern_engaging.sh
	export MNAME=murphy40
	sbatch --ntasks=${NCPUS} --job-name=murphy40_${i} ./run/run_perfvsacc_kern_engaging.sh
	export MNAME=murphy42
	sbatch --ntasks=${NCPUS} --job-name=murphy42_${i} ./run/run_perfvsacc_kern_engaging.sh
	export MNAME=murphy44
	sbatch --ntasks=${NCPUS} --job-name=murphy44_${i} ./run/run_perfvsacc_kern_engaging.sh
	export MNAME=murphy60
	sbatch --ntasks=${NCPUS} --job-name=murphy60_${i} ./run/run_perfvsacc_kern_engaging.sh
	export MNAME=murphy62
	sbatch --ntasks=${NCPUS} --job-name=murphy62_${i} ./run/run_perfvsacc_kern_engaging.sh
	
	# increment the counter
	((i++))
	
done


i=3
while [ $i -le 7 ]
do
	
	# get the tolerances
	export RTOL=$(echo "10^(-($i+0))" | bc -l)
	export CTOL=$(echo "10^(-($i+2))" | bc -l)
	
	# get the folder
	export RUN_DIR=$SCRATCH/ref_level$i
	
	# do not not-adapt
	export NO_ADAPT="--no-adapt"
	export ILEVEL=$i
	
	# submit
	#bash ./run/run_perfvsacc_kern_engaging.sh
	export MNAME=murphy42
	sbatch --ntasks=${NCPUS} --job-name=murphy_ref_${i} ./run/run_perfvsacc_kern_engaging.sh
	
	# increment the counter
	((i++))
	
done





