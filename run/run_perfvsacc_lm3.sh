# Let's gooooo

#remove the old stuffs
TAG=`date '+%Y-%m-%d'`-`uuidgen -t | head -c 8`
SCRATCH=$GLOBALSCRATCH/perf_vs_acc_lm3-${TAG}
echo "scratch file = ${SCRATCH}"

NCPUS=240

##==============================================================================
## WAVELET 6.x
##==============================================================================
## 6.0
#make destroy
#ARCH_FILE=make_arch/make.lm3 make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=0 -DBLOCK_GS=10"
#mv murphy murphy60
## 6.2
#make destroy
#ARCH_FILE=make_arch/make.lm3 make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=2 -DBLOCK_GS=12"
#mv murphy murphy62
##=============================
#i=2
#while [ $i -le 11 ]
#do
#	
#	# get the tolerances
#	export RTOL=$(echo "1e-$i")
#	export CTOL=$(echo "1e-$(($i+2))")
#	
#	# get the folder
#	export RUN_DIR=$SCRATCH/level$i
#	
#	# do not not-adapt
#	export NO_ADAPT=""
#	export ILEVEL=4
#	
#	# submit
#	export MNAME=murphy60
#	sbatch --ntasks=${NCPUS} --job-name=murphy60_${i} ./run/run_perfvsacc_kern_lm3.sh
#	export MNAME=murphy62
#	sbatch --ntasks=${NCPUS} --job-name=murphy62_${i} ./run/run_perfvsacc_kern_lm3.sh
#	
#	# increment the counter
#	((i++))
#	
#done

#==============================================================================
# WAVELET 2.x
#==============================================================================
# 2.0
make destroy
ARCH_FILE=make_arch/make.lm3 make -j12 OPTS="-DWAVELET_N=2 -DWAVELET_NT=0 -DBLOCK_GS=4"
mv murphy murphy20
# 2.2
make destroy
ARCH_FILE=make_arch/make.lm3 make -j12 OPTS="-DWAVELET_N=2 -DWAVELET_NT=2 -DBLOCK_GS=4"
mv murphy murphy22
#=============================
i=0
while [ $i -le 6 ]
do
	
	# get the tolerances
	export RTOL=$(echo "1e-$i")
	export CTOL=$(echo "1e-$(($i+2))")
	
	# get the folder
	export RUN_DIR=$SCRATCH/level$i
	
	# do not not-adapt
	export NO_ADAPT=""
	export ILEVEL=4
	
	# submit
	export MNAME=murphy20
	sbatch --ntasks=${NCPUS} --job-name=murphy20_${i} ./run/run_perfvsacc_kern_lm3.sh
	export MNAME=murphy22
	sbatch --ntasks=${NCPUS} --job-name=murphy22_${i} ./run/run_perfvsacc_kern_lm3.sh
	# increment the counter
	((i++))
	
done

#==============================================================================
# WAVELET 4.x
#==============================================================================
# 4.0
make destroy
ARCH_FILE=make_arch/make.lm3 make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy40
# 4.2
make destroy
ARCH_FILE=make_arch/make.lm3 make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=8"
mv murphy murphy42
#=============================
i=0
while [ $i -le 6 ]
do
	
	# get the tolerances
	export RTOL=$(echo "1e-$i")
	export CTOL=$(echo "1e-$(($i+2))")
	
	# get the folder
	export RUN_DIR=$SCRATCH/level$i
	
	# do not not-adapt
	export NO_ADAPT=""
	export ILEVEL=4
	
	# submit
	export MNAME=murphy40
	sbatch --ntasks=${NCPUS} --job-name=murphy40_${i} ./run/run_perfvsacc_kern_lm3.sh
	export MNAME=murphy42
	sbatch --ntasks=${NCPUS} --job-name=murphy42_${i} ./run/run_perfvsacc_kern_lm3.sh
	# increment the counter
	((i++))
	
done


##==============================================================================
## UNIFORM
##==============================================================================
#i=2
#while [ $i -le 8 ]
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
#	#bash ./run/run_perfvsacc_kern_lm3.sh
#	export MNAME=murphy20
#	sbatch --ntasks=${NCPUS} --job-name=murphy_ref_${i} ./run/run_perfvsacc_kern_lm3.sh
#	#export MNAME=murphy42
#	#sbatch --ntasks=${NCPUS} --job-name=murphy_ref_${i} ./run/run_perfvsacc_kern_lm3.sh
#	
#	# increment the counter
#	((i++))
#	
#done





