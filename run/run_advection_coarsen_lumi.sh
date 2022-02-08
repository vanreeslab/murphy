# Let's gooooo

#remove the old stuffs
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
MURPHY_SCRATCH=$SCRATCH/lumi_murphy-adv-coarsen-${TAG}
echo "scratch file = ${MURPHY_SCRATCH}"

#==============================================================================
# NCPUS
#NCPUS=512

##==============================================================================
# compilation
# 4.0
make clean
ARCH_FILE=make_arch/make.lumi make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy40
# 4.2
make clean
ARCH_FILE=make_arch/make.lumi make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=8"
mv murphy murphy42
# 6.0
make clean
ARCH_FILE=make_arch/make.lumi make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=0 -DBLOCK_GS=10"
mv murphy murphy60
#### 6.2
####make clean
##ARCH_FILE=make_arch/make.lumi make -j12 OPTS="-DWAVELET_N=6 -DWAVELET_NT=2 -DBLOCK_GS=12"
##mv murphy murphy62
#=============================
i=4
while [ $i -le 4 ]
do
    export NCPUS_TMP=$(echo "$(( ($i+1)*128 ))")
    export NCPUS=$(echo "$(( $NCPUS_TMP > 0 ? $NCPUS_TMP : 128))")
    #export NCPUS_TMP=$(echo "$(( ($i+1)*128*5 ))")
    #export NCPUS=$(echo "$(( $NCPUS_TMP > 0 ? $NCPUS_TMP : 384))")
	# get the folder
	export RUN_DIR=${MURPHY_SCRATCH}/level$i
    echo "submission in ${RUN_DIR}"

    #--------------------------------------------------------------------------
    # start the "with compression" -> W4.0
	export NO_ADAPT="--iter-adapt=6"
	export ILEVEL=$i
    for r in {1,2,4,16}
    do
        echo "submitting r = ${r}"
        export RTOL=$(bc -l <<< "scale=16; 10^(-${i}-2)")
        export CTOL=$(bc -l <<< "scale=16; 10^(-${i}-2-${r})")

        # 4.0
	    export MNAME=murphy40_c3_${i}_a6_t${r}
        cp murphy40 ${MNAME}
        export WENO="--weno=3 --fix-weno"
	    sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
        #export WENO="--weno=5 --fix-weno"
	    #sbatch --ntasks=${NCPUS} --job-name=${MNAME}_c5_${i}_a6_t${r} ./run/run_advection_kern_lumi.sh

        # 4.2
	    export MNAME=murphy42_c3_${i}_a6_t${r}
        cp murphy42 ${MNAME}
        export WENO="--weno=3 --fix-weno"
	    sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
        #export WENO="--weno=5 --fix-weno"
	    #sbatch --ntasks=${NCPUS} --job-name=${MNAME}_c5_${i}_a6_t${r} ./run/run_advection_kern_lumi.sh
    done
    ##--------------------------------------------------------------------------
    echo "submission 60"
	export NO_ADAPT="--iter-adapt=6"
	export ILEVEL=$i
    export RTOL=$(bc -l <<< "scale=16; 10^(-${i}-2)")
    export CTOL=$(bc -l <<< "scale=16; 10^(-${i}-4)")
	export MNAME=murphy60_c5_${i}_a6_t2
    cp murphy60 ${MNAME}
    export WENO="--weno=5 --fix-weno"
	sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
    
    ##--------------------------------------------------------------------------
	# refine only
	export NO_ADAPT="--iter-adapt=6"
	export ILEVEL=0
    export RTOL=$(bc -l <<< "scale=16; 10^(-$i-2)")
    export CTOL=-1.0
	# 4.0
	export MNAME=murphy40_c3_${i}_r6_t0
    cp murphy40 ${MNAME}
    export WENO="--weno=3 --fix-weno"
	sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
    #export WENO="--weno=5 --fix-weno"
    #sbatch --ntasks=${NCPUS} --job-name=${MNAME}_c5_${i}_r6_t0 ./run/run_advection_kern_lumi.sh

	# 4.2
	export MNAME=murphy42_c3_${i}_r6_t0 
    cp murphy42 ${MNAME}
    export WENO="--weno=3 --fix-weno"
	sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
    #export WENO="--weno=5 --fix-weno"
    #sbatch --ntasks=${NCPUS} --job-name=${MNAME}_c5_${i}_r6_t0 ./run/run_advection_kern_lumi.sh

    ##--------------------------------------------------------------------------
	# do not adapt
	export NO_ADAPT="--no-adapt"
	export ILEVEL=$i
    export RTOL=1.0
    export CTOL=0.0
	# 4.0
	export MNAME=murphy40_c3_${i}_a0_t0
    cp murphy40 ${MNAME}
    export WENO="--weno=3 --fix-weno"
	sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
    #export WENO="--weno=5 --fix-weno"
	#sbatch --ntasks=${NCPUS} --job-name=${MNAME}_w5_${i}_a0_t0 ./run/run_advection_kern_lumi.sh

    # 4.2
	export MNAME=murphy42_c3_${i}_a0_t0
    cp murphy42 ${MNAME}
    export WENO="--weno=3 --fix-weno"
	sbatch --ntasks=${NCPUS} --job-name=${MNAME} ./run/run_advection_kern_lumi.sh
    #export WENO="--weno=5 --fix-weno"
	#sbatch --ntasks=${NCPUS} --job-name=${MNAME}_w5_${i}_a0_t0 ./run/run_advection_kern_lumi.sh

    echo "done with i = ${i}"

	# increment the counter
	((i++))
	
done






