# Let's gooooo

#remove the old stuffs
TAG=`date '+%Y-%m-%d-%H%M'`-`uuidgen -t | head -c 8`
SCRATCH=/SCRATCH/acad/examples/tgillis/advection-${TAG}

#==============================================================================
# compilation
# 4.0
#make destroy
make clean
ARCH_FILE=make_arch/make.zenobe make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=0 -DBLOCK_GS=6"
mv murphy murphy40
make clean
ARCH_FILE=make_arch/make.zenobe make -j12 OPTS="-DWAVELET_N=4 -DWAVELET_NT=2 -DBLOCK_GS=8"
mv murphy murphy42

#==============================================================================
echo "scratch file = ${SCRATCH}"

#==============================================================================
i=0
while [ $i -le 5 ]
do
    # get the number of chunks of 24 cpus
    NCPUS_TMP=$(echo "$(($i*4))")
    NCPUS=$(echo "$(( $NCPUS_TMP > 0 ? $NCPUS_TMP : 4))")
    # get folder
	export RUN_DIR=$SCRATCH/level$i

    echo "submission in ${RUN_DIR}"

    #--------------------------------------------------------------------------
    echo "submission in ${RUN_DIR}: adapt"
    # start the "with compression"
	export NO_ADAPT="--iter-adapt=6"
	export ILEVEL=$i
    for r in {1,2,4,16}
    do
        export RTOL=$(printf "%e" $(bc -l <<< "scale=16; 10^(-${i}-2)"))
        export CTOL=$(printf "%e" $(bc -l <<< "scale=16; 10^(-${i}-2-${r})"))
        export WENO="--weno=3 --fix-weno"
        #run the shit

        #4.0
        export MNAME=murphy40
        echo "qsub -v MNAME=murphy40,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
             -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh"
        qsub -v MNAME=murphy40,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
             -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh

        #4.2
        export MNAME=murphy42
        echo "qsub -v MNAME=murphy42,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
             -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh"
        qsub -v MNAME=murphy42,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
             -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh
    done

    ###--------------------------------------------------------------------------
    echo "submission in ${RUN_DIR}: refinement only"
	# refine only
	export NO_ADAPT="--iter-adapt=6"
	export ILEVEL=$i
    export RTOL=$(printf "%e" $(bc -l <<< "scale=16; 10^(-${i}-2)"))
    export CTOL=-1.0
    export WENO="--weno=3 --fix-weno"
    #4.0
    export MNAME=murphy40
    echo "qsub -v MNAME=murphy40,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh"
    qsub -v MNAME=murphy40,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh

    #4.2
    export MNAME=murphy42
    echo "qsub -v MNAME=murphy42,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh"
    qsub -v MNAME=murphy42,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh

    ##--------------------------------------------------------------------------
    echo "submission in ${RUN_DIR}: no-adapt"
	# do not adapt
	export NO_ADAPT="--no-adapt"
	export ILEVEL=0
    export RTOL=1.0
    export CTOL=0.0
    export WENO="--weno=3 --fix-weno"

    #4.0
    export MNAME=murphy40
    echo "qsub -v MNAME=murphy40,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh"
    qsub -v MNAME=murphy40,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh

    #4.2
    export MNAME=murphy42
    echo "qsub -v MNAME=murphy42,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh"
    qsub -v MNAME=murphy42,RUN_DIR="${RUN_DIR}",WENO="${WENO}",RTOL=${RTOL},CTOL=${CTOL},NO_ADAPT="${NO_ADAPT}",ILEVEL=${ILEVEL} \
         -l select=${NCPUS}:ncpus=24:mem=63000mb:mpiprocs=24:ompthreads=1 -N ${MNAME}_c3_${i}_a6_t${r} ./run/run_advection_kern_zenobe.sh

	# increment the counter
	((i++))
	
done

