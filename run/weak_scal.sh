# loop over the 70 nodes

NODE_MAX=32

for INODE in $(seq 1 $NODE_MAX)
do
    # submit the job
    sbatch --ntasks=$(($INODE*64)) run/run_nic5.sh
done
