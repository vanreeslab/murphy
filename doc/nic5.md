## Build on NIC5


```
# modules
module load releases/2020b foss/2020b HDF5/1.10.7-gompi-2020b

# p4est
git clone https://github.com/cburstedde/p4est.git
cd p4est/
git checkout feature-mesh-edge
git submodule init
git submodule update
./bootstrap
CC=mpic++ ./configure --prefix=/home/ucl/tfl/tgillis/p4est-lib CFLAGS="-O3 -Wall" --enable-mpi --enable-openmp --with-blas=-lopenblas
make install -j
 ```


:warning: OpenMPI need some special configuration to make the most out of the UCX lib:

```
mpirun --mca pml ucx --mca osc ucx --mca btl ^uct
```
