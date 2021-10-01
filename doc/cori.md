
## Build on CORI (cray compilers)

```bash
#cray -> does not work due to a bug in GCC 8.1 https://github.com/TRIQS/triqs/issues/570
# module load cdt/20.10 PrgEnv-cray cray-hdf5-parallel craype-x86-skylake
#intel
module load PrgEnv-intel impi cray-hdf5-parallel

# p4est
wget https://p4est.github.io/release/p4est-2.3.2.tar.gz
tar -xvf p4est-2.3.2.tar.gz
cd p4est-2.3.2

# cray
#CC=cc CXX=CC F77=ftn FC=ftn ./configure --prefix=/global/homes/t/tgillis/p4est-2.3.2-cray6.0.9 CFLAGS="-Ofast -Wall -fopenmp" --enable-mpi --enable-openmp
#intel
CC=mpiicc CXX=mpiicpc F77=mpiifort FC=mpiifort ./configure --prefix=/global/homes/t/tgillis/p4est-2.3.2-intel19.0.3 CFLAGS="-Ofast -Wall -qopenmp" --enable-mpi --enable-openmp

make install -j
```
