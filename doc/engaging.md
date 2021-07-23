## Build on Engaging

OpenMPI + UCX [information page](https://github.com/openucx/ucx/wiki/OpenMPI-and-OpenSHMEM-installation-with-UCX)

```bash
# modules
module purge
module load gcc/9.3.0 OpenBLAS/0.2.0

# UCX
wget https://github.com/openucx/ucx/releases/download/v1.10.1/ucx-1.10.1.tar.gz
tar -xvf ucx-1.10.1.tar.gz 
cd ucx-1.10.1
contrib/configure-release --prefix=/home/tgillis/lib-ompi-4.0.6/ucx-1.10.1
make install -j12

## OpenMPI (disable the verbs)
./configure --prefix=/home/tgillis/lib-ompi-4.0.6/ompi-4.0.6 --with-ucx=/home/tgillis/lib-ompi-4.0.6/ucx-1.10.1 --without-verbs
make install -j12

#p4est
export MY_MPI_DIR=/home/tgillis/lib-ompi-4.0.6/ompi-4.0.6
F77=${MY_MPI_DIR}/bin/mpif77 FC=${MY_MPI_DIR}/bin/mpif90 CC=${MY_MPI_DIR}/bin/mpicc CXX=${MY_MPI_DIR}/bin/mpic++ ./configure --prefix=/home/tgillis/lib-ompi-4.0.6/p4est CFLAGS="-march=skylake -O3 -Wall" --enable-mpi --enable-openmp --with-blas=-lopenblas

# hdf5
CC=${MY_MPI_DIR}/bin/mpicc CXX=${MY_MPI_DIR}/bin/mpic++ CFLAGS="-march=skylake" CXXFLAGS="-march=skylake" ./configure --prefix=/home/tgillis/lib-ompi-4.0.6/hdf5-1.10.6 --enable-parallel --enable-optimization=high --enable-build-mode=production
make install -j
```