## Build on Engaging

OpenMPI + UCX [information page](https://github.com/openucx/ucx/wiki/OpenMPI-and-OpenSHMEM-installation-with-UCX)

```bash
# modules
module purge
module load gcc/11.2.0 OpenBLAS/0.2.0

export INSTALL_DIR=/home/tgillis/lib-ucx-1.11.0-ompi-4.0.6

# UCX
#wget https://github.com/openucx/ucx/releases/download/v1.10.1/ucx-1.10.1.tar.gz
#tar -xvf ucx-1.10.1.tar.gz 
#cd ucx-1.10.1
#contrib/configure-release --prefix=/home/tgillis/lib-ompi-4.0.6/ucx-1.10.1
#make install -j12
wget https://github.com/openucx/ucx/releases/download/v1.11.0/ucx-1.11.0.tar.gz
tar -xvf ucx-1.11.0.tar.gz
cd ucx-1.11.0
mkdir build
cd build
../configure --with-march=native --prefix=${INSTALL_DIR}/ucx-1.11.0
make install -j

## OpenMPI (disable the verbs)
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.6.tar.gz
tar -xvf openmpi-4.0.6.tar.gz
./configure --prefix=${INSTALL_DIR}/ompi-4.0.6 --with-ucx=${INSTALL_DIR}/ucx-1.11.0 --without-verbs
make install -j12

export MY_MPI_DIR=/home/tgillis/lib-ucx-1.11.0-ompi-4.0.6/ompi-4.0.6

# p4est 2.3.2
wget https://p4est.github.io/release/p4est-2.3.2.tar.gz
tar -xvf p4est-2.3.2.tar.gz
cd p4est-2.3.2/
F77=${MY_MPI_DIR}/bin/mpif77 FC=${MY_MPI_DIR}/bin/mpif90 CC=${MY_MPI_DIR}/bin/mpicc CXX=${MY_MPI_DIR}/bin/mpic++ ./configure --prefix=${INSTALL_DIR}/p4est-2.3.2 CFLAGS="-march=native -O3 -Wall" --enable-mpi --enable-openmp --with-blas=-lopenblas
make install -j

#p4est
#export MY_MPI_DIR=/home/tgillis/lib-ompi-4.0.6/ompi-4.0.6
#F77=${MY_MPI_DIR}/bin/mpif77 FC=${MY_MPI_DIR}/bin/mpif90 CC=${MY_MPI_DIR}/bin/mpicc CXX=${MY_MPI_DIR}/bin/mpic++ ./configure --prefix=/home/tgillis/lib-ompi-4.0.6/p4est CFLAGS="-march=skylake -O3 -Wall" --enable-mpi --enable-openmp --with-blas=-lopenblas


# hdf5
CC=${MY_MPI_DIR}/bin/mpicc CXX=${MY_MPI_DIR}/bin/mpic++ CFLAGS="-march=native" CXXFLAGS="-march=native" ./configure --prefix=${INSTALL_DIR}/hdf5-1.10.6 --enable-parallel --enable-optimization=high --enable-build-mode=production
make install -j
```