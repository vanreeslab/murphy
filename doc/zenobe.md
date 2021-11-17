## Zenobe


```bash
module purge
module load OpenMPI/.4.1.2-GCC-8.3.0

export PREFIX=$HOME/lib-OpenMPI-4.1.2rc3-GCC-8.3.0
export CFLAGS="-march=ivybridge"
export CXXFLAGS="-march=ivybridge"

# hdf5
gzip -d hdf5-1.12.1.tar.gz
tar -xvf hdf5-1.12.1.tar
cd hdf5-1.12.1/
./configure --prefix=${PREFIX} --enable-parallel --enable-optimization=high --enable-build-mode=production --with-default-api-version=v110
make install -j


# ------------------------------------------------------------------------------
# p4est
wget https://p4est.github.io/release/p4est-2.3.2.tar.gz
tar -xvf p4est-2.3.2.tar.gz
cd p4est-2.3.2/
./configure --prefix=${PREFIX} --enable-mpi --enable-openmp
make install -j


```

:warning: when running a job, use `mpirun --mca osc pt2pt`, otherwise it will not work


<!-- ### OpenMPI
```bash
# ------------------------------------------------------------------------------
# load the needed modules
module load compiler/gcc/7.2.0
module load compiler/intel/comp_and_lib/2019.3.199

# define some stuff for an easy life
export PATH=$HOME/local/bin:$PATH
export PREFIX=$HOME/local
export CFLAGS="-march=ivybridge"
export CXXFLAGS="-march=ivybridge"
export CC=icc
export CXX=icpc

# ------------------------------------------------------------------------------
# OpenUCX
# v1.11.2
# wget https://github.com/openucx/ucx/releases/download/v1.11.2/ucx-1.11.2.tar.gz
# tar -xvf ucx-1.11.2.tar.gz
# cd ucx-1.11.2
# v1.10.1
wget https://github.com/openucx/ucx/releases/download/v1.10.1/ucx-1.10.1.tar.gz
tar -xvf ucx-1.10.1.tar.gz
cd ucx-1.10.1

./configure --prefix=--prefix=$PREFIX
./contrib/configure-release --prefix=$PREFIX --with-avx --enable-mt
make install -j

# ------------------------------------------------------------------------------
# OpenMPI - dependences
mkdir ompi_dep
cd ompi_dep

# M4
wget https://mirrors.tripadvisor.com/gnu/m4/m4-latest.tar.gz
tar -xvf m4-latest.tar.gz
cd m4-1.4.19/
./configure --prefix=$PREFIX
make install -j

# autoconf
wget https://mirrors.tripadvisor.com/gnu/autoconf/autoconf-latest.tar.gz
tar -xvf autoconf-latest.tar.gz
cd autoconf-2.71
./configure --prefix=$PREFIX
make install -j

# automake
wget https://mirrors.tripadvisor.com/gnu/automake/automake-1.16.5.tar.gz
tar -xvf automake-1.16.5.tar.gz
cd automake-1.16.5
./configure --prefix=$PREFIX
make install -j


# libtool
wget https://mirrors.tripadvisor.com/gnu/libtool/libtool-2.4.tar.gz
tar -xvf libtool-2.4.tar.gz
cd libtool-2.4
./configure --prefix=$PREFIX
make install -j

# ------------------------------------------------------------------------------
# OpenMPI
git clone git@github.com:open-mpi/ompi.git
cd ompi
git checkout v4.1.x
git submodule update --init --recursive

./autogen.pl
./configure --prefix=${PREFIX} --with-ucx=${PREFIX} --enable-mca-no-build=btl-uct --without-verbs
make install #-j
```
<!-- 
### p4est
```
wget https://p4est.github.io/release/p4est-2.3.2.tar.gz
tar -xvf p4est-2.3.2.tar.gz
cd p4est-2.3.2

./configure --prefix=${PREFIX} --enable-mpi --enable-openmp
``` -->

this is the result
```bash
-----------------------
Version: 4.1.2rc3
Build MPI C bindings: yes
Build MPI C++ bindings (deprecated): no
Build MPI Fortran bindings: mpif.h, use mpi, use mpi_f08
MPI Build Java bindings (experimental): no
Build Open SHMEM support: yes
Debug build: no
Platform file: (none)

Miscellaneous
-----------------------
CUDA support: no
HWLOC support: internal
Libevent support: internal
PMIx support: Internal

Transports
-----------------------
Cisco usNIC: no
Cray uGNI (Gemini/Aries): no
Intel Omnipath (PSM2): no
Intel TrueScale (PSM): no
Mellanox MXM: no
Open UCX: yes
OpenFabrics OFI Libfabric: no
OpenFabrics Verbs: no
Portals4: no
Shared memory/copy in+copy out: yes
Shared memory/Linux CMA: yes
Shared memory/Linux KNEM: no
Shared memory/XPMEM: no
TCP: yes

Resource Managers
-----------------------
Cray Alps: no
Grid Engine: no
LSF: no
Moab: no
Slurm: yes
ssh/rsh: yes
Torque: no

OMPIO File Systems
-----------------------
DDN Infinite Memory Engine: no
Generic Unix FS: yes
IBM Spectrum Scale/GPFS: yes
Lustre: no
PVFS2/OrangeFS: no
``` -->






