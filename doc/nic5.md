
## Build on NIC5 (intel)

```bash
module load releases/2020b
module load intel/2020b
module load HDF5

# p4est
wget https://p4est.github.io/release/p4est-2.3.2.tar.gz
tar -xvf p4est-2.3.2.tar.gz
cd p4est-2.3.2
CC=mpiicc CXX=mpiicpc ./configure --prefix=/home/ucl/tfl/tgillis/p4est-2.3.2-intel CFLAGS="-O3 -Wall" --enable-mpi --enable-openmp --with-blas=-mkl
make install -j


```

## Build on NIC5 (custom made OpenMPI and UCX)


```
# modules
module load releases/2020b foss/2020b HDF5/1.10.7-gompi-2020b


# UCX
wget https://github.com/openucx/ucx/releases/download/v1.10.1/ucx-1.10.1.tar.gz
tar -xvf ucx-1.10.1.tar.gz 
cd ucx-1.10.1
contrib/configure-release --prefix=/home/ucl/tfl/tgillis/lib-ompi-4.0.6/ucx-1.10.1 --enable-mt --with-avx --enable-compiler-opt=3
make install -j12

## OpenMPI (disable the verbs)
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.6.tar.gz
tar -xvf openmpi-4.0.6.tar.gz
./configure --prefix=/home/ucl/tfl/tgillis/lib-ompi-4.0.6/ompi-4.0.6 --with-ucx=/home/ucl/tfl/tgillis/lib-ompi-4.0.6/ucx-1.10.1 --without-verbs
make install -j12

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



```
configure: =========================================================
configure: UCX build configuration:
configure:       Build prefix:   /home/ucl/tfl/tgillis/lib-ompi-4.0.6/ucx-1.10.1
configure: Preprocessor flags:   -DCPU_FLAGS="|avx" -I${abs_top_srcdir}/src -I${abs_top_builddir} -I${abs_top_builddir}/src
configure:         C compiler:   gcc -O3 -g -Wall -Werror -mavx -funwind-tables -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-label -Wno-long-long -Wno-endif-labels -Wno-sign-compare -Wno-multichar -Wno-deprecated-declarations -Winvalid-pch -Wno-pointer-sign -Werror-implicit-function-declaration -Wno-format-zero-length -Wnested-externs -Wshadow -Werror=declaration-after-statement
configure:       C++ compiler:   g++ -O3 -g -Wall -Werror -mavx -funwind-tables -Wno-missing-field-initializers -Wno-unused-parameter -Wno-unused-label -Wno-long-long -Wno-endif-labels -Wno-sign-compare -Wno-multichar -Wno-deprecated-declarations -Winvalid-pch
configure:       Multi-thread:   enabled
configure:          MPI tests:   disabled
configure:      Devel headers:   no
configure:           Bindings:   < >
configure:        UCT modules:   < ib rdmacm cma knem >
configure:       CUDA modules:   < >
configure:       ROCM modules:   < >
configure:         IB modules:   < cm >
configure:        UCM modules:   < >
configure:       Perf modules:   < >
configure: =========================================================
```

```
Open MPI configuration:
-----------------------
Version: 4.0.6
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
HWLOC support: external
Libevent support: external
PMIx support: Internal

Transports
-----------------------
Cisco usNIC: no
Cray uGNI (Gemini/Aries): no
Intel Omnipath (PSM2): no
Intel TrueScale (PSM): no
Mellanox MXM: no
Open UCX: yes
OpenFabrics OFI Libfabric: yes
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
Generic Unix FS: yes
Lustre: no
PVFS2/OrangeFS: no
```