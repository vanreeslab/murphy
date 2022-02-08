
## Build on LUMI

#### modules

load the modules:
```bash
module load cray-libsci
module load cray-hdf5-parallel/1.12.0.6
module load cray-mpich-ucx/8.1.8
module load cray-ucx/2.7.0-1
```

:warning: **NEVER** run `module purge`, logout/login instead

#### p4est

```bash
wget https://p4est.github.io/release/p4est-2.3.3.tar.gz
tar -xvf p4est-2.3.3.tar.gz
cd p4est-2.3.3
CC=cc CXX=CC F77=ftn FC=ftn ./configure --prefix=${HOME}/lib-cray-8.0.0-mpich-ucx-8.1.8 CFLAGS="-Ofast -Wall -fopenmp" --enable-mpi --enable-openmp
make install -j
```

:warning: add the following variable to you bashrc and your submission scripts: `MPICH_RMA_MAX_PENDING=512` :-)

#### murphy
murphy: change and/or update the `make_arch/make.lumi` arch file and then

```bash
ARCH_FILE=make_arch/make.lumi make -j
```

#### scalability runs
```bash
./murphy --weak-scal --profile
```

the domain created = the work load will have an aspect ratio of NCPU/32 x 1 x 1 and we run 50 times each operation (cfr the profiling)
