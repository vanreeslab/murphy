
## Build on LUMI

for p4est:
```bash
wget https://p4est.github.io/release/p4est-2.3.2.tar.gz
tar -xvf p4est-2.3.2.tar.gz
cd p4est-2.3.2
CC=cc CXX=CC F77=ftn FC=ftn ./configure --prefix=/global/homes/t/tgillis/p4est-2.3.2-gcc10.1.0 CFLAGS="-Ofast -Wall -fopenmp" --enable-mpi --enable-openmp
make install -j
```

:warning: don't forget to change the `prefix`.


murphy: change and/or update the `make_arch/make.lumi` arch file and then

```bash
ARCH_FILE=make_arch/make.lumi make -j
```

### scalability runs
```bash
./murphy --weak-scal --profile
```

the domain created = the work load will have an aspect ratio of NCPU/32 x 1 x 1 and we run 50 times each operation (cfr the profiling)