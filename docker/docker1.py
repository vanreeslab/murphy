## choose ubuntu
Stage0 += baseimage(image='ubuntu:20.04')

# install some basic stuffs for both stages
apt = apt_get(ospackages=['git', 'wget', 'apt-transport-https', 'ca-certificates',
                          'autotools-dev', 'autoconf', 'libtool', 'automake',
                          'vim',
                          'make', 'cmake',
                          'valgrind',
                          'm4',
                          'python'])
# clang
cc = llvm(version='12', openmp=True, extra_tools=True, toolset=True,
          environment=True)

# fotran compiler from gnu !only fortran and no c/c++!
fc = gnu(version='9', cc=False, cxx=False, fortran=True, environment=True)

# add all that to Stage0
Stage0 += apt
Stage0 += fc
Stage0 += cc

# openmpi 4.1.1 with ucx
Stage0 += ofed(toolchain=cc.toolchain)
Stage0 += ucx(version='1.11.2', cuda=False, toolchain=cc.toolchain)

omp = openmpi(version='4.1.1', cuda=False, ucx=True,
              toolchain=cc.toolchain, configure_opts=['--disable-mpi-fortran'])
Stage0 += omp

# hdf5
Stage0 += hdf5(version='1.10.5', toolchain=omp.toolchain,
               configure_opts=['--disable-fortran', '--enable-parallel', '--enable-optimization=high', '--enable-build-mode=production'])

#fftw
Stage0 += fftw(version='3.3.9', toolchain=omp.toolchain,
               configure_opts=['--disable-fortran', '--enable-openmp', '--enable-sse2', '--enable-avx'])

# blas and lapack
Stage0 += openblas(version='0.3.18', toolchain=omp.toolchain,
                   make_opts=['USE_OPENMP=1', 'DYNAMIC_ARCH=1', 'CC=clang'])

# build stage 1
Stage1 += baseimage(image='ubuntu:20.04')
Stage1 += apt
Stage1 += fc
Stage1 += cc
Stage1 += Stage0.runtime()
