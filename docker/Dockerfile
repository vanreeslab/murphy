FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && apt-get install -y \
    apt-utils \
    gcc g++ gfortran \
    make git vim wget unzip cmake valgrind \
    automake autoconf m4 libblas-dev liblapack-dev

# OpenMPI v4.0.5
RUN wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.5.tar.gz &&\
    tar -xvf openmpi-4.0.5.tar.gz &&\
    cd openmpi-4.0.5 &&\
    ./configure --prefix=/usr/local --enable-debug --enable-memchecker &&\
    make all install

RUN apt-get -y update && apt-get install -y \
    libhdf5-mpi-dev libfftw3-dev \
    python3 python3-pip

#p4est
RUN git clone --recursive https://github.com/cburstedde/p4est
RUN cd /p4est &&\
    git checkout feature-mesh-edge &&\
    git submodule init && git submodule update &&\
    ./bootstrap &&\
    CC=mpicc CXX=mpic++ ./configure --prefix=/soft/p4est-github CFLAGS="-O3 -Wall" --enable-mpi --enable-openmp --enable-valgrind-clean &&\
    make install -j
RUN rm -rf /p4est

#flups
# flups is downloaded locally mannually RUN git clone https://git.immc.ucl.ac.be/examples/flups.git
COPY flups /flups
RUN cd /flups &&\
    ARCH_FILE=make_arch/make.docker_valid_opt PREFIX=/soft/flups make install -j
RUN rm -rf /flups

#Google test
RUN git clone https://github.com/google/googletest.git
RUN cd /googletest &&\
    cmake . -DCMAKE_INSTALL_PREFIX=/soft/googletest &&\
    make install -j
RUN rm -rf /googletest


# add an openMPI setting to make it run                                         
ENV OMPI_MCA_btl_vader_single_copy_mechanism none  

#cpp lint
#RUN pip3 install cpplint
