FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

#RUN echo "nameserver 8.8.4.4" >> /etc/resolv.conf
#RUN apt update
#RUN apt upgrade -y
RUN apt-get -y update && apt-get install -y \
    make git gfortran gcc g++ linux-libc-dev vim wget unzip cmake valgrind \
    automake autoconf m4 libblas-dev liblapack-dev \
    libhdf5-mpi-dev libfftw3-dev libopenmpi-dev openmpi-bin \
    python3 python3-pip \
    gcovr

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

#cpp lint
#RUN pip3 install cpplint
