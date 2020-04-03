FROM ubuntu:18.04


RUN echo "nameserver 8.8.4.4" >> /etc/resolv.conf
#RUN apt update
#RUN apt upgrade -y
RUN apt-get -y update && apt-get install -y \
make git g++ linux-libc-dev vim wget unzip cmake \
libhdf5-mpi-dev libfftw3-dev libopenmpi-dev openmpi-bin \
python3 python3-pip

RUN pip3 install cpplint

VOLUME ["/etc/gitlab-runner"]

RUN useradd -ms /bin/bash validator
USER validator
WORKDIR /home/validator

