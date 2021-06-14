# How to build the containers needed by murphy

Two containers are involved in the process: a base one build with the basics (`C++` compilers, `MPI`,...) and a murphy-dedicated one build on top of it using some special libs.

### base container
To build the first container you need HPCCM.
To install it, simply use `sudo pip install hpccm` (or `pip3` instead of `pip`).

Then, get the Dockerfile and build the container:
```bash
hpccm --recipe mpi_bandwidth.py > Dockerfile_base
docker build -f Dockerfile_base -t vanreeslab/murphy_base .

# register and put is on DockerHub
docker login
docker push vanreeslab/murphy_base

# to test is
docker run --rm -it vanreeslab/murphy_base
```

### Continuous integration container
The continuous integration container relies on the base container and adds the needed dependencies (`p4est`, `flups`, etc).
```bash
docker build -f Dockerfile_ci -t vanreeslab/murphy_ci:v1.10 .

# register and put is on DockerHub
docker login
docker push vanreeslab/murphy_ci:v1.10

# to test it
docker run --rm -it vanreeslab/murphy_ci:v1.10
```

### Daily used container
Finally for the beauty of it we add a user to the CI container and install `clangd` (the vscode integrator).

:warning: don't forget to update the version in the first line of the `Dockerfile_daily` accrodingly

```bash
docker build -f Dockerfile_daily -t vanreeslab/murphy:v1.10 .

# register and put is on DockerHub
docker login
docker push vanreeslab/murphy:v1.10

# to test it
docker run --rm -it vanreeslab/murphy:v1.10
```


<!-- 
### 
```
#download flups locally
git clone git@git.immc.ucl.ac.be:examples/flups.git
# build the container
docker build --no-cache -t vanreeslab/murphy:v1.8 .
#push it
docker login
docker push vanreeslab/murphy:v1.8
```
to test the container locally:
```
docker run --rm -it vanreeslab/murphy:v1.8
```
or to permanently build it and mount the murphy folder
```
docker run -it --name murphy -v /Users/tgillis/Dropbox/research/codes/murphy:/murphy vanreeslab/murphy:v1.8
``` -->
