# Build instructions for *murphy*


### Dependences
Murphy relies the following external libraries:
- **OpenMPI** or **Intel MPI** as mpi library
- **HDF5** for the I/O operations (version >= 1.10),
- **P4est** for the parallel gestion of the tree (branch `feature-mesh-edge`).

### Compilation
To compile Murphy, you have three different options:
1. **Local libraries**: install the depencies yourself (you should use a `make_arch/make.YourArchitecture` to indicate the library path, see `make.docker_gcc` as an example),
2. **Docker**: use the `Dockerfile` in the folder `docker/` or you can also download the image from DockerHub (run `docker pull immc/murphy-ci:v1.1`),
3. **VSCode**: You can use the remote container extension to open, build and run the code direction into a Docker container:
    - install the *Remote Container* extension in VSCode,
    - get the DockerHub image (run `docker pull immc/murphy-ci:v1.1`),
    - open the murphy folder in VSCode and accept to open in a container. For more information about compiling in a remote

We recommend to use the last one as it provides the fastest hand-on approach.


### Possible compilation flag
You can modify the login level by adding the following compilations flags in your `make_arch` file:
- ```-DLOG_ALLRANKS``` will enable log on every processor. By default, only the master logs
- ```-DVERBOSE``` enable extended logs
- ```-DNDEBUG``` disable the assertion checks and the debug comments