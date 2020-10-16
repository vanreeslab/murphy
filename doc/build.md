# Building MURPHY on your architecture

---------------------
### Dependences
Murphy relies the following external libraries:
- **OpenMPI** or **Intel MPI** as mpi library
- **HDF5** for the I/O operations (version >= 1.10),
- **P4est** for the parallel gestion of the tree (branch `feature-mesh-edge`).

You can also use Google Test library to run the tests (optional).

---------------------
### Compilation
#### Get the dependencies and the `ARCH` file
To compile Murphy you have three different options. Each of them will give you a different `ARCH` file that is used to compile the library.
1. **Local libraries**: install the depencies yourself (you should use a `make_arch/make.YourArchitecture` to indicate the library paths, see `make.docker_gcc` as an example). Then, use the your created file for the `ARCH` file: `ARCH_FILE=make_arch/make.YourArchitecture`.
2. **Docker**: use the `Dockerfile` in the folder `docker/` or you can also download the image from DockerHub (run `docker pull vanreeslab/murphy:v1.7`). Then, the default `ARCH` file (`make_arch/make.docker_gcc`) is valid, no need to change.
3. **VSCode**: You can use the remote container extension to open, build and run the code direction into a Docker container:
    - install the *Remote Container* extension in VSCode,
    - get the DockerHub image (run `docker pull vanreeslab/murphy:v1.7`),
    - open the murphy folder in VSCode and accept to open in a container. For more information about compiling in a remote

#### Docker exemple
We recommend to use the last one as it provides the fastest hand-on approach. However, for advance debugging, the first option is usually the more reliable choice.
To download the image, simply use
```
docker pull vanreeslab/murphy:v1.7
```
To build the container and sync your murphy folder, `MY_MURPHY`, run
```
docker run --name murphy -it -v MY_MURPHY:/murphy vanreeslab/murphy:v1.7
```
You are now within the container and you can access the folder in `/murphy`, which is synced with your local machine.
To exit, simply type `exit` and to relaunch it, enter `docker start -i murphy`.

#### The first build
The first time you compiler, you need to create the `build` directory:
```
mkdir build
```
To check the values (libs, paths,...) that `make` will use, you can do a dry-run:
```
ARCH_FILE=make_arch/make.MyArchFile make info
```
To compile the library, given the `arch` file, use
```
ARCH_FILE=make_arch/make.MyArchFile make -j
```
You will now have an executable `murphy` that you can call.

---------------------
### Building your own ARCH file
By default, the make file uses some default values. If an `ARCH` file is given, we use the values defined in the `ARCH` file.
You can specify the following variables:
- flags
    - `CXX_FLAGS`: the compilation flags you want to use (default `none`)
- FFTW library
    - `FFTW_INC`: the folder with the fftw headers (default `/usr/include`)
    - `FFTW_LIB`: the folder with the fftw libraries (default `/usr/lib`)
    - `FFTW_LIBNAME`: the name of the fftw library on your system (default: `-lfftw3_omp -lfftw3`)
- HDF5 library
    - `HDF5_INC`: the folder with the hdf5 headers (default `/usr/include`)
    - `HDF5_LIB`: the folder with the hdf5 libraries (default `/usr/lib`)
    - `HDF5_LIBNAME`: the name of the hdf5 library on your system (default: `-lhdf5`)
- P4est library
    - `P4EST_INC`: the folder with the p4est headers (default `/usr/include`)
    - `P4EST_LIB`: the folder with the p4est libraries (default `/usr/lib`)
    - `P4EST_LIBNAME`: the name of the p4est library on your system (default: `-lsc -lp4est`)
- Google Test library
    - `GTEST_INC`: the folder with the google test headers (default `/usr/include`)
    - `GTEST_LIB`: the folder with the google test libraries (default `/usr/lib`)
    - `GTEST_LIBNAME`: the name of the google test library on your system (default: `-lgtest`)

---------------------
### Possible compilation flag
Some compilations flags are available to change the behavior of the code:
- `-DBLOCK_GS=X` sets the number of ghost points to use (default `X=2`)
- `-DWAVELET_N=X` dictates the interpolation order of the wavelets (default `X=2`)
- `-DWAVELET_NT=X` dictates the moment order of the wavelets (default `X=2`)
- `-DLOG_ALLRANKS` will enable log on every processor. By default, only the master logs
- `-DVERBOSE` enable extended logs
- `-DNDEBUG` disable the assertion checks and the other debuging sections
- `-DLOG_MUTE` disable every logs
<!-- - ```-DMG_GAUSSSEIDEL``` uses the gauss-seidel smoother instead of the Jacobi one -->

To use them, you can append the make command, e.g. to change the wavelet behavior
```sh
make OPTS="-DDBLOCK_GS=4 -DWAVELET_N=4"
```
or add it to the `ARCH` file:
```makefile
override OPTS += -DDBLOCK_GS=4 -DWAVELET_N=2
```

In doubts, you should always run `make info` to test the different options.

:warning: the order of the wavelet must be compatible with the number of ghost points!

---------------------
### Testing
See the details about automatic testing in [this page](doc/contribute.md).

To run the tests, you need to recompile all the sources (it's important to delete the old ones if you want to retrieve)
```
make destroy
ARCH_FILE=make_arch/make.docker_gcc_test make test -j
```
Then you can lunch the tests
```
./murphy_test
```

---------------------
### Documentation
To get the doxygen documentation, go to the main folder and enter
```
doxygen doc/Doxyfile
```
You will need `dot` to get a visual graphs (see `HAVE_DOT` option in the Doxyfile).
On MacOS, you can install it using homebrew: `brew install graphviz`.

<!-- 
-----------------------
### Installation from scratch
We here assume that you have nothing installed and we showcase how to build OpenMPI, HDF5 and p4est all together. 

:warning: we here perform a debug build, do not use this one for production use
```sh
# Open MPI
wget https://download.open-mpi.org/release/open-mpi/v4.0/openmpi-4.0.4.tar.gz
tar -xvf openmpi-4.0.4.tar.gz 
cd openmpi-4.0.4
./configure --prefix=${HOME}/ompi_4.0.4 --enable-debug --enable-memchecker
make install -j12
``` -->

<!-- *to be continued* -->