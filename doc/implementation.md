x# Implementation choices

Here we go through the implementation strategy that guided us to the current version.

### 1. Grid initialization
```
lid_t  L[3]      = {2, 1, 3};
bool periodic[3] = {false, false, false};
Grid* grid       = new Grid(1, periodic, L, MPI_COMM_WORLD, NULL);
```
We create a distributed (MPI+OpenMP) forest of octrees, organized in a rectangular shape. The number of tree in each direction is given by `L[0]`x`L[1]`x`L[2]` (each octree is a unit cube of physical size `1`x`1`x`1`);
The tree management operations are performed by the library `p4est`. Additionally to the forest size, the periodicity of its boundary is given to p4est and cannot change during the simulation. Other boundary conditions are given afterwards. Each leaf of the tree contains a `GridBlock`, a continuous memory space of fixed size. The block size in direction is given by the number of unknowns in its core, `M_N`, and the number of ghost points, `M_GS`. The memory stride is then given by `M_STRIDE = 2 * M_GS + M_N`. (Those contants are fixed in the file `defs.hpp`.)


### 2. Field creation
```
Field* vort = new Field("vorticity", 3);
```
One defines a field by a unique `name`, `vorticity`, and a dimensionality `lda`, `3`. The field does not contain its own memory and is almost an empty shell. However, its `name` and `lda` are used each time one wants to perform an opeation on it.


### 3. Grid-Field association
```
grid->AddField(vort);
```
The association of a field to a grid will perform the real memory allocation. A field can be associated to multiple grids and a grid contains a map of all the fields it is associated with.
Each field added will trigger a loop on the blocks to initialize the memory. The memory is allocated continuously, each dimension separated by `M_STRIDE * M_STRIDE * M_STRIDE`. This means that, in memory, we have `[Field(ida=0) Field(ida=1) Field(ida=2)]` in one continuous array.


### 4. Boundary conditions
```
vort->bctype(M_BC_ODD);
```
The imposition of a boudnary condition is done on the field. Each dimension, `ida`, will get 6 boundary conditions: `x`-, `x`+, `y`-, `y`+, `z`-, `z`+. The boundary conditions will be used any time ghost points are needed for the field. If a direction is periodic, the boundary condition will be discarded automatically. 


### 5. Ghost points
```
grid->GhostPull(vort);
```
The ghost points computation is done using the wavelets to reconstruct missing information in the case of a level mismatch. The ghost reconstruction is done dimension by dimension, allowing to overlap wavelet reconstruction with the communication for the next dimension. Each field owns a boolean to indicate if the ghost points are up-to-date. Hence, calling the `GhostPull` with an already up-to-date field will return immediatly. To ensure consistency of this boolean, please use carefully the operators abstract classes.


### 6. Operators
To perfom operations on the blocks, we inherite from abstact classes. There are currently the following different operators
<!-- 1. `OperatorS`: a simple operator, does not interact with any field -->
1. `OperatorF`: operates on a given field, ghost status is automatically changed to `false` afterwards.
1. `ConstOperatorF`: act on a given field, without changing its content. The ghost status is unchanged afterwards.
1. `OperatorF2F`: operates on two fields, a source and a target. The ghost status of the target is automatically changed to `false` afterwards
1. `OperatorFF2F`: operates on three fields, two sources and one target. The ghost status of the target is automatically changed to `false` afterwards
1. `ConstOperatorFF`: operates on two fields, without changing their content. No ghost status change afterwards.

E.g. the error computation between two fields is a `ConstOperatorFF` operator, etc. It is possible to inherit from multiple at the same time as they don't have a common ancestor.

All operators overload the corresponding function `ApplyMyOperator()` which is triggered by the executation of the `operator()` function on the entire field.
Given the operator type, the arguments of the `ApplyMyOperator` function changes. In any case, the `ApplyMyOperator` function, as a member function, has access to the content of the current object.

:warning: The `ApplyMyOperator` function is automatically processed in a multi-threaded section (using OpenMP). This means that the threads will execute the functions on a different blocks __at the same time__, each of them with a copy of the adress of the current object (`this` pointer). Hence, any operation performed in this function has to be __thread-safe__! (use `#pragma omp critical`, `#pragma omp single`, `#pragma omp atomic`,... if needed).


### 7. Stencils
The computation of a stencil is done using the `Stencil` class. This class implements the overlap between the stencil computation and the required ghost exchange.
To improve the computation/communication overlap, the ghost exchange is done dimension by dimension. This means that we follow the following approach, for each dimension:

1. we start the ghost exchange in the dimension `ida_` of the source field
2. we compute the inner part of the stencil that depends on the dimension `ida_` of the source field
3. we receive the ghost for the dimension `ida_` of the source field
4. we compute the outer part of the stencil, that depends on the ghost just recevied.

Additionnaly to this routine, we also intertwine the other dimension's send/receive MPI calls.

Hence, the innner and outer application of the stencil **must** we written with respect to the current available dimension of the source field.
It is possible to write any stencil like that, even the cross-products, where the other dimensions can be accessed with `(ida_+1)%3` and `(ida_+2)%3`.

:warning: As the `Stencil` class is an operator, the inner and outer functions are automatically processed in a multi-threaded section (using OpenMP).

