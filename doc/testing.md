# Testing in Murphy

Tests in MURPHY use the Google Test framework, which organizes individual tests into a series of test suites. By convention, all MURPHY test suite names (and instantiation names for parametrized tests) are `CamelCase`, while individual tests are named in `snake_case`. This page provides descriptions of each test and is organized according to the Google Test hierarchy. The source file containing any given test is listed in its description. 

## Test Hierarchy 
A high level overview of all test can be obtained by running `./murphy_test --gtest_list_tests`. A compressed version of the output is below, with the format

```
(InstiationName/)TestSuiteName.
    test_1_name(0-N)
    test_2_name
```
Here the `InstantiationName` prefix is present only for parametrized tests, and `(0-N)` represents a parametrized test with `N+1` parameters. 

#### Currently implemented tests (updated October 19th, 2021):
```
BlockTypeDeathTest.
  doop

ValidStencilUniform.
  weno_extrap_cosinus

ValidRK.
  rk3_tvd

ValidWaveletKernel.
  filter_length

ValidPartitioner/TestPartitioner.
  nblocks_leq_ncpus/0 
  nblocks_geq_ncpus/0 

ValidGhost/ValidWaveletInterpolation.
  ghost_reconstruction_periodic_sin/0-6
  ghost_reconstruction_periodic_cos/0-6
  ghost_reconstruction_extrap_cos/0-6
  ghost_reconstruction_perper_dirichlet0_polynom/0-6
  ghost_reconstruction_perper_neuman0_cos/0-6

ValidWavelet/TwoLevel.
  periodic/0-255
  flipfop/0-255

ValidStencil/Adapt.
  weno_periodic_cosinus/0-255

ValidWavelet/Epsilon.
  periodic/0-2
  extrap/0-2
```

## Test Descriptions

#### `BlockTypeDeathTest.doop`
From `block_types.cpp`. Asserts that the `DoOp` interface works properly with the different `BlockDataType` and throws an error when the interface is misused. 

#### `ValidStencilUniform.weno_extrap_cosinus`
From `valid_stencil.cpp`. Chooses a random uniform velocity field, sets up a sinusoidal scalar field with two periods in each direction, and computes an advection RHS using four different stencils (WENO 3, WENO 5, conservative 3, conservative 5). For the WENO cases, asserts that the stencil is conservative up to floating point error. For all cases, runs at two different spatial resolutions, and asserts that the stencils achieve their expected convergence order.

#### `ValidRK.rk3_tvd` 
From `valid_rk.cpp`. Sets up an advection test with uniform velocity and a Gaussian scalar field. RHS is computed analytically, but the time integration is done with the RK3 implementation. The test does this computation using several different time steps and asserts that the observed convergence order is better than 2.9.

#### `ValidWaveletKernel.filter_length`
From `valid_wavelet_kernel.cpp`. Asserts that members of `InterpolatingWavelet` related to ghosting size have the correct values.

#### `ValidPartitioner/TestPartitioner.nblocks_leq_ncpus/0`
From `valid_partitioner.cpp`. Sets up a grid with 16 blocks, adapts (and partitions) it to reduce the number of blocks to two. With the current tests configuration, it change the grid from a  number of blocks by CPUs > 1 to a number of blocks per CPUs <= 1. Checks that there is no _nan_ or _segfault_. Adapts (and partitions) the grid and retrieves its initial configuration. Performs the check again.

#### `ValidPartitioner/TestPartitioner.nblocks_geq_ncpus/0`
From `valid_partitioner.cpp`. Sets up a grid with 16 blocks, adapts (and partitions) it to increase the number of blocks to 128. Checks that there is no _nan_ or _segfault_. Adapts (and partitions) the grid and retrieves its initial configuration. Performs the check again.

#### `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_periodic_sin/0-6`
From `valid_wavelet_convergence.cpp`. Sets up a periodic two-level grid with a sinusoidal scalar field, executes a ghost pull, and measures the error in the field (including ghost values). The calculation is done at two different spatial resolutions, and the test asserts that the observed convergence rate matches the interpolation order of the wavelet to within a prescribed tolerance. For wavelet N.0, the test also asserts that there is no error in the coarse-level ghosts. The test parameter represents the number of ghosts in each dimension.

#### `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_periodic_cos/0-6`
From `valid_wavelet_convergence.cpp`. Same as `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_periodic_sin/0-6`, but the scalar field is a product of cosines.

#### `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_extrap_cos/0-6`
From `valid_wavelet_convergence.cpp`. Same as `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_periodic_co/0-6`, but the grid is set up with an extrapolation boundary condition in all directions, and the assertion relating to wavelets N.0 is dropped.

#### `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_perper_dirichlet0_polynom/0-6`
From `valid_wavelet_convergence.cpp`. Sets up a two-level grid that is periodic in two directions and Dirichlet in the third, evaluates a polynomial that respects the Dirichlet boundary condition, executes a ghost pull, and measures the error in the field (including ghost values). The calculation is done at two different spatial resolutions, and the test asserts that the observed convergence rate matches the interpolation order of the wavelet to within a prescribed tolerance. This is repeated for all three grid configurations that have a Dirichlet condition in only one direction. The test parameter represents the number of ghosts in each dimension.

#### `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_perper_neuman0_cos/0-6`
From `valid_wavelet_convergence.cpp`. Same as `ValidGhost/ValidWaveletInterpolation.ghost_reconstruction_perper_dirichlet0_polynom/0-6`, but with a Neumann boundary condition.

#### `ValidWavelet/TwoLevel.periodic/0-255`
From `valid_wavelet_twolevel.cpp`. Sets up a Gaussian scalar field on a grid with periodic boundary conditions, forces one level of coarsening on a set of blocks, and then refines back to the original level. The test asserts that the resulting error is less than the constant factor `InterpolatingWavelet::eps_const()` times the largest detail coefficient on the original grid. For wavelets N.2, the test also asserts that the 0th and 1st moments are conserved up to floating point error. The parameter represents all 256 possible subsets of eight octants.

#### `ValidWavelet/TwoLevel.flipfop/0-255`
From `valid_wavelet_twolevel.cpp`. Same as `ValidWavelet/TwoLevel.periodic/0-255`, but with a flip-flop scalar field. For wavelets N.2, the test asserts that the 0th moments (but not the 1st) is conserved to floating point precision.

#### `ValidStencil/Adapt.weno_periodic_cosinus/0-255`
From `valid_stencil.cpp`. Same setup as `ValidStencilUniform.weno_extrap_cosinus`, but performed on an adapted grid. The parameter represents all 256 possible configurations of eight octants on two resolution levels. Runs at two spatial resolutions, and asserts that the stencils achieve their expected convergence order.

#### `ValidWavelet/Epsilon.periodic/0-2`
From `valid_wavelet_epsilon.cpp`. Sets up a Gaussian scalar field on a grid with periodic boundary conditions, then repeatedly calls `Grid::Coarsen()` with a preset coarsening criterion (epsilon) until the grid reaches a desired coarse resolution level. Once this is done, the grid is refined to its original level and the test asserts that the resulting error is less than three times epsilon. For wavelets N.2, the test also asserts that the 0th and 1st moments are conserved up to floating point error. The test parameter represents different values of epsilon.

#### `ValidWavelet/Epsilon.extrap/0-2`
From `valid_wavelet_epsilon.cpp`. Same setup as `ValidStencilUniform.weno_extrap_cosinus`, but on a grid with extrapolation boundary conditions. Becuase of the extrapolation, this test does not assert that moments are conserved.