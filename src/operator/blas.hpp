#ifndef SRC_BLAS_HPP_
#define SRC_BLAS_HPP_

#include "blockoperator.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "doop.hpp"
#include "gridblock.hpp"
#include "forestgrid.hpp"

/**
 * @brief perform the dset operation on a block
 *
 * x = value
 *
 */
class Dset : public BlockOperator {
    real_t value_ = 0.0;

   public:
    explicit Dset();
    explicit Dset(const Wavelet*  interp);

    void operator()(const ForestGrid*  grid, const real_t value, Field*  fid_x);
    void ComputeDsetGridBlock(const qid_t*  qid, const GridBlock*  block, Field*  fid_x);
};

/**
 * @brief perform the dcopy operation on a block
 *
 * y = x
 *
 */
class Dcopy : public BlockOperator {
   public:
    explicit Dcopy();
    explicit Dcopy(const Wavelet*  interp);

    void operator()(const ForestGrid*  grid, const Field*  fid_x, Field*  fid_y);
    void ComputeDcopyGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid_x, Field*  fid_y);
};

/**
 * @brief perform the daxpy operation on a block
 *
 * z = alpha * x + y
 *
 */
class Daxpy : public BlockOperator {
   protected:
    real_t alpha_ = 0.0;

   public:
    explicit Daxpy();
    explicit Daxpy(const Wavelet*  interp);

    void operator()(const ForestGrid*  grid, const real_t alpha, const Field*  fid_x, const Field*  fid_y, Field*  fid_z);
    void ComputeDaxpyGridBlock(const qid_t*  qid, const GridBlock*  block, const Field*  fid_x, const Field*  fid_y, Field*  fid_z) ;
};

/**
 * @brief perform the scale operation on a block
 *
 * x = alpha * x
 *
 */
class Dscale : public BlockOperator {
   protected:
    real_t alpha_ = 0.0;

   public:
    explicit Dscale();
    explicit Dscale(const Wavelet*  interp);

    void operator()(const ForestGrid*  grid, const real_t alpha, Field*  fid_x);
    void ComputeDscaleGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid_x);
};

#endif  // SRC_BLAS_HPP_
