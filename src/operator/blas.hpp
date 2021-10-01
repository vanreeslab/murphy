#ifndef SRC_BLAS_HPP_
#define SRC_BLAS_HPP_

#include "blockoperator.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "doop.hpp"
#include "forestgrid.hpp"
#include "gridblock.hpp"

/**
 * @brief perform the dset operation on a block
 *
 * x = value
 *
 */
class Dset : public BlockOperator {
    real_t value_ = 0.0;

   public:
    explicit Dset() noexcept;
    explicit Dset(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const real_t value, Field* fid_x);
    void ComputeDsetGridBlock(const qid_t* qid, const CartBlock* block, Field* fid_x);
};

/**
 * @brief perform the dcopy operation on a block
 *
 * y = x
 *
 */
class Dcopy : public BlockOperator {
   public:
    explicit Dcopy() noexcept;
    explicit Dcopy(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const Field* fid_x, Field* fid_y);
    void ComputeDcopyGridBlock(const qid_t* qid, CartBlock* block, const Field* fid_x, Field* fid_y);
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
    explicit Daxpy() noexcept;
    explicit Daxpy(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const real_t alpha, const Field* fid_x, const Field* fid_y, Field* fid_z);
    void ComputeDaxpyGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid_x, const Field* fid_y, Field* fid_z);
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
    explicit Dscale() noexcept;
    explicit Dscale(const bidx_t* ghost_len) noexcept;

    void operator()(const ForestGrid* grid, const real_t alpha, Field* fid_x);
    void ComputeDscaleGridBlock(const qid_t* qid, CartBlock* block, Field* fid_x);
};

#endif  // SRC_BLAS_HPP_
