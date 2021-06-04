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
    explicit Dset(m_ptr<const Wavelet> interp);

    void operator()(m_ptr<const ForestGrid> grid, const real_t value, m_ptr<Field> fid_x);
    void ComputeDsetGridBlock(m_ptr<const qid_t> qid, m_ptr<const GridBlock> block, m_ptr<Field> fid_x);
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
    explicit Dcopy(m_ptr<const Wavelet> interp);

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, m_ptr<Field> fid_y);
    void ComputeDcopyGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x, m_ptr<Field> fid_y);
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
    explicit Daxpy(m_ptr<const Wavelet> interp);

    void operator()(m_ptr<const ForestGrid> grid, const real_t alpha, m_ptr<const Field> fid_x, m_ptr<const Field> fid_y, m_ptr<Field> fid_z);
    void ComputeDaxpyGridBlock(m_ptr<const qid_t> qid, m_ptr<const GridBlock> block, m_ptr<const Field> fid_x, m_ptr<const Field> fid_y, m_ptr<Field> fid_z) ;
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
    explicit Dscale(m_ptr<const Wavelet> interp);

    void operator()(m_ptr<const ForestGrid> grid, const real_t alpha, m_ptr<Field> fid_x);
    void ComputeDscaleGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid_x);
};

#endif  // SRC_BLAS_HPP_
