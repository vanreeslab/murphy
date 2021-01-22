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

/**
 * @brief perform the max(fabs()) operation on a block, i.e. return the infinite norm of a field
 *
 * when the values are asked back, do a AllReduce MPI call
 *
 */
class Dmax : public BlockOperator {
   protected:
    real_t max_;

   public:
    explicit Dmax();
    explicit Dmax(m_ptr<const Wavelet> interp);

    real_t operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x);
    void   ComputeDmaxGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
};

/**
 * @brief perform the min and max operation on a block and return the two values after and AllReduce
 *
 */
class Dminmax : public BlockOperator {
   protected:
    lda_t ida_;
    real_t max_, min_;

   public:
    explicit Dminmax();
    explicit Dminmax(m_ptr<const Wavelet> interp);

    void operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, real_t* min, real_t* max);
    void ComputeDminmaxGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x);
};

#endif  // SRC_BLAS_HPP_
