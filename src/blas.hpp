#ifndef SRC_BLAS_HPP_
#define SRC_BLAS_HPP_

#include "blockoperator.hpp"
#include "doop.hpp"
#include "gridblock.hpp"
#include "murphy.hpp"

/**
 * @brief perform the dcopy operation on a block
 * 
 * y = x
 * 
 */
class Dcopy : public BlockOperator {
   public:
    explicit Dcopy();
    explicit Dcopy(const Wavelet* interp);
    void ComputeDcopyGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y);
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
    explicit Daxpy(real_t alpha);
    void ComputeDaxpyGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z);
};

/**
 * @brief perform the scale operation on a block
 * 
 * x = alpha * x
 * 
 */
class Scale : public BlockOperator {
   protected:
    real_t alpha_ = 0.0;

   public:
    explicit Scale(real_t alpha);
    void ComputeScaleGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x);
};

#endif  // SRC_BLAS_HPP_