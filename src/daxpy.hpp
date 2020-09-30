#ifndef SRC_DAXPY_HPP_
#define SRC_DAXPY_HPP_

#include "blockoperator.hpp"
#include "doop.hpp"
#include "gridblock.hpp"
#include "murphy.hpp"

/**
 * @brief perform the daxpy operation on a block
 * 
 * z = alpha * x + y
 * 
 */
class Daxpy : BlockOperator {
   protected:
    real_t alpha_ = 0.0;

   public:
    explicit Daxpy(real_t alpha);

    void ComputeDaxpyGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z);
};

#endif  // SRC_DAXPY_HPP_