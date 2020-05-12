#ifndef SRC_DAXPY_HPP_
#define SRC_DAXPY_HPP_

#include "murphy.hpp"
#include "operator.hpp"

class Daxpy : public OperatorFF2F {
   protected:
    real_t alpha_ = 0.0;

   public:
    explicit Daxpy(real_t alpha);
    void ApplyOpFF2F(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z) override;
};

#endif  // SRC_DAXPY_HPP_