#ifndef SRC_SETVALUE_HPP_
#define SRC_SETVALUE_HPP_

#include "murphy.hpp"
#include "operator.hpp"

class SetGaussian : public OperatorF {
   protected:
    real_t sigma_     = 0.1;
    real_t center_[3] = {0, 0, 0};

    void apply(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetGaussian(real_t sigma, real_t center[3]);
};

// class SetValueSinus : public Operator {
//    public:

// }

#endif  // SRC_SETVALUE_HPP_
