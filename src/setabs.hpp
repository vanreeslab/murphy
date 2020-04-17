#ifndef SRC_SETABS_HPP_
#define SRC_SETABS_HPP_

#include "murphy.hpp"
#include "operator.hpp"

class SetAbs : public OperatorF {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3] = {0.0, 0.0, 0.0};

    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetAbs(real_t sigma, real_t center[3]);
};

#endif  // SRC_SETABS_HPP_
