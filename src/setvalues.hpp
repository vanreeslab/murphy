#ifndef SRC_SETVALUE_HPP_
#define SRC_SETVALUE_HPP_

#include "murphy.hpp"
#include "operator.hpp"

class SetGaussian : public OperatorF {
   protected:
    real_t sigma_     = 0.1;
    real_t center_[3] = {0, 0, 0};

    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetGaussian(real_t sigma, real_t center[3]);
};


class SetAbs : public OperatorF {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3] = {0.0, 0.0, 0.0};

    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetAbs(real_t alpha[3], real_t center[3]);
};



class SetJump : public OperatorF {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3] = {0.0, 0.0, 0.0};

    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetJump(real_t alpha[3], real_t center[3]);
};



class SetSinus : public OperatorF {
   protected:
    real_t freq_[3] = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetSinus(real_t length[3], real_t freq[3]);
};


class SetCosinus : public OperatorF {
   protected:
    real_t freq_[3] = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void ApplyOperatorF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetCosinus(real_t length[3], real_t freq[3]);
};



#endif  // SRC_SETVALUE_HPP_
