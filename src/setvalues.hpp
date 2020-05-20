/**
 * @file setvalues.hpp
 * @brief defines several Operators to set various functions in the blocks
 * @version
 */
#ifndef SRC_SETVALUES_HPP_
#define SRC_SETVALUES_HPP_

#include "murphy.hpp"
#include "operator.hpp"

/**
 * @brief sets a gaussian in every dimension of a field
 */
class SetGaussian : public OperatorF {
   protected:
    real_t sigma_     = 0.1;
    real_t center_[3] = {0, 0, 0};

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetGaussian(real_t sigma, real_t center[3]);
};

/**
 * @brief sets a the absolute value function in every dimension of a field
 */
class SetAbs : public OperatorF {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3]  = {0.0, 0.0, 0.0};

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetAbs(real_t alpha[3], real_t center[3]);
};

/**
 * @brief sets a discontinuity in every dimension of a field
 */
class SetJump : public OperatorF {
   protected:
    real_t center_[3] = {0, 0, 0};
    real_t alpha_[3]  = {0.0, 0.0, 0.0};

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetJump(real_t alpha[3], real_t center[3]);
};

/**
 * @brief sets a sinus in every dimension of a field
 */
class SetSinus : public OperatorF {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetSinus(real_t length[3], real_t freq[3]);
};

/**
 * @brief sets a cosinus in every dimension of a field
 */
class SetCosinus : public OperatorF {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetCosinus(real_t length[3], real_t freq[3]);
};

class SetLaplaCosinus : public OperatorF {
   protected:
    real_t freq_[3]   = {0, 0, 0};
    real_t length_[3] = {0.0, 0.0, 0.0};

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetLaplaCosinus(real_t length[3], real_t freq[3]);
};

/**
 * @brief sets a polynomial in every dimension of a field
 */
class SetPolynom : public OperatorF {
   protected:
    lid_t  deg_[3] = {0, 0, 0};        //!< the degree of the polynomial
    real_t dir_[3] = {0.0, 0.0, 0.0};  //!< the direction concerned: 1.0 means involved, 0.0 means not involved

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetPolynom(lid_t degree[3], real_t direction[3]);
};


class SetExponential : public OperatorF {
   protected:
    real_t center_[3] = {0.0, 0.0, 0.0};
    real_t sigma_[3]  = {0.0, 0.0, 0.0};
    real_t alpha_ = 1.0;

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetExponential(real_t center[3], real_t sigma[3],real_t alpha);
};

class SetErf : public OperatorF {
   protected:
    real_t center_[3] = {0.0, 0.0, 0.0};
    real_t sigma_[3]  = {0.0, 0.0, 0.0};
    real_t alpha_ = 1.0;

    void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

   public:
    SetErf(real_t center[3], real_t sigma[3],real_t alpha);
};

// /**
//  * @brief sets a cosinus multiplied by an exponential in every dimension of a field
//  */
// class SetExpoCosinus : public OperatorF {
//    protected:
//     real_t center_[3] = {0, 0, 0};
//     real_t sigma_[3]  = {0, 0, 0};
//     real_t freq_[3]   = {0, 0, 0};
//     real_t length_[3] = {0.0, 0.0, 0.0};

//     void ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override;

//    public:
//     SetExpoCosinus(real_t center[3], real_t sigma[3], real_t length[3], real_t freq[3]);
// };

#endif  // SRC_SETVALUES_HPP_
