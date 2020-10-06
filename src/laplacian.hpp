#ifndef SRC_LAPLACIAN_HPP_
#define SRC_LAPLACIAN_HPP_

#include "murphy.hpp"
#include "stencil.hpp"

/**
 * @brief Implements a stencil laplacian in cross of order (length+1)/2
 * 
 * e.g. if the length is 3, the order is 2
 * 
 * @tparam length 
 */
template <sid_t length>
class LaplacianCross : public Stencil {
   protected:
    real_t coef_[length];  //!< coefficients of the laplacian to apply

   public:
    explicit LaplacianCross() : Stencil() {
        // get the stencil
        if (length == 3) {
            coef_[0] = +1.0;
            coef_[1] = -2.0;
            coef_[2] = +1.0;
        } else if (length == 5) {
            coef_[0] = -1.0 / 12.0;
            coef_[1] = +4.0 / 3.0;
            coef_[2] = -5.0 / 2.0;
            coef_[3] = +4.0 / 3.0;
            coef_[4] = -1.0 / 12.0;
        } else {
            m_assert(false, "not coded yet");
        }
    }

    virtual void ApplyStencilInner(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) override;
    virtual void ApplyStencilOuter(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) override;
};

#endif  // SRC_LAPLACIAN_HPP_

#include "laplacian.inc"
