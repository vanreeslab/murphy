#ifndef SRC_LAPLACIAN_HPP_
#define SRC_LAPLACIAN_HPP_

#include "murphy.hpp"
#include "operator.hpp"
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
    real_t coef_[length];  //!< coefficients of the laplacian to apply
   public:
    LaplacianCross(Grid* grid) : Stencil(grid) {
        if (length == 3) {
            coef_[0] = +1.0;
            coef_[1] = -2.0;
            coef_[2] = +1.0;
        } else {
            m_assert(false, "not coded yet");
        }
    }

   protected:
    void ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) override;
    void ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) override;
};

#endif  // SRC_LAPLACIAN_HPP_

#include "laplacian.ipp"
