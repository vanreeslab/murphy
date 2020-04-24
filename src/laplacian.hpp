#ifndef SRC_LAPLACIAN_HPP
#define SRC_LAPLACIAN_HPP

#include "murphy.hpp"
#include "operator.hpp"
#include "stencil.hpp"

template <sid_t length>
class LaplacianCross : public Stencil {
    sid_t  ida_ = 0;       //!< current working dimension
    real_t coef_[length];  //!< coefficients of the stencils to apply
   public:
    LaplacianCross(Grid* grid) : Stencil(grid) {
        if (length == 3) {
            coef_[0] = +1.0;
            coef_[1] = -2.0;
            coef_[2] = +1.0;
        } else {
            m_assert(false, "not coded yet");
        }
    };

   protected:
    void ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) override;
    void ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_src, Field* fid_trg) override;
};

#endif  // SRC_LAPLACIAN_HPP
#include "laplacian.ipp"