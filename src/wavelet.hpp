#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

template <int order>
class Wavelet : public Interpolator {
   public:
    Wavelet() {
        if (order == 3) {
            hs_[0]  = 0.125;
            hs_[1]  = 1.0;
            hs_[2]  = -0.125;
            sgn_[0] = -1.0;
            sgn_[1] = 1.0;
            sgn_[2] = -1.0;
        } else if (order == 5) {
            hs_[0]  = -3.0 / 128.0;
            hs_[1]  = 11.0 / 64.0;
            hs_[2]  = 1.0;
            hs_[3]  = -11.0 / 64.0;
            hs_[4]  = 3.0 / 128.0;
            sgn_[0] = -1.0;
            sgn_[1] = -1.0;
            sgn_[2] = 1.0;
            sgn_[3] = -1.0;
            sgn_[4] = -1.0;
        } else {
            m_assert(false, "order %d not implemented yet", order);
        }
    }

   protected:
    real_t hs_[order];
    real_t sgn_[order];

    void Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const;
    void Refine_(const interp_ctx_t* ctx) const;
    void Copy_(const interp_ctx_t* ctx) const;
};

#include "wavelet.ipp"

#endif  // SRC_WAVELET_HPP_