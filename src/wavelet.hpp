#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

template <int order>
class Wavelet : public Interpolator {
   public:
    Wavelet() {
        if (order == 3) {
            real_t hs_[3]  = {0.125, 1.0, -0.125};
            real_t sgn_[3] = {-1.0, 1.0, -1.0};
        } else if (order == 5) {
            real_t hs_[5]  = {-3.0 / 128.0, 11.0 / 64.0, 1.0, -11.0 / 64.0, 3.0 / 128.0};
            real_t sgn_[5] = {-1.0, -1.0, 1.0, -1.0, -1.0};
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