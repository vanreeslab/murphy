#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

template <int order>
class Wavelet : public Interpolator {
   protected:
    real_t hs_[order];
    real_t sgn_hs_[order];

    real_t ga_[order];
    real_t sgn_ga_[order];

   public:
    Wavelet() {
        if (order == 3) {
            // h synthesis - refinement
            hs_[0] = 0.125;
            hs_[1] = 1.0;
            hs_[2] = -0.125;
            // sign convention
            sgn_hs_[0] = -1.0;
            sgn_hs_[1] = 1.0;
            sgn_hs_[2] = -1.0;
            // g analysis - analysis
            ga_[0] = -0.0625;
            ga_[1] = 0.5;
            ga_[2] = 0.0625;
            // sign convention
            sgn_ga_[0] = 1.0;
            sgn_ga_[1] = -1.0;
            sgn_ga_[2] = 1.0;
        } else if (order == 5) {
            // h synthesis
            hs_[0] = -3.0 / 128.0;
            hs_[1] = 11.0 / 64.0;
            hs_[2] = 1.0;
            hs_[3] = -11.0 / 64.0;
            hs_[4] = 3.0 / 128.0;
            // sign convention
            sgn_hs_[0] = -1.0;
            sgn_hs_[1] = -1.0;
            sgn_hs_[2] = 1.0;
            sgn_hs_[3] = -1.0;
            sgn_hs_[4] = -1.0;
            // g analysis
            ga_[0] = 0.01171875;
            ga_[1] = -0.171875;
            ga_[2] = 0.5;
            ga_[3] = 0.171875;
            ga_[4] = -0.01171875;
            // sign convention
            sgn_ga_[0] = 1.0;
            sgn_ga_[1] = 1.0;
            sgn_ga_[2] = -1.0;
            sgn_ga_[3] = 1.0;
            sgn_ga_[4] = 1.0;
        } else {
            m_assert(false, "order %d not implemented yet", order);
        }
    }

    void Criterion(MemLayout* block, real_p data, real_t* criterion) override;

   protected:
    void Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const override;
    void Refine_(const interp_ctx_t* ctx) const override;
    void Copy_(const interp_ctx_t* ctx) const override;
    void Detail_(const interp_ctx_t* ctx, real_t* details_inf_norm) const;
};

#endif  // SRC_WAVELET_HPP_

#include "wavelet.ipp"