#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

class Wavelet : public Interpolator {
   protected:
    lid_t  hslen_  = 1;
    real_t hs_[1]  = {0.0};
    real_t sgn_[1] = {0.0};

    void Coarsen_(const lid_t dlvl, const real_p data_src, real_p data_trg) const;
    void Refine_(const real_p data_src, real_p data_trg) const;
    void Copy_(const real_p data_src, real_p data_trg) const;
};

class Wavelet3 : public Wavelet {
   protected:
    lid_t  hslen_  = 3;
    real_t hs_[3]  = {0.125, 1.0, -0.125};
    real_t sgn_[3] = {-1.0, 1.0, -1.0};
};

#endif  // SRC_WAVELET_HPP_