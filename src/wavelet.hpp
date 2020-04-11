#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

template <int order>
class Wavelet : public Interpolator {
   protected:
    real_t hs_[order];
    real_t sgn_[order];

    void Coarsen_(const lid_t dlvl, const real_p data_src, real_p data_trg) const;
    void Refine_(const real_p data_src, real_p data_trg) const;
    void Copy_(const real_p data_src, real_p data_trg) const;
};

template<>
class Wavelet<3> {
   protected:
    real_t hs_[3]  = {0.125, 1.0, -0.125};
    real_t sgn_[3] = {-1.0, 1.0, -1.0};
};

template<>
class Wavelet<5> {
   protected:
    real_t hs_[5]  = {-3.0 / 128.0, 11.0 / 64.0, 1.0, -11.0 / 64.0, 3.0 / 128.0};
    real_t sgn_[5] = {-1.0,-1.0, 1.0, -1.0,-1.0};
};


#include "wavelet.ipp"

#endif  // SRC_WAVELET_HPP_