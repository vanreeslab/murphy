#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

// check if the compilation defines the order of the wavelet. if not, we do it
#ifndef WAVELET_N
#define M_WAVELET_N 2
#else
#define M_WAVELET_N WAVELET_N
#endif

#ifndef WAVELET_NT
#define M_WAVELET_NT 2
#else
#define M_WAVELET_NT WAVELET_NT
#endif

/**
 * @brief Implement the different wavelet types
 * The wavelet application considered here is based on the lifting/dual lifting scheme.
 * More details on the wavelet approach can be found in the wavelet tutorial (https://github.com/van-Rees-Lab/wavelet_tutorial).
 * 
 * For any couple N,Nt (with Nt <= N) the filter coefficients are computed by (https://github.com/van-Rees-Lab/wavelet_tutorial/blob/master/code/wavelet_xy.py)
 * 
 * @tparam N 
 * @tparam NT 
 */
template <lda_t N = M_WAVELET_N, lda_t NT = M_WAVELET_NT>
class Wavelet : public InterpolatingWavelet {
   public:
    const sid_t N() const { return 0; }
    const sid_t Nt() const { return 0; }

    // filter
    const sid_t   len_ha() const override { return 0; };
    const sid_t   len_ga() const override { return 0; };
    const sid_t   len_gs() const override { return 0; };
    const real_t* ha() const override { return 0; };
    const real_t* ga() const override { return 0; };
    const real_t* gs() const override { return 0; };
};

template <>
class Wavelet<2, 2> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 5;
    const lid_t  len_ga_ = 3;
    const lid_t  len_gs_ = 2;
    const real_t ha_[5]  = {-0.125, 0.25, 0.75, 0.25, -0.125};
    const real_t ga_[3]  = {-0.5, 1.0, -0.5};
    const real_t gs_[2]  = {0.5, 0.5};

    const sid_t N() const { return 2; }
    const sid_t Nt() const { return 2; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};
template <>
class Wavelet<4, 0> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 1;
    const lid_t  len_ga_ = 7;
    const lid_t  len_gs_ = 4;
    const real_t ha_[1]  = {1.0};
    const real_t ga_[7]  = {0.0625, 0.0, -0.5625, 1.0, -0.5625, 0.0, 0.0625};
    const real_t gs_[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};

    const sid_t N() const { return 4; }
    const sid_t Nt() const { return 0; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};
template <>
class Wavelet<4, 2> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 9;
    const lid_t  len_ga_ = 7;
    const lid_t  len_gs_ = 4;
    const real_t ha_[9]  = {0.015625, 0.0, -0.125, 0.25, 0.71875, 0.25, -0.125, 0.0, 0.015625};
    const real_t ga_[7]  = {0.0625, 0.0, -0.5625, 1.0, -0.5625, 0.0, 0.0625};
    const real_t gs_[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};

    const sid_t N() const { return 4; }
    const sid_t Nt() const { return 2; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};

// template <> Wavelet<4,4>::len_ha_=13;

template <>
class Wavelet<4, 4> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 13;
    const lid_t  len_ga_ = 7;
    const lid_t  len_gs_ = 4;
    const real_t ha_[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
    const real_t ga_[7]  = {0.0625, 0.0, -0.5625, 1.0, -0.5625, 0.0, 0.0625};
    const real_t gs_[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};

    const sid_t N() const { return 4; }
    const sid_t Nt() const { return 4; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};

template <>
class Wavelet<6, 0> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 1;
    const lid_t  len_ga_ = 11;
    const lid_t  len_gs_ = 6;
    const real_t ha_[1]  = {1.0};
    const real_t ga_[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
    const real_t gs_[6]  = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

    const sid_t N() const { return 6; }
    const sid_t Nt() const { return 0; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};
template <>
class Wavelet<6, 2> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 13;
    const lid_t  len_ga_ = 11;
    const lid_t  len_gs_ = 6;
    const real_t ha_[13] = {-3.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -125.0 / 1024.0, 1.0 / 4.0, 181.0 / 256.0, 1.0 / 4.0, -125.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -3.0 / 1024.0};
    const real_t ga_[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
    const real_t gs_[6]  = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

    const sid_t N() const { return 6; }
    const sid_t Nt() const { return 2; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};
template <>
class Wavelet<6, 4> : public InterpolatingWavelet {
    const lid_t  len_ha_ = 17;
    const lid_t  len_ga_ = 11;
    const lid_t  len_gs_ = 6;
    const real_t ha_[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
    const real_t ga_[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
    const real_t gs_[6]  = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

    const sid_t N() const { return 6; }
    const sid_t Nt() const { return 4; }

    // filter
   public:
    const sid_t   len_ha() const override { return len_ha_; };
    const sid_t   len_ga() const override { return len_ga_; };
    const sid_t   len_gs() const override { return len_gs_; };
    const real_t* ha() const override { return ha_; };
    const real_t* ga() const override { return ga_; };
    const real_t* gs() const override { return gs_; };
};

#endif  // SRC_WAVELET_HPP_