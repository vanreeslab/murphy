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

template <lda_t N, lda_t NT>
static constexpr sid_t len_ha_ = 0;

template <lda_t N, lda_t NT>
static constexpr sid_t len_gs_ = 0;

template <lda_t N, lda_t NT>
static constexpr real_t ha_[1] = {0.0};
template <lda_t N, lda_t NT>
static constexpr real_t gs_[1] = {0.0};

//-----------------------------------------------------------------------------
// Wavelet 2.2
template <>
constexpr sid_t len_ha_<2, 2> = 5;
template <>
constexpr sid_t len_gs_<2, 2> = 2;

// static constexpr real_t ha_2_2_[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
// static constexpr real_t gs_2_2_[2] = {0.5, 0.5};

template <>
constexpr real_t ha_<2, 2>[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
template <>
constexpr real_t gs_<2, 2>[2] = {0.5, 0.5};

//-----------------------------------------------------------------------------
// Wavelet 4.0
template <>
constexpr sid_t len_ha_<4, 0> = 1;
template <>
constexpr sid_t len_gs_<4, 0> = 4;

// static constexpr real_t ha_4_0_[1] = {1.0};
// static constexpr real_t gs_4_0_[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
template <>
constexpr real_t ha_<4, 0>[1] = {1.0};
template <>
constexpr real_t gs_<4, 0>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};

//-----------------------------------------------------------------------------
// Wavelet 4.2
template <>
constexpr sid_t len_ha_<4, 2> = 9;
template <>
constexpr sid_t len_gs_<4, 2> = 4;

// static constexpr real_t ha_4_2_[9] = {0.015625, 0.0, -0.125, 0.25, 0.71875, 0.25, -0.125, 0.0, 0.015625};
// static constexpr real_t gs_4_2_[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
template <>
constexpr real_t ha_<4, 2>[9] = {0.015625, 0.0, -0.125, 0.25, 0.71875, 0.25, -0.125, 0.0, 0.015625};
template <>
constexpr real_t gs_<4, 2>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};

//-----------------------------------------------------------------------------
// Wavelet 4.4
template <>
constexpr sid_t len_ha_<4, 4> = 13;
template <>
constexpr sid_t len_gs_<4, 4> = 4;

// static constexpr real_t ha_4_4_[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
// static constexpr real_t gs_4_4_[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};
template <>
constexpr real_t ha_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
template <>
constexpr real_t gs_<4, 4>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};

//-----------------------------------------------------------------------------
// Wavelet 6.0
template <>
constexpr sid_t len_ha_<6, 0> = 1;
template <>
constexpr sid_t len_gs_<6, 0> = 6;

// static constexpr real_t ha_6_0_[1] = {1.0};
// static constexpr real_t gs_6_0_[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
template <>
constexpr real_t ha_<6, 0>[1] = {1.0};
template <>
constexpr real_t gs_<6, 0>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

//-----------------------------------------------------------------------------
// Wavelet 6.2
template <>
constexpr sid_t len_ha_<6, 2> = 13;
template <>
constexpr sid_t len_gs_<6, 2> = 6;

// static constexpr real_t ha_6_2_[13] = {-3.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -125.0 / 1024.0, 1.0 / 4.0, 181.0 / 256.0, 1.0 / 4.0, -125.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -3.0 / 1024.0};
// static constexpr real_t gs_6_2_[6]  = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
template <>
constexpr real_t ha_<6, 2>[13] = {-3.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -125.0 / 1024.0, 1.0 / 4.0, 181.0 / 256.0, 1.0 / 4.0, -125.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -3.0 / 1024.0};
template <>
constexpr real_t gs_<6, 2>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

//-----------------------------------------------------------------------------
// Wavelet 6.4
template <>
constexpr sid_t len_ha_<6, 4> = 17;
template <>
constexpr sid_t len_gs_<6, 4> = 6;

// static constexpr real_t ha_6_4_[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
// static constexpr real_t gs_6_4_[6]  = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
template <>
constexpr real_t ha_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
template <>
constexpr real_t gs_<6, 4>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

/**
 * @brief Implement the interpolating wavelets
 * The wavelet application considered here is based on the lifting/dual lifting scheme.
 * More details on the wavelet approach can be found in the wavelet tutorial (https://github.com/van-Rees-Lab/wavelet_tutorial).
 * 
 * For any couple N,Nt (with Nt <= N) the filter coefficients are computed by (https://github.com/van-Rees-Lab/wavelet_tutorial/blob/master/code/wavelet_xy.py)
 * 
 * @tparam N 
 * @tparam NT
 */
template <lda_t TN = M_WAVELET_N, lda_t TNT = M_WAVELET_NT>
class InterpolatingWavelet : public Wavelet {
   public:
    // lengths
    const sid_t N() const override { return TN; }
    const sid_t Nt() const override { return TNT; }

    // filters
    const sid_t len_ha() const override { return len_ha_<TN, TNT>; };
    const sid_t len_gs() const override { return len_gs_<TN, TNT>; };

    // // filters
    // const real_t* filter_ha() const override {
    //     if constexpr (TN == 2, TNT == 2) {
    //         return ha_2_2_;
    //     } else if (TN == 4, TNT == 0) {
    //         return ha_4_0_;
    //     } else if (TN == 4, TNT == 2) {
    //         return ha_4_2_;
    //     } else if (TN == 4, TNT == 4) {
    //         return ha_4_4_;
    //     } else if (TN == 6, TNT == 0) {
    //         return ha_6_0_;
    //     } else if (TN == 6, TNT == 2) {
    //         return ha_6_2_;
    //     } else if (TN == 6, TNT == 4) {
    //         return ha_6_4_;
    //     }
    // };
    // const real_t* filter_gs() const override {
    //     if constexpr (TN == 2, TNT == 2) {
    //         return gs_2_2_;
    //     } else if (TN == 4, TNT == 0) {
    //         return gs_4_0_;
    //     } else if (TN == 4, TNT == 2) {
    //         return gs_4_2_;
    //     } else if (TN == 4, TNT == 4) {
    //         return gs_4_4_;
    //     } else if (TN == 6, TNT == 0) {
    //         return gs_6_0_;
    //     } else if (TN == 6, TNT == 2) {
    //         return gs_6_2_;
    //     } else if (TN == 6, TNT == 4) {
    //         return gs_6_4_;
    //     }
    // };

   protected:
    /**
     * @brief coarsen the values of the source memory to gather them in the target memory.
     * 
     * @tparam N the number of vanishing moment
     * @tparam Nt the order of interpolation
     * @param ctx the interpolation context
     */
    void Coarsen_(const interp_ctx_t* ctx) const override {
        //-------------------------------------------------------------------------
        // assure alignment for the target, the source, the constant and the temp data
        // m_assume_aligned(ctx->tdata);
        m_assume_aligned(ctx->sdata);
        // m_assume_aligned(ctx->cdata);

        // const lid_t   ha_lim = ha_half_lim();
        // const_mem_ptr ha     = ha_filter();

        constexpr sid_t ha_lim = len_ha_<TN, TNT> / 2;
        const real_t    alpha  = ctx->alpha;

        // const_mem_ptr ha = filter_ha() + ha_lim;
        const_mem_ptr ha = ha_<TN,TNT> + ha_lim;

        for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ++ik2) {
            for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ++ik1) {
                for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ++ik0) {
                    // do some checks
                    m_assert(((2 * ik0 - ha_lim) >= (ctx->srcstart[0])) && ((2 * ik0 + ha_lim) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", 2 * ik0 - ha_lim, ctx->srcstart[0], 2 * ik0 + ha_lim, ctx->srcend[0]);
                    m_assert(((2 * ik1 - ha_lim) >= (ctx->srcstart[1])) && ((2 * ik1 + ha_lim) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", 2 * ik1 - ha_lim, ctx->srcstart[1], 2 * ik1 + ha_lim, ctx->srcend[1]);
                    m_assert(((2 * ik2 - ha_lim) >= (ctx->srcstart[2])) && ((2 * ik2 + ha_lim) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", 2 * ik2 - ha_lim, ctx->srcstart[2], 2 * ik2 + ha_lim, ctx->srcend[2]);
                    //get the local adress of the source, the target and the constant
                    data_ptr       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                    const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                    const data_ptr lsdata = ctx->sdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);

                    // add the constant
                    ltdata[0] = alpha * lcdata[0];
                    // apply the filter
                    for (sid_t id2 = -ha_lim; id2 <= ha_lim; id2++) {
                        for (sid_t id1 = -ha_lim; id1 <= ha_lim; id1++) {
                            for (sid_t id0 = -ha_lim; id0 <= ha_lim; id0++) {
                                const real_t fact = ha[id0] * ha[id1] * ha[id2];
                                ltdata[0] += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                                // check for nan's
                                m_assert(lsdata[0] == lsdata[0], "nan detected");
                                m_assert(fact == fact, "nan detected");
                            }
                        }
                    }
                    m_assert(ltdata[0] == ltdata[0], "nan detected");
                }
            }
        }
        //-------------------------------------------------------------------------
    }

    /**
     * @brief refine the source memory to get the associated target memory information
     * 
     * Here, we assume that the detail coefficients are null.
     * Hence, the values of the function are the scaling coefficient and we simply apply the dual-lifting scheme to obtain the missing information
     * 
     * @tparam N the number of vanishing moment
     * @tparam Nt the order of interpolation
     * @param ctx the interpolation context
     */
    void Refine_(const interp_ctx_t* ctx) const override {
        //-------------------------------------------------------------------------
        // assure alignment for the target, the source, the constant and the temp data
        constexpr sid_t gs_lim = len_gs_<TN, TNT> / 2 - 1;
        const real_t    one    = 1.0;
        const_mem_ptr   gs     = gs_<TN,TNT> + gs_lim;
        const real_t    alpha  = ctx->alpha;

        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        // for each of the data for the needed target
        for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
            for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
                for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
                    // get 0 if odd, 1 if even (even if negative!!)
                    const sid_t iy = m_sign(ik1) * (ik1 % 2);
                    const sid_t ix = m_sign(ik0) * (ik0 % 2);
                    const sid_t iz = m_sign(ik2) * (ik2 % 2);
                    m_assert(ix == 0 || ix == 1, "this are the two possible values");
                    m_assert(iy == 0 || iy == 1, "this are the two possible values");
                    m_assert(iz == 0 || iz == 1, "this are the two possible values");

                    // get the target location
                    data_ptr       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                    const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);

                    //get the local adress of the source, a bit more complicated to handle the negative numbers
                    const lid_t    ik0_s  = (ik0 - ix) / 2;
                    const lid_t    ik1_s  = (ik1 - iy) / 2;
                    const lid_t    ik2_s  = (ik2 - iz) / 2;
                    const data_ptr lsdata = ctx->sdata + m_sidx(ik0_s, ik1_s, ik2_s, 0, ctx->srcstr);
                    m_assert((ik0_s * 2) <= ik0, "if not, we made something wrong...: source = %d, target = %d", ik0_s, ik0);
                    m_assert((ik1_s * 2) <= ik1, "if not, we made something wrong...: source = %d, target = %d", ik1_s, ik1);
                    m_assert((ik2_s * 2) <= ik2, "if not, we made something wrong...: source = %d, target = %d", ik2_s, ik2);

                    // get the filter, depending on if I am odd or even
                    const_mem_ptr gs_x         = (ix == 1) ? (gs) : (&one);
                    const_mem_ptr gs_y         = (iy == 1) ? (gs) : (&one);
                    const_mem_ptr gs_z         = (iz == 1) ? (gs) : (&one);
                    const sid_t   lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
                    const sid_t   lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

                    m_assert(((ik0 / 2 - lim_start[0]) >= ctx->srcstart[0]) && ((ik0 / 2 + lim_end[0]) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                    m_assert(((ik1 / 2 - lim_start[1]) >= ctx->srcstart[1]) && ((ik1 / 2 + lim_end[1]) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                    m_assert(((ik2 / 2 - lim_start[2]) >= ctx->srcstart[2]) && ((ik2 / 2 + lim_end[2]) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);

                    // add the constant array
                    ltdata[0] = alpha * lcdata[0];

                    // if one dim is even, id = 0, -> gs[0] = 1 and that's it
                    // if one dim is odd, id = 1, -> we loop on gs, business as usual
                    for (sid_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                        for (sid_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                            for (sid_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                                const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                                // m_assert(lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)] == lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)], "the error cannot be nan: block @ %d %d %d: %f", ik0 / 2 + id0, ik1 / 2 + id1, ik2 / 2 + id2, lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)]);
                                ltdata[0] += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                            }
                        }
                    }
                    m_assert(ltdata[0] == ltdata[0], "the error cannot be nan: block @ %d %d %d: %f", ik0, ik1, ik2, ltdata[0]);
                    // if (ix == 0 && ix == 1 && ix == 2) {
                    //     m_assert(ltdata[0] == lsdata[0], "ouuups: %e vs %e",ltdata[0], lsdata[0]);
                    // }
                }
            }
        }
        //-------------------------------------------------------------------------
    }

    /**
     * @brief gets the detail coefficients of the wavelet. This approximates the local slope of the data
     * 
     * @tparam order 
     * @param ctx only the trgdata information are used, the source is considered empty
     * @param details_inf_norm the maximum of the local detail coefficients
     */
    void Detail_(const interp_ctx_t* ctx, real_t* details_max) const override {
        //-------------------------------------------------------------------------
        constexpr sid_t  gs_lim = (len_gs_<TN, TNT> / 2 - 1);
        constexpr real_t one    = 1.0;

        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        const_mem_ptr tdata = ctx->tdata;
        const_mem_ptr sdata = ctx->sdata;
        const_mem_ptr gs    = gs_<TN,TNT>  + gs_lim;

        // for each of the data for the considered children
        (*details_max) = 0.0;
        for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
            for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
                for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
                    // get 0 if odd, 1 if even (even if negative!!)
                    const lda_t iy = m_sign(ik1) * (ik1 % 2);
                    const lda_t ix = m_sign(ik0) * (ik0 % 2);
                    const lda_t iz = m_sign(ik2) * (ik2 % 2);
                    m_assert(ix == 0 || ix == 1, "this are the two possible values");
                    m_assert(iy == 0 || iy == 1, "this are the two possible values");
                    m_assert(iz == 0 || iz == 1, "this are the two possible values");

                    // get the nearest even data
                    // const lid_t ik0_s = (ik0 - ix) / 2;
                    // const lid_t ik1_s = (ik1 - iy) / 2;
                    // const lid_t ik2_s = (ik2 - iz) / 2;
                    // const_mem_ptr lsdata = sdata + m_sidx(ik0_s, ik1_s, ik2_s, 0, ctx->srcstr);
                    // const_mem_ptr ltdata = tdata + m_sidx(ik0_s * 2, ik1_s * 2, ik2_s * 2, 0, ctx->trgstr);
                    // m_assert(lsdata[0] == lsdata[0], "the value MUST be the same");

                    const lid_t   ik0_s  = (ik0 - ix);
                    const lid_t   ik1_s  = (ik1 - iy);
                    const lid_t   ik2_s  = (ik2 - iz);
                    const_mem_ptr ltdata = tdata + m_sidx(ik0_s, ik1_s, ik2_s, 0, ctx->trgstr);

                    // get it's location

                    // get the filter, depending on if I am odd or even
                    const_mem_ptr gs_x         = (ix == 1) ? (gs) : (&one);
                    const_mem_ptr gs_y         = (iy == 1) ? (gs) : (&one);
                    const_mem_ptr gs_z         = (iz == 1) ? (gs) : (&one);
                    const sid_t   lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
                    const sid_t   lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

                    // m_assert(((ik0_s - lim_start[0]) >= ctx->srcstart[0]) && ((ik0_s + lim_end[0]) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0_s - lim_start[0], ctx->srcstart[0], ik0_s + lim_end[0], ctx->srcend[0]);
                    // m_assert(((ik1_s - lim_start[1]) >= ctx->srcstart[1]) && ((ik1_s + lim_end[1]) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1_s - lim_start[1], ctx->srcstart[1], ik1_s + lim_end[1], ctx->srcend[1]);
                    // m_assert(((ik2_s - lim_start[2]) >= ctx->srcstart[2]) && ((ik2_s + lim_end[2]) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2_s - lim_start[2], ctx->srcstart[2], ik2_s + lim_end[2], ctx->srcend[2]);

                    // if one dim is even, id = 0, -> gs[0] = 1 and that's it
                    // if one dim is odd, id = 1, -> we loop on gs, business as usual
                    real_t interp = 0.0;
                    for (sid_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                        for (sid_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                            for (sid_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                                const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                                interp += fact * ltdata[m_sidx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->trgstr)];
                                // interp += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                                // m_assert(lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)] == ltdata[m_sidx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->trgstr)], "this must be @%d %d %d + %d %d %d: %e vs %e", ik0, ik1, ik2, id0, id1, id2, lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)], ltdata[m_sidx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->trgstr)]);
                            }
                        }
                    }
                    real_t detail = ctx->tdata[m_sidx(ik0, ik1, ik2, 0, ctx->trgstr)] - interp;

                    // check the maximum
                    (*details_max) = m_max(fabs(detail), (*details_max));
                }
            }
        }
        //-------------------------------------------------------------------------
    }
};

#endif  // SRC_WAVELET_HPP_