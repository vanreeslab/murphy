#ifndef SRC_INTERPOLATING_WAVELET_HPP_
#define SRC_INTERPOLATING_WAVELET_HPP_

#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "gridblock.hpp"
#include "wavelet.hpp"

template <lda_t N, lda_t NT>
static constexpr short_t len_ha_ = 0;  //!< length of the Ha filter
template <lda_t N, lda_t NT>
static constexpr short_t len_ga_ = 0;  //!< length of the S dual filter
template <lda_t N, lda_t NT>
static constexpr short_t len_js_ = 0;  //!< length of the Js filter
template <lda_t N, lda_t NT>
static constexpr short_t len_ks_ = 0;  //!< length of the Ks filter

/** @brief Ha filter, used to compute the coarser scaling coefficients */
template <lda_t N, lda_t NT>
static constexpr real_t ha_[1] = {0.0};

/** @brief S dual filter, substitutes the Ga filter and the Ks filter when details are null
 * 
 * - used to compute the detail coefficients Ga = 1 - S_dual
 * - when the details are null, the finer scaling can also be obtained using the S_dual
 */
template <lda_t N, lda_t NT>
static constexpr real_t ga_[1] = {0.0};

/**
 * @brief Js filter: used to compute the EVEN scaling coefficients from the coarse value + details
 * 
 * - in practice, only needed for the smoothing, to compute the influence of a detail modification to the fine scaling
 */
template <lda_t N, lda_t NT>
static constexpr real_t js_[1] = {0.0};

/**
 * @brief Ks filter: used to compute the ODD scaling coefficients from the coarse value + details
 * 
 * - in practice, only needed for the smoothing, to compute the influence of a detail modification to the fine scaling
 */
template <lda_t N, lda_t NT>
static constexpr real_t ks_[1] = {0.0};

/**
 * @brief constant in front of the epsilon 
 */
template <lda_t N, lda_t NT>
static constexpr real_t eps_c = 1.0;

//-----------------------------------------------------------------------------
// Wavelet 2.0
template <>
inline constexpr short_t len_ha_<2, 0> = 1;
template <>
inline constexpr short_t len_ga_<2, 0> = 3;
template <>
inline constexpr short_t len_js_<2, 0> = 1;
template <>
inline constexpr short_t len_ks_<2, 0> = 3;
template <>
inline constexpr real_t ha_<2, 0>[1] = {1.0};
template <>
inline constexpr real_t ga_<2, 0>[3] = {-0.5, 1.0, -0.5};
template <>
inline constexpr real_t js_<2, 0>[1] = {1.0};
template <>
inline constexpr real_t ks_<2, 0>[3] = {0.5, 1.0, 0.5};
template <>
inline constexpr real_t eps_c<2, 0> = 7.0;

//-----------------------------------------------------------------------------
// Wavelet 2.2
template <>
inline constexpr short_t len_ha_<2, 2> = 5;
template <>
inline constexpr short_t len_ga_<2, 2> = 3;
template <>
inline constexpr short_t len_js_<2, 2> = 3;
template <>
inline constexpr short_t len_ks_<2, 2> = 5;
template <>
inline constexpr real_t ha_<2, 2>[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
template <>
inline constexpr real_t ga_<2, 2>[3] = {-0.5, 1.0, -0.5};
template <>
inline constexpr real_t js_<2, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
inline constexpr real_t ks_<2, 2>[5] = {-1.0 / 8.0, 1.0 / 2.0, 3.0 / 4.0, 1.0 / 2.0, -1.0 / 8.0};
template <>
inline constexpr real_t eps_c<2, 2> = 7.0;

//-----------------------------------------------------------------------------
// Wavelet 4.0
template <>
inline constexpr short_t len_ha_<4, 0> = 1;
template <>
inline constexpr short_t len_ga_<4, 0> = 7;
template <>
inline constexpr short_t len_js_<4, 0> = 1;
template <>
inline constexpr short_t len_ks_<4, 0> = 7;
template <>
inline constexpr real_t ha_<4, 0>[1] = {1.0};
template <>
inline constexpr real_t ga_<4, 0>[7] = {1.0 / 16.0, 0.0, -9.0 / 16.0, 1.0, -9.0 / 16.0, 0.0, 1.0 / 16.0};
template <>
inline constexpr real_t js_<4, 0>[1] = {1.0};
template <>
inline constexpr real_t ks_<4, 0>[7] = {-1.0 / 16.0, 0.0, 9.0 / 16.0, 1.0, 9.0 / 16.0, 0.0, -1.0 / 16.0};
template <>
inline constexpr real_t eps_c<4, 0> = 9.4375;

//-----------------------------------------------------------------------------
// Wavelet 4.2
template <>
inline constexpr short_t len_ha_<4, 2> = 9;
template <>
inline constexpr short_t len_ga_<4, 2> = 7;
template <>
inline constexpr short_t len_js_<4, 2> = 3;
template <>
inline constexpr short_t len_ks_<4, 2> = 9;
template <>
inline constexpr real_t ha_<4, 2>[9] = {1.0 / 64.0, 0.0, -1.0 / 8.0, 1.0 / 4.0, 23.0 / 32.0, 1.0 / 4.0, -1.0 / 8.0, 0.0, 1.0 / 64.0};
template <>
inline constexpr real_t ga_<4, 2>[7] = {1.0 / 16.0, 0.0, -9.0 / 16.0, 1.0, -9.0 / 16.0, 0.0, 1.0 / 16.0};
template <>
inline constexpr real_t js_<4, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
inline constexpr real_t ks_<4, 2>[9] = {1.0 / 64.0, -1.0 / 16.0, -1.0 / 8.0, 9.0 / 16.0, 23.0 / 32.0, 9.0 / 16.0, -1.0 / 8.0, -1.0 / 16.0, 1.0 / 64.0};
template <>
inline constexpr real_t eps_c<4, 2> = 9.4375;

// //-----------------------------------------------------------------------------
// // Wavelet 4.4
// template <>
// constexpr short_t len_ha_<4, 4> = 13;
// template <>
// constexpr short_t len_ga_<4, 4> = 4;
// template <>
// constexpr short_t len_js_<4, 4> = 7;
// template <>
// constexpr short_t len_ks_<4, 4> = 13;
// template <>
// constexpr real_t ha_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
// template <>
// constexpr real_t ga_<4, 4>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
// template <>
// constexpr real_t js_<4, 4>[7] = {1.0 / 32.0, 0.0, -9.0 / 32.0, 1.0, -9.0 / 32.0, 0.0, 1.0 / 32.0};
// template <>
// constexpr real_t ks_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 16.0, -63.0 / 512.0, 9 / 16.0, 87.0 / 128.0, 9.0 / 16.0, -63.0 / 512.0, -1.0 / 16.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};

//-----------------------------------------------------------------------------
// Wavelet 6.0
template <>
inline constexpr short_t len_ha_<6, 0> = 1;
template <>
inline constexpr short_t len_ga_<6, 0> = 11;
template <>
inline constexpr short_t len_js_<6, 0> = 1;
template <>
inline constexpr short_t len_ks_<6, 0> = 11;
template <>
inline constexpr real_t ha_<6, 0>[1] = {1.0};
template <>
inline constexpr real_t ga_<6, 0>[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
template <>
inline constexpr real_t js_<6, 0>[1] = {1.0};
template <>
inline constexpr real_t ks_<6, 0>[11] = {3.0 / 256.0, 0.0, -25.0 / 256.0, 0.0, 75.0 / 128.0, 1.0, 75.0 / 128.0, 0.0, -25.0 / 256.0, 0.0, 3.0 / 256.0};
template <>
inline constexpr real_t eps_c<6, 0> = 10.973388671875;

//-----------------------------------------------------------------------------
// Wavelet 6.2
template <>
inline constexpr short_t len_ha_<6, 2> = 13;
template <>
inline constexpr short_t len_ga_<6, 2> = 11;
template <>
inline constexpr short_t len_js_<6, 2> = 3;
template <>
inline constexpr short_t len_ks_<6, 2> = 13;
template <>
inline constexpr real_t ha_<6, 2>[13] = {-3.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -125.0 / 1024.0, 1.0 / 4.0, 181.0 / 256.0, 1.0 / 4.0, -125.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -3.0 / 1024.0};
template <>
inline constexpr real_t ga_<6, 2>[11] = {-3.0 / 256.0, 0.0, 25.0 / 256.0, 0.0, -75.0 / 128.0, 1.0, -75.0 / 128.0, 0.0, 25.0 / 256.0, 0.0, -3.0 / 256.0};
template <>
inline constexpr real_t js_<6, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
inline constexpr real_t ks_<6, 2>[13] = {-3.0 / 1024.0, 3.0 / 256.0, 11.0 / 512.0, -25.0 / 256.0, -125.0 / 1024.0, 75.0 / 128.0, 181.0 / 256.0, 75.0 / 128.0, -125.0 / 1024.0, -25.0 / 256.0, 11.0 / 512.0, 3.0 / 256.0, -3.0 / 1024.0};
template <>
inline constexpr real_t eps_c<6, 2> = 10.973388671875;

// //-----------------------------------------------------------------------------
// // Wavelet 6.4
// template <>
// constexpr short_t len_ha_<6, 4> = 17;
// template <>
// constexpr short_t len_ga_<6, 4> = 6;
// template <>
// constexpr short_t len_js_<6, 4> = 7;
// template <>
// constexpr short_t len_ks_<6, 4> = 17;
// template <>
// constexpr real_t ha_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
// template <>
// constexpr real_t ga_<6, 4>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
// template <>
// constexpr real_t js_<6, 4>[7] = {1.0 / 32.0, 0.0, -9.0 / 32.0, 1.0, -9.0 / 32.0, 0.0, 1.0 / 32.0};
// template <>
// constexpr real_t ks_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 3.0 / 256.0, 87.0 / 2048.0, -25.0 / 256.0, -243.0 / 2048.0, 75.0 / 128.0, 2721.0 / 4096.0, 75.0 / 128.0, -243.0 / 2048.0, -25.0 / 256.0, 87.0 / 2048.0, 3.0 / 256.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};

/**
 * @brief Implement the interpolating wavelets from the familly N.Nt (with Nt <= N)
 * 
 * The wavelet application considered here is based on the lifting/dual lifting scheme.
 * More details on the wavelet approach can be found in the wavelet tutorial (https://github.com/van-Rees-Lab/wavelet_tutorial).
 * For any couple N,Nt (with Nt <= N) the filter coefficients are computed by (https://github.com/van-Rees-Lab/wavelet_tutorial/blob/master/code/wavelet_xy.py)
 * 
 * @tparam TN  the order of the interpolation
 * @tparam TNT the number of conserved moments
 */
template <short_t TN = M_WAVELET_N, short_t TNT = M_WAVELET_NT>
class InterpolatingWavelet : public Wavelet {
   public:
    const short_t N() const override { return TN; }
    const short_t Nt() const override { return TNT; }
    const short_t len_ha() const override { return len_ha_<TN, TNT>; };
    const short_t len_ga() const override { return len_ga_<TN, TNT>; };
    const short_t len_js() const override { return len_js_<TN, TNT>; };
    const short_t len_ks() const override { return len_ks_<TN, TNT>; };
    const real_t  eps_const() const override { return eps_c<TN, TNT>; };

   protected:
    /**
     * @brief coarsen the values of the source memory to gather them in the target memory.
     * 
     * The details are destroyed whatever value they have. we use the ga filter
     * 
     * @param ctx the interpolation context
     */
    void Coarsen_(const interp_ctx_t* const  ctx) const override {
        //-------------------------------------------------------------------------
        // filter crap
        constexpr short_t   ha_lim = len_ha_<TN, TNT> / 2;
        const real_t* const ha     = ha_<TN, TNT> + ha_lim;

        // get pointers
        real_t* const       tdata = ctx->tdata.Write();
        const real_t* const sdata = ctx->sdata.Read();

        // define the lambda, only tdata is changed
        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // do some checks on the bounds
            m_assert(((2 * i0 - ha_lim) >= (ctx->srcstart[0])) && ((2 * i0 + ha_lim) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", 2 * i0 - ha_lim, ctx->srcstart[0], 2 * i0 + ha_lim, ctx->srcend[0]);
            m_assert(((2 * i1 - ha_lim) >= (ctx->srcstart[1])) && ((2 * i1 + ha_lim) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", 2 * i1 - ha_lim, ctx->srcstart[1], 2 * i1 + ha_lim, ctx->srcend[1]);
            m_assert(((2 * i2 - ha_lim) >= (ctx->srcstart[2])) && ((2 * i2 + ha_lim) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", 2 * i2 - ha_lim, ctx->srcstart[2], 2 * i2 + ha_lim, ctx->srcend[2]);

            //get the local adress of the source, which is twice finer compared to the target
            const real_t* const lsdata = sdata + m_idx(2 * i0, 2 * i1, 2 * i2, 0, ctx->srcstr);

            // apply the filter
            real_t value = 0.0;
            for (short_t id2 = -ha_lim; id2 <= ha_lim; ++id2) {
                for (short_t id1 = -ha_lim; id1 <= ha_lim; ++id1) {
                    for (short_t id0 = -ha_lim; id0 <= ha_lim; ++id0) {
                        // add the info
                        const real_t fact = ha[id0] * ha[id1] * ha[id2];
                        value += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];

                        // check for nan's
                        m_assert(std::isfinite(lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)]), "nan detected");
                        m_assert(std::isfinite(fact), "nan detected");
                    }
                }
            }
            // check for nan's
            m_assert(std::isfinite(value), "the value cannot be nan: block @ %d %d %d: %f", i0, i1, i2, value);

            tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)] = value;
        };

        // get the start and end
        const bidx_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const bidx_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        // run that shit
        for_loop(&op, start, end);
        //-------------------------------------------------------------------------
    }

    /**
     * @brief refine the source memory to get the associated target memory information, assuming no details
     * 
     * we assume that the detail coefficients are null and therefore we skip every non-needed coefficient.
     * It also means that the lifting step is useless and the scalings are simply copied
     *
     * @param ctx the interpolation context
     */
    void RefineZeroDetails_(const interp_ctx_t* const  ctx) const override {
        m_assert(ctx->alpha == 0.0, "the alpha = %e must be = 0.0", ctx->alpha);
        //-------------------------------------------------------------------------
        constexpr short_t   ks_lim = len_ks_<TN, TNT> / 2;
        const real_t* const ks     = ks_<TN, TNT> + ks_lim;
        constexpr short_t   js_lim = len_js_<TN, TNT> / 2;
        const real_t* const js     = js_<TN, TNT> + js_lim;

        // get the pointers
        real_t* const       tdata = ctx->tdata.Write();
        const real_t* const sdata = ctx->sdata.Read();

        // define the lambda, only tdata is changed
        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get if we are odd or even in the current location
            const short_t odd_x = m_sign(i0) * (i0 % 2);
            const short_t odd_y = m_sign(i1) * (i1 % 2);
            const short_t odd_z = m_sign(i2) * (i2 % 2);
            m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
            m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
            m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");

            //get the start adress of the source, we need to shift by one backward if we are odd (see the filter loop)
            const lid_t i0_s = (i0 - odd_x) / 2;
            const lid_t i1_s = (i1 - odd_y) / 2;
            const lid_t i2_s = (i2 - odd_z) / 2;
            m_assert((i0_s * 2) <= i0, "if not, we made something wrong...: source = %d, target = %d", i0_s, i0);
            m_assert((i1_s * 2) <= i1, "if not, we made something wrong...: source = %d, target = %d", i1_s, i1);
            m_assert((i2_s * 2) <= i2, "if not, we made something wrong...: source = %d, target = %d", i2_s, i2);

            const real_t* const lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);

            // get the filter, depending on if I am odd or even
            const real_t* const ks_x = (odd_x) ? (ks) : (js);
            const real_t* const ks_y = (odd_y) ? (ks) : (js);
            const real_t* const ks_z = (odd_z) ? (ks) : (js);

            // we start at the last scaling:
            //  - if we are odd = the last odd number
            //  - if we are even = the last even number
            const bidx_t lim[3] = {
                (odd_x) ? m_max(0, ks_lim - (1 - ks_lim % 2)) : m_max(0, js_lim - (js_lim % 2)),
                (odd_y) ? m_max(0, ks_lim - (1 - ks_lim % 2)) : m_max(0, js_lim - (js_lim % 2)),
                (odd_z) ? m_max(0, ks_lim - (1 - ks_lim % 2)) : m_max(0, js_lim - (js_lim % 2))};

            m_assert(((i0_s + (-lim[0] + 1) / 2) >= ctx->srcstart[0]) && ((i0_s + (lim[0] + 1) / 2) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", (i0_s + (-lim[0] + 1) / 2), ctx->srcstart[0], (i0_s + (lim[0] + 1) / 2), ctx->srcend[0]);
            m_assert(((i1_s + (-lim[1] + 1) / 2) >= ctx->srcstart[1]) && ((i1_s + (lim[1] + 1) / 2) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", (i1_s + (-lim[1] + 1) / 2), ctx->srcstart[1], (i1_s + (lim[1] + 1) / 2), ctx->srcend[1]);
            m_assert(((i2_s + (-lim[2] + 1) / 2) >= ctx->srcstart[2]) && ((i2_s + (lim[2] + 1) / 2) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", (i2_s + (-lim[2] + 1) / 2), ctx->srcstart[2], (i2_s + (lim[2] + 1) / 2), ctx->srcend[2]);

            // if we are a detail, we never pass by the index 0 as we start from an odd number
            // the source is taken as (i+1)/2 so that when i=-1, we access the position 0 and when i=1, we access the position 1
            real_t value = 0.0;
            for (bidx_t id2 = -lim[2]; id2 <= lim[2]; id2 += 2) {
                for (bidx_t id1 = -lim[1]; id1 <= lim[1]; id1 += 2) {
                    for (bidx_t id0 = -lim[0]; id0 <= lim[0]; id0 += 2) {
                        const real_t fact   = ks_x[id0] * ks_y[id1] * ks_z[id2];
                        const real_t svalue = lsdata[m_idx((id0 + 1) / 2, (id1 + 1) / 2, (id2 + 1) / 2, 0, ctx->srcstr)];
                        value += fact * svalue;

                        // check for nan's
                        m_assert(std::isfinite(svalue), "nan detected");
                        m_assert(std::isfinite(fact), "nan detected");
                    }
                }
            }
            m_assert(std::isfinite(value), "the value cannot be nan: block @ %d %d %d: %f", i0, i1, i2, value);
            tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)] = value;
        };

        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
        for_loop(&op, start, end);
        //-------------------------------------------------------------------------
    };

    /**
     * @brief overwrites the data to ensure 0-detail coefficients (it uses the dual lifting only!)
     * 
     * @warning Due to the interdependence between the detail coefficients, we need to update dimension by dimension:
     * first dx, dy, dz, then dxy, dyz and dxz and finally dxyz
     * 
     * @param ctx 
     */
    void OverwriteDetailsDualLifting_(const interp_ctx_t* const  ctx) const override {
        m_assert(ctx->alpha == 0.0, "the alpha = %e must be = 0.0", ctx->alpha);
        m_assert((len_ga_<TN, TNT> % 2) == 1, "the filter size must be odd");
        //-------------------------------------------------------------------------

        // get the filter of dual lifting = filter of lifting, the "-" is done in the loop
        constexpr real_t    one    = 1.0;
        constexpr short_t   ga_lim = len_ga_<TN, TNT> / 2;
        const real_t* const ga     = ga_<TN, TNT> + ga_lim;
        const real_t* const iga    = ga;  // + ga_lim;

        // get the pointers
        real_t* const       tdata = ctx->tdata.Write();
        const real_t* const sdata = ctx->sdata.Read();

        // count tracks how many dimension is done at the same time.
        // if count = 1, we do the dx, dy, dz
        // if count = 2, we do the dxy, dyz and dxz
        // if count = 3, we do the dxyz
        short_t count = 0;

        // define the lambda, only tdata is changed and count must be taken by ref (the copy always has count = 0)
        auto op = [=, &tdata, &count](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get if we are odd or even in the current location
            const short_t odd_x = m_sign(i0) * (i0 % 2);
            const short_t odd_y = m_sign(i1) * (i1 % 2);
            const short_t odd_z = m_sign(i2) * (i2 % 2);
            m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
            m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
            m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");

            // if we do not process the current case, return
            if ((odd_x + odd_y + odd_z) != count) {
                return;
            }

            // get the filter, depending on if I am odd or even
            const real_t* const iga_x = (odd_x) ? (iga) : (&one);
            const real_t* const iga_y = (odd_y) ? (iga) : (&one);
            const real_t* const iga_z = (odd_z) ? (iga) : (&one);

            // get the limit for the different filters
            // if even we use a lim of 0 as the bounds in the for loops are "="
            const bidx_t lim[3] = {
                (odd_x) ? m_max(0, ga_lim) : 0,
                (odd_y) ? m_max(0, ga_lim) : 0,
                (odd_z) ? m_max(0, ga_lim) : 0};

            m_assert(((i0 - lim[0]) >= ctx->srcstart[0]) && ((i0 + lim[0]) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", (i0 - lim[0]), ctx->srcstart[0], (i0 + lim[0]), ctx->srcend[0]);
            m_assert(((i1 - lim[1]) >= ctx->srcstart[1]) && ((i1 + lim[1]) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", (i1 - lim[1]), ctx->srcstart[1], (i1 + lim[1]), ctx->srcend[1]);
            m_assert(((i2 - lim[2]) >= ctx->srcstart[2]) && ((i2 + lim[2]) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", (i2 - lim[2]), ctx->srcstart[2], (i2 + lim[2]), ctx->srcend[2]);

            const real_t* const lsdata = sdata + m_idx(i0, i1, i2, 0, ctx->srcstr);
            real_t              value  = 0.0;

            for (bidx_t id2 = -lim[2]; id2 <= lim[2]; ++id2) {
                for (bidx_t id1 = -lim[1]; id1 <= lim[1]; ++id1) {
                    for (bidx_t id0 = -lim[0]; id0 <= lim[0]; ++id0) {
                        // if I consider myself, don't take me into account
                        const bool to_skip = (fabs(id0) + fabs(id1) + fabs(id2)) == 0.0;
                        // minus to inverse the dual lifting
                        const real_t fact   = -(iga_x[id0] * iga_y[id1] * iga_z[id2]) * (!to_skip);
                        const real_t svalue = lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];
                        // sum-up
                        value += fact * svalue;

                        // check for nan's
                        m_assert(std::isfinite(svalue), "nan detected");
                        m_assert(std::isfinite(fact), "nan detected");
                    }
                }
            }
            m_assert(std::isfinite(value), "the value cannot be nan: block @ %d %d %d: %f", i0, i1, i2, value);
            tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)] = value;
        };

        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
        // count = 1: dx, dy and dz
        count = 1;
        for_loop(&op, start, end);
        // count = 2: dxy, dxz and dyz
        count = 2;
        for_loop(&op, start, end);
        // count = 3: dxyz
        count = 3;
        for_loop(&op, start, end);
        //-------------------------------------------------------------------------
    };

    /**
     * @brief compute the details over a memory block and store the maximum infinite norm of the details
     * 
     * Depending on the values of @ref interp_ctx_t::tdata and @ref interp_ctx_t::alpha, the detail coefficients are stored.
     *     - if the @ref interp_ctx_t::tdata is empty, nothing is stored
     *     - if the @ref interp_ctx_t::tdata is not empty and @ref interp_ctx_t::alpha >= 0, the value is stored if the abs(detail) < ctx->alpha
     *     - if the @ref interp_ctx_t::tdata is not empty and @ref interp_ctx_t::alpha < 0, the value is stored as tdata * detail
     * 
     * @param ctx the interpolation context
     * @param details_max the maximum of the local detail coefficients
     */
    void Detail_(const interp_ctx_t* const  ctx, real_t* const  details_max) const override {
        m_assert(details_max[0] == 0.0, "the value must be 0.0");
        //-------------------------------------------------------------------------
        constexpr short_t   ha_lim = len_ha_<TN, TNT> / 2;
        constexpr short_t   ga_lim = len_ga_<TN, TNT> / 2;
        const real_t* const ha     = ha_<TN, TNT> + ha_lim;
        const real_t* const ga     = ga_<TN, TNT> + ga_lim;

        // check if we need to store some value -> get the target pointer
        bool          store         = !(ctx->tdata.IsEmpty());
        bool          store_on_mask = (ctx->alpha < 0.0) && store;
        real_t        temp          = 0.0;
        real_t* const tdata         = (store) ? (ctx->tdata.Write()) : (&temp);

        // m_log("store? %d, store on mask? %d", store, store_on_mask);
        // m_log("computing details using %d %d %d to %d %d %d", ctx->srcstart[0], ctx->srcstart[1], ctx->srcstart[2],
        //       ctx->srcend[0], ctx->srcend[1], ctx->srcend[2]);
        // m_log("and computing details from %d %d %d to %d %d %d", ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2],
        //       ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]);

        // get the source pointer
        const real_t* const sdata = ctx->sdata.Read();

        // go
        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get if we are odd or even in the cupgridblockrrent location
            const short_t odd_x = m_sign(i0) * (i0 % 2);
            const short_t odd_y = m_sign(i1) * (i1 % 2);
            const short_t odd_z = m_sign(i2) * (i2 % 2);
            m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
            m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
            m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");

            const real_t* const lsdata = sdata + m_idx(i0, i1, i2, 0, ctx->srcstr);

            // get the filters
            const real_t* const f_x = (odd_x == 1) ? (ga) : (ha);
            const real_t* const f_y = (odd_y == 1) ? (ga) : (ha);
            const real_t* const f_z = (odd_z == 1) ? (ga) : (ha);

            // get the limits, if we are a scaling coefficent, we don't care as we know it's 0.0 -> bypass the loop
            const bool   is_scaling = (!odd_x) && (!odd_y) && (!odd_z);
            const bidx_t lim[3]     = {
                (is_scaling) ? (-1) : ((odd_x) ? (ga_lim) : (ha_lim)),
                (is_scaling) ? (-1) : ((odd_y) ? (ga_lim) : (ha_lim)),
                (is_scaling) ? (-1) : ((odd_z) ? (ga_lim) : (ha_lim))};

            m_assert(((i0 - lim[0]) >= ctx->srcstart[0]) && ((i0 + lim[0]) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0 - lim[0], ctx->srcstart[0], i0 + lim[0], ctx->srcend[0]);
            m_assert(((i1 - lim[1]) >= ctx->srcstart[1]) && ((i1 + lim[1]) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1 - lim[1], ctx->srcstart[1], i1 + lim[1], ctx->srcend[1]);
            m_assert(((i2 - lim[2]) >= ctx->srcstart[2]) && ((i2 + lim[2]) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2 - lim[2], ctx->srcstart[2], i2 + lim[2], ctx->srcend[2]);

            // if one dim is even, id = 0, -> gs[0] = 1 and that's it
            // if one dim is odd, id = 1, -> we loop on gs, business as usual
            real_t detail = 0.0;
            for (bidx_t id2 = -lim[2]; id2 <= lim[2]; ++id2) {
                for (bidx_t id1 = -lim[1]; id1 <= lim[1]; ++id1) {
                    for (bidx_t id0 = -lim[0]; id0 <= lim[0]; ++id0) {
                        const real_t fact = f_x[id0] * f_y[id1] * f_z[id2];
                        detail += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];
                    }
                }
            }
            m_assert(std::isfinite(detail), "the data in %d %d %d is nan %e", i0, i1, i2, detail);

            // check the maximum
            (*details_max) = m_max(fabs(detail), (*details_max));

            // store if needed, the index is 0 if not store
            const bidx_t store_id = store * m_idx(i0, i1, i2, 0, ctx->trgstr);
            const real_t value    = store_on_mask ? (detail * tdata[store_id]) : (detail * (fabs(detail) < ctx->alpha));
            tdata[store_id]       = value;

            // if (fabs(value) > 1e-13) {
            //     m_log("registering value for %d %d %d", i0, i1, i2);
            // }
        };

        // reset the detail max (to be sure)
        *details_max = 0.0;
        // let's go
        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
        for_loop(&op, start, end);
    };

    void ForwardWaveletTransform_(const interp_ctx_t* const  ctx, real_t* const  details_max) const override {
        m_assert(*details_max == 0.0, "the value must be 0.0");
        //-------------------------------------------------------------------------
        constexpr short_t   ha_lim = len_ha_<TN, TNT> / 2;
        constexpr short_t   ga_lim = len_ga_<TN, TNT> / 2;
        const real_t* const ha     = ha_<TN, TNT> + ha_lim;
        const real_t* const ga     = ga_<TN, TNT> + ga_lim;

        // check if we need to store some value -> get the target pointer

        // get the source pointer
        const real_t* const sdata = ctx->sdata.Read();
        real_t* const       tdata = ctx->tdata.Write();

        // go
        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get if we are odd or even in the current location
            const short_t odd_x = m_sign(i0) * (i0 % 2);
            const short_t odd_y = m_sign(i1) * (i1 % 2);
            const short_t odd_z = m_sign(i2) * (i2 % 2);
            m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
            m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
            m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");

            const real_t* const lsdata = sdata + m_idx(i0, i1, i2, 0, ctx->srcstr);

            // get the filters
            const real_t* const f_x = (odd_x == 1) ? (ga) : (ha);
            const real_t* const f_y = (odd_y == 1) ? (ga) : (ha);
            const real_t* const f_z = (odd_z == 1) ? (ga) : (ha);

            // get the limits of this
            const bidx_t lim[3] = {((odd_x) ? (ga_lim) : (ha_lim)),
                                   ((odd_y) ? (ga_lim) : (ha_lim)),
                                   ((odd_z) ? (ga_lim) : (ha_lim))};

            m_assert(((i0 - lim[0]) >= ctx->srcstart[0]) && ((i0 + lim[0]) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0 - lim[0], ctx->srcstart[0], i0 + lim[0], ctx->srcend[0]);
            m_assert(((i1 - lim[1]) >= ctx->srcstart[1]) && ((i1 + lim[1]) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1 - lim[1], ctx->srcstart[1], i1 + lim[1], ctx->srcend[1]);
            m_assert(((i2 - lim[2]) >= ctx->srcstart[2]) && ((i2 + lim[2]) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2 - lim[2], ctx->srcstart[2], i2 + lim[2], ctx->srcend[2]);

            real_t wave_data = 0.0;
            for (bidx_t id2 = -lim[2]; id2 <= lim[2]; ++id2) {
                for (bidx_t id1 = -lim[1]; id1 <= lim[1]; ++id1) {
                    for (bidx_t id0 = -lim[0]; id0 <= lim[0]; ++id0) {
                        const real_t fact = f_x[id0] * f_y[id1] * f_z[id2];
                        wave_data += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];
                    }
                }
            }
            m_assert(std::isfinite(wave_data), "the data in %d %d %d is nan %e", i0, i1, i2, wave_data);

            // store the value
            tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)] = wave_data;

            // Store the max if we are a detail
            const bool   is_scaling = (!odd_x) && (!odd_y) && (!odd_z);
            const real_t detail     = (is_scaling) ? 0.0 : wave_data;
            (*details_max)          = m_max(fabs(detail), (*details_max));
        };

        // reset the detail max (to be sure)
        *details_max = 0.0;

        // compute the start and end indexes, we need more details than the number of scalings
        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
        for_loop(&op, start, end);
    };

    /**
     * @brief Remove from the target the non-zero detail coefficients using the inverse wavelet transform
     * 
     * It first computes the correction as the inverse wavelet transform of the detail coefficients,
     * and then remove the correction from the target
     * 
     * @param ctx the interpolation context
     */
    void Smooth_(const interp_ctx_t* const  ctx) const override {
        //-------------------------------------------------------------------------
        // filters
        constexpr short_t   js_lim = len_js_<TN, TNT> / 2;
        constexpr short_t   ks_lim = len_ks_<TN, TNT> / 2;
        const real_t* const js     = js_<TN, TNT> + js_lim;
        const real_t* const ks     = ks_<TN, TNT> + ks_lim;

        // get the pointers
        real_t* const       tdata = ctx->tdata.Write();
        const real_t* const ddata = ctx->sdata.Read();

        // go, only tada is changed
        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const bool odd_x = m_sign(i0) * (i0 % 2);
            const bool odd_y = m_sign(i1) * (i1 % 2);
            const bool odd_z = m_sign(i2) * (i2 % 2);
            m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
            m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
            m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");

            // get the filters
            const real_t* const f_x = (odd_x) ? (ks) : (js);
            const real_t* const f_y = (odd_y) ? (ks) : (js);
            const real_t* const f_z = (odd_z) ? (ks) : (js);

            // get the limits
            const bidx_t lim[3] = {
                (js_lim) * (!odd_x) + (ks_lim) * (odd_x),
                (js_lim) * (!odd_y) + (ks_lim) * (odd_y),
                (js_lim) * (!odd_z) + (ks_lim) * (odd_z)};

            m_assert(((i0 - lim[0]) >= ctx->srcstart[0]) && ((i0 + lim[0]) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0 - lim[0], ctx->srcstart[0], i0 + lim[0], ctx->srcend[0]);
            m_assert(((i1 - lim[1]) >= ctx->srcstart[1]) && ((i1 + lim[1]) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1 - lim[1], ctx->srcstart[1], i1 + lim[1], ctx->srcend[1]);
            m_assert(((i2 - lim[2]) >= ctx->srcstart[2]) && ((i2 + lim[2]) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2 - lim[2], ctx->srcstart[2], i2 + lim[2], ctx->srcend[2]);

            // get the local datassss
            real_t* const       ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
            const real_t* const lddata = ddata + m_idx(i0, i1, i2, 0, ctx->srcstr);

            // let's go tocard
            real_t corr = 0.0;
            for (bidx_t id2 = -lim[2]; id2 <= lim[2]; ++id2) {
                for (bidx_t id1 = -lim[1]; id1 <= lim[1]; ++id1) {
                    for (bidx_t id0 = -lim[0]; id0 <= lim[0]; ++id0) {
                        const real_t fact = (f_x[id0]) * (f_y[id1]) * (f_z[id2]);
                        corr += fact * lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)];
                        // sanity check
                        m_assert(!((id0 + i0) % 2 == 0 && (id1 + i1) % 2 == 0 && (id2 + i2) % 2 == 0 && fabs(lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]) > 1e-16), "detail is wrong: %e, it should be null at even locations: %d %d %d", lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)], i0 + id0, i1 + id1, i2 + id2);
                    }
                }
            }
            // // get the target location
            // if (0 <= i0 && i0 < M_N && 0 <= i1 && i1 < M_N && 0 <= i2 && i2 < M_N) {
            //     m_log("@ %d %d %d : smoothed %e -> %e", i0, i1, i2, ltdata[0], ltdata[0] - corr);
            // }
            m_assert(std::isfinite(corr), "the data in %d %d %d is nan %e", i0, i1, i2, corr);
            ltdata[0] -= corr;

            // if (fabs(corr) > 1e-3) {
            //     m_log("correction is %e at %d %d %d", corr, i0, i1, i2);
            // }
        };

        const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
        for_loop(&op, start, end);
    };

    /**
     * @brief Overwrites from the target all non-zero detail coefficients
     * 
     * The simply overwrite makes sure that no small un-wanted detail coefficients have been created
     * 
     * @param ctx the interpolation context
     */
    void Clear_(const interp_ctx_t* const  ctx) const override{
        // //-------------------------------------------------------------------------
        // // filters, need the forward ones
        // constexpr short_t   ha_lim = len_ha_<TN, TNT> / 2;
        // constexpr short_t   ga_lim = len_ga_<TN, TNT> / 2;
        // const real_t* const ha     = ha_<TN, TNT> + ha_lim;
        // const real_t* const ga     = ga_<TN, TNT> + ga_lim;

        // // get the pointers
        // real_t* const       tdata = ctx->tdata.Write();
        // const real_t* const ddata = ctx->sdata.Read();

        // // go, only tada is changed
        // auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        //     // get 0 if odd, 1 if even (even if negative!!)
        //     const bool odd_x = m_sign(i0) * (i0 % 2);
        //     const bool odd_y = m_sign(i1) * (i1 % 2);
        //     const bool odd_z = m_sign(i2) * (i2 % 2);
        //     m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
        //     m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
        //     m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");

        //     // get the filters: ga for details, ha for scalings
        //     const real_t* const f_x = (odd_x) ? (ga) : (ha);
        //     const real_t* const f_y = (odd_y) ? (ga) : (ha);
        //     const real_t* const f_z = (odd_z) ? (ga) : (ha);

        //     const real_t fact = f_x[0] * f_y[0] * fz[0];

        //     m_assert((i0 >= ctx->srcstart[0]) && (i0 < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0, ctx->srcstart[0], i0, ctx->srcend[0]);
        //     m_assert((i1 >= ctx->srcstart[1]) && (i1 < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1, ctx->srcstart[1], i1, ctx->srcend[1]);
        //     m_assert((i2 >= ctx->srcstart[2]) && (i2 < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2, ctx->srcstart[2], i2, ctx->srcend[2]);

        //     // get the local datassss
        //     real_t* const       ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
        //     const real_t* const lddata = ddata + m_idx(i0, i1, i2, 0, ctx->srcstr);

        //     // let's go tocard
        //     real_t corr = fact * lddata[0];
        //     ltdata[0] -= corr;
        //     m_assert(corr == corr, "the data in %d %d %d is nan %e", i0, i1, i2, corr);

        //     // if (fabs(corr) > 1e-3) {
        //     //     m_log("correction is %e at %d %d %d", corr, i0, i1, i2);
        //     // }
        // };

        // const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        // const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
        // for_loop(&op, start, end);
    };
};

#endif  // SRC_INTERPOLATING_WAVELET_HPP_