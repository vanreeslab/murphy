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
static constexpr short_t len_sd_ = 0;  //!< length of the S dual filter
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
static constexpr real_t sd_[1] = {0.0};

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

//-----------------------------------------------------------------------------
// Wavelet 2.2
template <>
constexpr short_t len_ha_<2, 2> = 5;
template <>
constexpr short_t len_sd_<2, 2> = 2;
template <>
constexpr short_t len_js_<2, 2> = 3;
template <>
constexpr short_t len_ks_<2, 2> = 5;
template <>
constexpr real_t ha_<2, 2>[5] = {-0.125, 0.25, 0.75, 0.25, -0.125};
template <>
constexpr real_t sd_<2, 2>[2] = {0.5, 0.5};
template <>
constexpr real_t js_<2, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
// constexpr real_t js_<2, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
constexpr real_t ks_<2, 2>[5] = {-1.0 / 8.0, 1.0 / 2.0, 3.0 / 4.0, 1.0 / 2.0, -1.0 / 8.0};
// constexpr real_t ks_<2, 2>[5] = {-1.0 / 8.0, 1.0 / 2.0, 3.0 / 4.0, 1.0 / 2.0, -1.0 / 8.0};

//-----------------------------------------------------------------------------
// Wavelet 4.0
template <>
constexpr short_t len_ha_<4, 0> = 1;
template <>
constexpr short_t len_sd_<4, 0> = 4;
template <>
constexpr short_t len_js_<4, 0> = 1;
template <>
constexpr short_t len_ks_<4, 0> = 7;
template <>
constexpr real_t ha_<4, 0>[1] = {1.0};
template <>
constexpr real_t sd_<4, 0>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
template <>
constexpr real_t js_<4, 0>[1] = {0.0};
// constexpr real_t js_<4, 0>[1] = {1.0};
template <>
constexpr real_t ks_<4, 0>[7] = {-1.0 / 16.0, 0.0, 9.0 / 16.0, 0.0, 9.0 / 16.0, 0.0, -1.0 / 16.0};
// constexpr real_t ks_<4, 0>[7] = {-1.0 / 16.0, 0.0, 9.0 / 16.0, 1.0, 9.0 / 16.0, 0.0, -1.0 / 16.0};

//-----------------------------------------------------------------------------
// Wavelet 4.2
template <>
constexpr short_t len_ha_<4, 2> = 9;
template <>
constexpr short_t len_sd_<4, 2> = 4;
template <>
constexpr short_t len_js_<4, 2> = 3;
template <>
constexpr short_t len_ks_<4, 2> = 9;
template <>
constexpr real_t ha_<4, 2>[9] = {0.015625, 0.0, -0.125, 0.25, 0.71875, 0.25, -0.125, 0.0, 0.015625};
template <>
constexpr real_t sd_<4, 2>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
template <>
constexpr real_t js_<4, 2>[3] = {-1.0 / 4.0, 0.0, -1.0 / 4.0};
// constexpr real_t js_<4, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
constexpr real_t ks_<4, 2>[9] = {1.0 / 64.0, -1.0 / 16.0, -1.0 / 8.0, 9.0 / 16.0, -9.0 / 32.0, 9.0 / 16.0, -1.0 / 8.0, -1.0 / 16.0, 1.0 / 64.0};
// constexpr real_t ks_<4, 2>[9] = {1.0 / 64.0, -1.0 / 16.0, -1.0 / 8.0, 9.0 / 16.0, 23.0 / 32.0, 9.0 / 16.0, -1.0 / 8.0, -1.0 / 16.0, 1.0 / 64.0};

//-----------------------------------------------------------------------------
// Wavelet 4.4
template <>
constexpr short_t len_ha_<4, 4> = 13;
template <>
constexpr short_t len_sd_<4, 4> = 4;
template <>
constexpr short_t len_js_<4, 4> = 7;
template <>
constexpr short_t len_ks_<4, 4> = 13;
template <>
constexpr real_t ha_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
template <>
constexpr real_t sd_<4, 4>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};
template <>
constexpr real_t js_<4, 4>[7] = {1.0 / 32.0, 0.0, -9.0 / 32.0, 1.0, -9.0 / 32.0, 0.0, 1.0 / 32.0};
template <>
constexpr real_t ks_<4, 4>[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 16.0, -63.0 / 512.0, 9 / 16.0, 87.0 / 128.0, 9.0 / 16.0, -63.0 / 512.0, -1.0 / 16.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};

//-----------------------------------------------------------------------------
// Wavelet 6.0
template <>
constexpr short_t len_ha_<6, 0> = 1;
template <>
constexpr short_t len_sd_<6, 0> = 6;
template <>
constexpr short_t len_js_<6, 0> = 1;
template <>
constexpr short_t len_ks_<6, 0> = 11;
template <>
constexpr real_t ha_<6, 0>[1] = {1.0};
template <>
constexpr real_t sd_<6, 0>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
template <>
constexpr real_t js_<6, 0>[1] = {1.0};
template <>
constexpr real_t ks_<6, 0>[11] = {3.0 / 256.0, 0.0, -25.0 / 256.0, 0.0, 75.0 / 128.0, 1.0, 75.0 / 128.0, 0.0, -25.0 / 256.0, 0.0, 3.0 / 256.0};

//-----------------------------------------------------------------------------
// Wavelet 6.2
template <>
constexpr short_t len_ha_<6, 2> = 13;
template <>
constexpr short_t len_sd_<6, 2> = 6;
template <>
constexpr short_t len_js_<6, 2> = 3;
template <>
constexpr short_t len_ks_<6, 2> = 13;
template <>
constexpr real_t ha_<6, 2>[13] = {-3.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -125.0 / 1024.0, 1.0 / 4.0, 181.0 / 256.0, 1.0 / 4.0, -125.0 / 1024.0, 0.0, 11.0 / 512.0, 0.0, -3.0 / 1024.0};
template <>
constexpr real_t sd_<6, 2>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
template <>
constexpr real_t js_<6, 2>[3] = {-1.0 / 4.0, 1.0, -1.0 / 4.0};
template <>
constexpr real_t ks_<6, 2>[13] = {-3.0 / 1024.0, 3.0 / 256.0, 11.0 / 512.0, -25.0 / 256.0, -125.0 / 1024.0, 75.0 / 128.0, 181.0 / 256.0, 75.0 / 128.0, -125.0 / 1024.0, -25.0 / 256.0, 11.0 / 512.0, 3.0 / 256.0, -3.0 / 1024.0};

//-----------------------------------------------------------------------------
// Wavelet 6.4
template <>
constexpr short_t len_ha_<6, 4> = 17;
template <>
constexpr short_t len_sd_<6, 4> = 6;
template <>
constexpr short_t len_js_<6, 4> = 7;
template <>
constexpr short_t len_ks_<6, 4> = 17;
template <>
constexpr real_t ha_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
template <>
constexpr real_t sd_<6, 4>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};
template <>
constexpr real_t js_<6, 4>[7] = {1.0 / 32.0, 0.0, -9.0 / 32.0, 1.0, -9.0 / 32.0, 0.0, 1.0 / 32.0};
template <>
constexpr real_t ks_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 3.0 / 256.0, 87.0 / 2048.0, -25.0 / 256.0, -243.0 / 2048.0, 75.0 / 128.0, 2721.0 / 4096.0, 75.0 / 128.0, -243.0 / 2048.0, -25.0 / 256.0, 87.0 / 2048.0, 3.0 / 256.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};

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
    const short_t len_sd() const override { return len_sd_<TN, TNT>; };

   protected:
    /**
     * @brief coarsen the values of the source memory to gather them in the target memory.
     * 
     * The details are destroyed whatever value they have. we use the ga filter
     * 
     * @param ctx the interpolation context
     */
    void Coarsen_(const m_ptr<const interp_ctx_t>& ctx) const override {
        //-------------------------------------------------------------------------
        // the size is know at compilation
        constexpr short_t ha_lim = len_ha_<TN, TNT> / 2;

        const real_t        alpha = ctx->alpha;
        const real_t* const ha    = ha_<TN, TNT> + ha_lim;

        // restrict the pointers
        real_t*       tdata = ctx->tdata.Write();
        const real_t* sdata = ctx->sdata.Read();

        // assume alignment on the source
        // // m_assume_aligned(ctx->sdata);
        const bidx_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const bidx_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
            // do some checks
            m_assert(((2 * i0 - ha_lim) >= (ctx->srcstart[0])) && ((2 * i0 + ha_lim) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", 2 * i0 - ha_lim, ctx->srcstart[0], 2 * i0 + ha_lim, ctx->srcend[0]);
            m_assert(((2 * i1 - ha_lim) >= (ctx->srcstart[1])) && ((2 * i1 + ha_lim) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", 2 * i1 - ha_lim, ctx->srcstart[1], 2 * i1 + ha_lim, ctx->srcend[1]);
            m_assert(((2 * i2 - ha_lim) >= (ctx->srcstart[2])) && ((2 * i2 + ha_lim) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", 2 * i2 - ha_lim, ctx->srcstart[2], 2 * i2 + ha_lim, ctx->srcend[2]);

            //get the local adress of the source, the target and the constant
            real_t*       ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
            const real_t* lsdata = sdata + m_idx(2 * i0, 2 * i1, 2 * i2, 0, ctx->srcstr);

            // add the constant
            ltdata[0] = 0.0;

            // apply the filter
            for (short_t id2 = -ha_lim; id2 <= ha_lim; ++id2) {
                for (short_t id1 = -ha_lim; id1 <= ha_lim; ++id1) {
                    for (short_t id0 = -ha_lim; id0 <= ha_lim; ++id0) {
                        // add the info
                        const real_t fact = ha[id0] * ha[id1] * ha[id2];
                        ltdata[0] += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];

                        // check for nan's
                        m_assert(lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)] == lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)], "nan detected");
                        m_assert(fact == fact, "nan detected");
                    }
                }
            }
            // check for nan's
            m_assert(ltdata[0] == ltdata[0], "nan detected");
        };

        for_loop(&op, start, end);
        //-------------------------------------------------------------------------
    }

    /**
     * @brief refine the source memory to get the associated target memory information
     * 
     * @warning we assume that the detail coefficients are null, and we use the dual filter
     * Hence, the values of the function are the scaling coefficient and we simply apply the dual-lifting scheme to obtain the missing information
     *
     * @param ctx the interpolation context
     */
    void Refine_(const m_ptr<const interp_ctx_t>& ctx) const override {
        m_assert(ctx->alpha == 0.0, "the alpha = %e must be = 0.0", ctx->alpha);
        //-------------------------------------------------------------------------
        // the gs_lim is know @ compile time
        constexpr short_t sd_lim = len_sd_<TN, TNT> / 2 - 1;
        // constexpr short_t ks_lim = len_ks_<TN, TNT> / 2 - 1;

        const real_t one      = 1.0;
        const real_t alpha    = ctx->alpha;
        const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        // const real_t* ks = ks_<TN, TNT> + ks_lim;

        // const short  ga_lim = 1;
        // const real_t ga[3]  = {-0.5, 1.0, -0.5};
        // const short  ha_lim = 2;
        // const real_t ha[5]  = {-1.0 / 8.0, 1.0 / 4.0, 3.0 / 4.0, 1.0 / 4.0, -1.0 / 8.0};


        // restrict the pointers
        real_t*       tdata = ctx->tdata.Write();
        const real_t* sdata = ctx->sdata.Read();

        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const short_t ix = m_sign(i0) * (i0 % 2);
            const short_t iy = m_sign(i1) * (i1 % 2);
            const short_t iz = m_sign(i2) * (i2 % 2);
            m_assert(ix == 0 || ix == 1, "this are the two possible values");
            m_assert(iy == 0 || iy == 1, "this are the two possible values");
            m_assert(iz == 0 || iz == 1, "this are the two possible values");

            // get the target location
            real_t* ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);

            //get the local adress of the source, a bit more complicated to handle the negative numbers
            const lid_t   i0_s   = (i0 - ix) / 2;
            const lid_t   i1_s   = (i1 - iy) / 2;
            const lid_t   i2_s   = (i2 - iz) / 2;
            const real_t* lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);
            m_assert((i0_s * 2) <= i0, "if not, we made something wrong...: source = %d, target = %d", i0_s, i0);
            m_assert((i1_s * 2) <= i1, "if not, we made something wrong...: source = %d, target = %d", i1_s, i1);
            m_assert((i2_s * 2) <= i2, "if not, we made something wrong...: source = %d, target = %d", i2_s, i2);

            // get the filter, depending on if I am odd or even
            const real_t* sd_x         = (ix == 1) ? (sd) : (&one);
            const real_t* sd_y         = (iy == 1) ? (sd) : (&one);
            const real_t* sd_z         = (iz == 1) ? (sd) : (&one);
            const bidx_t  lim_start[3] = {sd_lim * ix, sd_lim * iy, sd_lim * iz};
            const bidx_t  lim_end[3]   = {(sd_lim + 1) * ix, (sd_lim + 1) * iy, (sd_lim + 1) * iz};

            m_assert(((i0 / 2 - lim_start[0]) >= ctx->srcstart[0]) && ((i0 / 2 + lim_end[0]) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0 - sd_lim, ctx->srcstart[0], i0 + sd_lim, ctx->srcend[0]);
            m_assert(((i1 / 2 - lim_start[1]) >= ctx->srcstart[1]) && ((i1 / 2 + lim_end[1]) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1 - sd_lim, ctx->srcstart[1], i1 + sd_lim, ctx->srcend[1]);
            m_assert(((i2 / 2 - lim_start[2]) >= ctx->srcstart[2]) && ((i2 / 2 + lim_end[2]) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2 - sd_lim, ctx->srcstart[2], i2 + sd_lim, ctx->srcend[2]);

            // add the constant array
            ltdata[0] = 0.0;

            // if one dim is even, id = 0, -> gs[0] = 1 and that's it
            // if one dim is odd, id = 1, -> we loop on gs, business as usual
            for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                    for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                        const real_t fact = sd_x[id0] * sd_y[id1] * sd_z[id2];
                        ltdata[0] += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];
                    }
                }
            }
            m_assert(ltdata[0] == ltdata[0], "the value cannot be nan: block @ %d %d %d: %f", i0, i1, i2, ltdata[0]);
        };

        for_loop(&op, start, end);
        //-------------------------------------------------------------------------
    };

    /**
     * @brief gets the maximum infinite norm of the detail coefficients over the block, no change to the block memory is performed
     * 
     * we use the Sd filter as a substitute to the Ga filter (it's cheaper in memory and in computations)
     * 
     * if the ctx->tdata is empty, we return the max detail
     * if the ctx->tdata is not empty, we additionally store it's value in the field if the abs(detail) < tol
     * 
     * @param ctx the interpolation context
     * @param details_max the maximum of the local detail coefficients
     */
    void Detail_(const m_ptr<const interp_ctx_t>& ctx, const m_ptr<real_t>& details_max) const override {
        m_assert(*details_max == 0.0, "the value must be 0.0");  // assert we can max on it
        //-------------------------------------------------------------------------
        // the size is know @ compiler time
        constexpr short_t sd_lim = (len_sd_<TN, TNT> / 2 - 1);

        const real_t one      = 1.0;
        const real_t zero      = 0.0;
        const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        const real_t*       sdata = ctx->sdata.Read();
        const real_t* const sd    = sd_<TN, TNT> + sd_lim;

        // check if we need to store some value
        bool    store = !(ctx->tdata.IsEmpty());
        real_t  temp  = 0.0;
        real_t* tdata = (store) ? (ctx->tdata.Write()) : (&temp);

        const short  ga_lim = 1;
        const real_t ga[3]  = {-0.5, 1.0, -0.5};
        const short  ha_lim = 2;
        const real_t ha[5]  = {-1.0 / 8.0, 1.0 / 4.0, 3.0 / 4.0, 1.0 / 4.0, -1.0 / 8.0};

        // m_log("we store? %d",store);

        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const lda_t ix = m_sign(i0) * (i0 % 2);
            const lda_t iy = m_sign(i1) * (i1 % 2);
            const lda_t iz = m_sign(i2) * (i2 % 2);
            m_assert(ix == 0 || ix == 1, "this are the two possible values");
            m_assert(iy == 0 || iy == 1, "this are the two possible values");
            m_assert(iz == 0 || iz == 1, "this are the two possible values");

            const lid_t   i0_s   = (i0 - ix);
            const lid_t   i1_s   = (i1 - iy);
            const lid_t   i2_s   = (i2 - iz);
            const real_t* lsdata = sdata + m_idx(i0, i1, i2, 0, ctx->srcstr);

            const bool is_scaling = (ix == 0) && (iy == 0) && (iz == 0);

            // get the filter, depending on if I am odd or even
            // const real_t* const sd_x         = is_scaling ? (&zero) : ((ix == 1) ? (ga + ga_lim) : (&one));
            // const real_t* const sd_y         = is_scaling ? (&zero) : ((iy == 1) ? (ga + ga_lim) : (&one));
            // const real_t* const sd_z         = is_scaling ? (&zero) : ((iz == 1) ? (ga + ga_lim) : (&one));
            // const bidx_t        lim_start[3] = {(ga_lim)*ix, (ga_lim)*iy, (ga_lim)*iz};
            // const bidx_t        lim_end[3]   = {(ga_lim)*ix, (ga_lim)*iy, (ga_lim)*iz};

            const real_t* const sd_x         = (ix == 1) ? (ga + ga_lim) : (ha + ha_lim);
            const real_t* const sd_y         = (iy == 1) ? (ga + ga_lim) : (ha + ha_lim);
            const real_t* const sd_z         = (iz == 1) ? (ga + ga_lim) : (ha + ha_lim);
            const bidx_t        lim_start[3] = {
                (ga_lim) * (ix == 1) + (ha_lim) * (ix == 0),
                (ga_lim) * (iy == 1) + (ha_lim) * (iy == 0),
                (ga_lim) * (iz == 1) + (ha_lim) * (iz == 0)};
            const bidx_t lim_end[3] = {
                (ga_lim) * (ix == 1) + (ha_lim) * (ix == 0),
                (ga_lim) * (iy == 1) + (ha_lim) * (iy == 0),
                (ga_lim) * (iz == 1) + (ha_lim) * (iz == 0)};

            // if one dim is even, id = 0, -> gs[0] = 1 and that's it
            // if one dim is odd, id = 1, -> we loop on gs, business as usual
            real_t interp = 0.0;
            for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                    for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                        // const real_t fact = sd_x[id0] * sd_y[id1] * sd_z[id2];
                        const real_t fact = sd_x[id0] * sd_y[id1] * sd_z[id2];
                        interp += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];

                        // if (i1 == 4 & i2 == 4) {
                        //     printf("we add %e * %e \n", i0, fact, lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)]);
                        // }
                    }
                }
            }
            real_t detail = interp * (!is_scaling);  // sdata[m_idx(i0, i1, i2, 0, ctx->srcstr)] - interp;

            // if (i1 == 4 & i2 == 4) {
            //     printf("detail in %d = %e\n", i0, detail);
            //     // m_assert(false,"coucou");
            //     m_log("limits where %d %d %d",lim_start[0],lim_start[1],lim_start[2]);
            // }

            // check the maximum
            (*details_max) = m_max(fabs(detail), (*details_max));

            // store if needed, the index is 0 if not store
            tdata[store * m_idx(i0, i1, i2, 0, ctx->trgstr)] = detail * (fabs(detail) < ctx->alpha);
            m_assert(tdata[store * m_idx(i0, i1, i2, 0, ctx->trgstr)] == tdata[store * m_idx(i0, i1, i2, 0, ctx->trgstr)],"the data in %d %d %d is nan %e",i0,i1,i2,detail * (fabs(detail) < ctx->alpha));
        };

        // run that
        *details_max = 0.0;
        for_loop(&op, start, end);
    };

    /**
     * @brief Use the computed details to discard the values
     * 
     * 
     * 
     * @param ctx the interpolation context
     * @param details_max return the max of the detail on the current block
     
     */
    void Smooth_(const m_ptr<const interp_ctx_t>& ctx) const override {
        //-------------------------------------------------------------------------
        const real_t one      = 1.0;
        const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        // const real_t* sdata    = ctx->sdata.Read();
        const real_t* ddata = ctx->sdata.Read();
        real_t* tdata = ctx->tdata.Write();

        // const real_t one = 1.0;

        // filter characteristics
        constexpr short_t js_lim = len_js_<TN, TNT> / 2;
        constexpr short_t ks_lim = len_ks_<TN, TNT> / 2;
        const real_t*     js     = js_<TN, TNT> + js_lim;
        const real_t*     ks     = ks_<TN, TNT> + ks_lim;

        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const bool odd_x = m_sign(i0) * (i0 % 2);
            const bool odd_y = m_sign(i1) * (i1 % 2);
            const bool odd_z = m_sign(i2) * (i2 % 2);
            m_assert(odd_x == 0 || odd_x == 1, "this are the two possible values");
            m_assert(odd_y == 0 || odd_y == 1, "this are the two possible values");
            m_assert(odd_z == 0 || odd_z == 1, "this are the two possible values");
            // const lid_t i0_s = (i0 - ix);
            // const lid_t i1_s = (i1 - iy);
            // const lid_t i2_s = (i2 - iz);

            const bool is_scaling = (!odd_x && !odd_y && !odd_z);

            // get the filter, depending on if I am odd or even
            // const real_t* const f_x          = (is_scaling)? (js) :( (odd_x) ? (ks) : (&one));
            // const real_t* const f_y          = (is_scaling)? (js) :( (odd_y) ? (ks) : (&one));
            // const real_t* const f_z          = (is_scaling)? (js) :( (odd_z) ? (ks) : (&one));
            const real_t* f_x    = (odd_x) ? (ks) : (js);
            const real_t* f_y    = (odd_y) ? (ks) : (js);
            const real_t* f_z    = (odd_z) ? (ks) : (js);
            const bidx_t  lim[3] = {
                (js_lim) * (!odd_x) + (ks_lim) * (odd_x),
                (js_lim) * (!odd_y) + (ks_lim) * (odd_y),
                (js_lim) * (!odd_z) + (ks_lim) * (odd_z)};
            // const bidx_t        lim[3] = {
            //     (is_scaling)? (js_lim) :( (odd_x) ? (ks_lim) : (0)),
            //     (is_scaling)? (js_lim) :( (odd_y) ? (ks_lim) : (0)),
            //     (is_scaling)? (js_lim) :( (odd_z) ? (ks_lim) : (0))};

            m_assert(((i0 - lim[0]) >= ctx->srcstart[0]) && ((i0 + lim[0]) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0 - lim[0], ctx->srcstart[0], i0 + lim[0], ctx->srcend[0]);
            m_assert(((i1 - lim[1]) >= ctx->srcstart[1]) && ((i1 + lim[1]) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1 - lim[1], ctx->srcstart[1], i1 + lim[1], ctx->srcend[1]);
            m_assert(((i2 - lim[2]) >= ctx->srcstart[2]) && ((i2 + lim[2]) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2 - lim[2], ctx->srcstart[2], i2 + lim[2], ctx->srcend[2]);

            // apply the filter
            real_t*       ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
            const real_t* lddata = ddata + m_idx(i0, i1, i2, 0, ctx->srcstr);

            // if I am odd in x or y or z, I need to add the detail value
            real_t corr = 0.0; //(odd_x || odd_y || odd_z) * lddata[0];

            // if (i1 == 5 & i2 == 4) {
            //     printf("\n================================= \n");
            //     real_t fact = (odd_x || odd_y || odd_z);
            //     printf("%e -> corr = %e \n", fact, corr);
            // }
            // corr = 0.0;

            // real_t sum = 0.0;

            for (bidx_t id2 = -lim[2]; id2 <= lim[2]; ++id2) {
                // if (i1 == 5 & i2 == 4)  {
                //     printf("\n---------- corr = %e\n", corr);
                // }
                for (bidx_t id1 = -lim[1]; id1 <= lim[1]; ++id1) {
                    // if (i1 == 5 & i2 == 4) {
                    //     printf("\n");
                    // }
                    for (bidx_t id0 = -lim[0]; id0 <= lim[0]; ++id0) {
                        // m_log("reading %d %d %d",i0+id0,i1+id1,i2+id2);
                        // const real_t fact_x = f_x[id0] + (id0 == 0) * (odd_x);
                        // const real_t fact_y = f_y[id1] + (id1 == 0) * (odd_y);
                        // const real_t fact_z = f_z[id2] + (id2 == 0) * (odd_z);
                        // const real_t fact   = fact_x * fact_y * fact_z;

                        m_assert(!((id0 + i0) % 2 == 0 && (id1 + i1) % 2 == 0 && (id2 + i2) % 2 == 0 && fabs(lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]) > 1e-16), "detail is wrong: %e", lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]);
                        // m_assert(std::fabs(lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]) <= 1e-2, "the details is too big");

                        // const real_t dx   = (id0 == 0) && (odd_x);
                        // const real_t dy   = (id1 == 0) && (odd_y);
                        // const real_t dz   = (id2 == 0) && (odd_z);
                        // const real_t sign = m_min(m_sign(f_x[id0]), m_min(m_sign(f_y[id1]), m_sign(f_z[id2])));
                        // const real_t fact = fabs(f_x[id0]) * fabs(f_y[id1]) * fabs(f_z[id2]);
                        const real_t fact = (f_x[id0]) * (f_y[id1]) * (f_z[id2]);

                        // // if(i0==10 && i1==10 & i2==10){ printf(" %e * %e ,",fact,lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]);}
                        // if (i1 == 5 & i2 == 4) {
                        //     // printf(" %e ,", fact);
                        //     // printf(" %e * %e ,", sign* fact, lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]);
                        //     printf(" %e * %e ,", fact, lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]);
                        // }

                        // sum += ((i0 % 2 + i1 % 2 + i2 % 2) != 0) * fact;

                        // sum += fact;

                        // if(i0==5 && i1==5 & i2==15){
                        //     m_log("we are odd? %d %d %d -> in pos %d %d %d, filter value = %e * detail = %e",odd_x,odd_y,odd_z,id0,id1,id2,fact,lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)]);
                        // }

                        // m_log("we are odd? %d %d %d -> in pos %d %d %d, filter value = %e",odd_x,odd_y,odd_z,id0,id1,id2,fact);

                        // const real_t fact = f_x[id0] * f_y[id1] * f_z[id2] + (id0 == 0) * (id1 == 0) * (id2 == 0);
                        // const real_t fact = dx * dy * dz + f_x[id0] * dy * dz + dx * f_y[id1] * dz + dx * dy * f_z[id2] +
                        // f_x[id0] * f_y[id1] * dz + dx * f_y[id1] * f_z[id2] + f_x[id0] * dy * f_z[id2] + f_x[id0] * f_y[id1] * f_z[id2];

                        //  f_x[id0] * f_y[id1] * f_z[id2] + (id0 == 0) * (id1 == 0) * (id2 == 0);
                        // m_assert(fact == 1.0 || fact == 0.0, "fact = %e is not ok in %d %d %d for %d %d %d", fact, id0, id1, id2, i0, i1, i2);
                        corr += fact * lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)];

                        m_assert(lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)] == lddata[m_idx(id0, id1, id2, 0, ctx->srcstr)], "no nan is allowed: %d %d %d with sub = %d %d %d", i0, i1, i2, id0, id1, id2);
                        m_assert(ltdata[0] == ltdata[0], "no nan is allowed: %d %d %d with sub = %d %d %d", i0, i1, i2, id0, id1, id2);
                    }
                }
            }
            // if (i1 == 5 & i2 == 4) {
            //     m_log("\nfinal corr in %d %d %d  = %e", i0, i1, i2, corr);
            // }
            // // m_assert(sum == 1.0,"the sum MUST be one");
            // if (i0 % 2 == 0 && i1 % 2 == 0 && i2 == 0) {
            // m_assert(std::fabs(corr) <= 1e-2, "the correction performed in %d %d %d = %e is too big, detail was %e, value %e", i0, i1, i2, corr, lddata[0], ltdata[0]);
            ltdata[0] = ltdata[0] - corr;
            // }
            // if ((i0 % 2 + i1 % 2 + i2 % 2) == 1) {
            //     m_assert(std::fabs(corr) <= 1e-2, "the correction performed in %d %d %d = %e is too big, detail was %e, value %e", i0, i1, i2, corr, lddata[0], ltdata[0]);
            //     ltdata[0] = ltdata[0] - corr;
            // }

            // m_assert(true, "coucou");
        };
        for_loop(&op, start, end);
    };

        //     /**
        //      * @brief Compute the details using the source field and write them in the target field
        //      *
        //      * @param ctx the interpolation context
        //      */
        //     void WriteDetail_(m_ptr<const interp_ctx_t> ctx) const override {
        //         m_assert(ctx->srcstart[0] == 0 && ctx->srcend[0] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->srcstart[0], ctx->srcend[0]);
        //         m_assert(ctx->srcstart[1] == 0 && ctx->srcend[1] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->srcstart[1], ctx->srcend[1]);
        //         m_assert(ctx->srcstart[2] == 0 && ctx->srcend[2] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->srcstart[2], ctx->srcend[2]);
        //         m_assert(ctx->trgstart[0] == 0 && ctx->trgend[0] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->trgstart[0], ctx->trgend[0]);
        //         m_assert(ctx->trgstart[1] == 0 && ctx->trgend[1] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->trgstart[1], ctx->trgend[1]);
        //         m_assert(ctx->trgstart[2] == 0 && ctx->trgend[2] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->trgstart[2], ctx->trgend[2]);
        //         //-------------------------------------------------------------------------
        //         // the size is know @ compiler time
        //         constexpr short_t gs_lim = (len_gs_<TN, TNT> / 2 - 1);

        //         const real_t        one = 1.0;
        //         const real_t* const gs  = gs_<TN, TNT> + gs_lim;

        //         // get the target + criterion field
        //         const real_t* sdata = ctx->sdata.Read();
        //         real_t*       tdata = ctx->tdata.Write();

        //         // declare the lambda to run
        //         auto lambda = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        //             // get 0 if odd, 1 if even (even if negative!!)
        //             const lda_t iy = m_sign(i1) * (i1 % 2);
        //             const lda_t ix = m_sign(i0) * (i0 % 2);
        //             const lda_t iz = m_sign(i2) * (i2 % 2);
        //             m_assert(ix == 0 || ix == 1, "this are the two possible values");
        //             m_assert(iy == 0 || iy == 1, "this are the two possible values");
        //             m_assert(iz == 0 || iz == 1, "this are the two possible values");

        //             const lid_t   i0_s   = (i0 - ix);
        //             const lid_t   i1_s   = (i1 - iy);
        //             const lid_t   i2_s   = (i2 - iz);
        //             const real_t* lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);

        //             // get the filter, depending on if I am odd or even
        //             const real_t* const gs_x         = (ix == 1) ? (gs) : (&one);
        //             const real_t* const gs_y         = (iy == 1) ? (gs) : (&one);
        //             const real_t* const gs_z         = (iz == 1) ? (gs) : (&one);
        //             const bidx_t        lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
        //             const bidx_t        lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

        //             // if one dim is even, id = 0, -> gs[0] = 1 and that's it
        //             // if one dim is odd, id = 1, -> we loop on gs, business as usual
        //             real_t interp = 0.0;
        //             for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
        //                 for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
        //                     for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
        //                         const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
        //                         interp += fact * lsdata[m_idx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->srcstr)];
        //                     }
        //                 }
        //             }
        //             real_t* ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
        //             ltdata[0]      = sdata[m_idx(i0, i1, i2, 0, ctx->srcstr)] - interp;

        //             // check that we retrieve the original value if we are a scaling coef
        //             m_assert(!(ix == 0 && iy == 0 && iz == 0 && tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)] != 0.0), "the target value should be 0.0 instead of %e", tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)]);
        //         };

        //         // do the loop
        //         for_loop<0, M_N>(&lambda);
        //         //-------------------------------------------------------------------------
        //     };
        // };
};

#endif  // SRC_INTERPOLATING_WAVELET_HPP_