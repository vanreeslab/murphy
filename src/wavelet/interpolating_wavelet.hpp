#ifndef SRC_INTERPOLATING_WAVELET_HPP_
#define SRC_INTERPOLATING_WAVELET_HPP_

#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "gridblock.hpp"
#include "wavelet.hpp"

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
template <>
constexpr real_t ha_<4, 0>[1] = {1.0};
template <>
constexpr real_t gs_<4, 0>[4] = {-0.0625, 0.5625, 0.5625, -0.0625};

//-----------------------------------------------------------------------------ai un
// Wavelet 4.2
template <>
constexpr sid_t len_ha_<4, 2> = 9;
template <>
constexpr sid_t len_gs_<4, 2> = 4;
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
template <>
constexpr real_t ha_<6, 4>[17] = {3.0 / 8192.0, 0.0, -13.0 / 2048.0, 0.0, 87.0 / 2048.0, -1.0 / 32.0, -243.0 / 2048.0, 9.0 / 32.0, 2721.0 / 4096.0, 9.0 / 32.0, -243.0 / 2048.0, -1.0 / 32.0, 87.0 / 2048.0, 0.0, -13.0 / 2048.0, 0.0, 3.0 / 8192.0};
template <>
constexpr real_t gs_<6, 4>[6] = {3.0 / 256.0, -25.0 / 256.0, 75.0 / 128.0, 75.0 / 128.0, -25.0 / 256.0, 3.0 / 256.0};

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
template <lda_t TN = M_WAVELET_N, lda_t TNT = M_WAVELET_NT>
class InterpolatingWavelet : public Wavelet {
   public:
    const sid_t N() const override { return TN; }
    const sid_t Nt() const override { return TNT; }
    const sid_t len_ha() const override { return len_ha_<TN, TNT>; };
    const sid_t len_gs() const override { return len_gs_<TN, TNT>; };

   protected:
    /**
     * @brief coarsen the values of the source memory to gather them in the target memory.
     * 
     * @param ctx the interpolation context
     */
    void Coarsen_(m_ptr<const interp_ctx_t> ctx) const override {
        //-------------------------------------------------------------------------
        // the size is know at compilation
        constexpr sid_t ha_lim = len_ha_<TN, TNT> / 2;

        // assure alignment for the target, the source, the constant and the temp data
        const real_t        alpha = ctx->alpha;
        const real_t* const ha    = ha_<TN, TNT> + ha_lim;

        // restrict the pointers
        real_t*       tdata = ctx->tdata.Write();
        const real_t* cdata = ctx->cdata.Read();
        const real_t* sdata = ctx->sdata.Read();

        // assume alignment on the source
        // // m_assume_aligned(ctx->sdata);
        const bidx_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const bidx_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
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
            for (sid_t id2 = -ha_lim; id2 <= ha_lim; ++id2) {
                for (sid_t id1 = -ha_lim; id1 <= ha_lim; ++id1) {
                    for (sid_t id0 = -ha_lim; id0 <= ha_lim; ++id0) {
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

        // for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ++ik2) {
        //     for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ++ik1) {
        //         for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ++ik0) {
        //             //get the local adress of the source, the target and the constant
        //             data_ptr ltdata = tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
        //             m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
        //             // const data_ptr lcdata = cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
        //             const data_ptr lsdata = sdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);

        //             // do some checks
        //             m_assert(((2 * ik0 - ha_lim) >= (ctx->srcstart[0])) && ((2 * ik0 + ha_lim) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", 2 * ik0 - ha_lim, ctx->srcstart[0], 2 * ik0 + ha_lim, ctx->srcend[0]);
        //             m_assert(((2 * ik1 - ha_lim) >= (ctx->srcstart[1])) && ((2 * ik1 + ha_lim) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", 2 * ik1 - ha_lim, ctx->srcstart[1], 2 * ik1 + ha_lim, ctx->srcend[1]);
        //             m_assert(((2 * ik2 - ha_lim) >= (ctx->srcstart[2])) && ((2 * ik2 + ha_lim) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", 2 * ik2 - ha_lim, ctx->srcstart[2], 2 * ik2 + ha_lim, ctx->srcend[2]);

        //             // add the constant
        //             // ltdata[0] = alpha * lcdata[0];
        //             m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
        //             ltdata[0] = 0.0;

        //             // apply the filter
        //             for (sid_t id2 = -ha_lim; id2 <= ha_lim; ++id2) {
        //                 for (sid_t id1 = -ha_lim; id1 <= ha_lim; ++id1) {
        //                     for (sid_t id0 = -ha_lim; id0 <= ha_lim; ++id0) {
        //                         // add the info
        //                         const real_t fact = ha[id0] * ha[id1] * ha[id2];
        //                         ltdata[0] += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];

        //                         // check for nan's
        //                         m_assert(lsdata[0] == lsdata[0], "nan detected");
        //                         m_assert(fact == fact, "nan detected");
        //                     }
        //                 }
        //             }
        //             m_assert(ltdata[0] == ltdata[0], "nan detected");
        //         }
        //     }
        // }
        //-------------------------------------------------------------------------
    }

    /**
     * @brief refine the source memory to get the associated target memory information
     * 
     * @warning we assume that the detail coefficients are null.
     * Hence, the values of the function are the scaling coefficient and we simply apply the dual-lifting scheme to obtain the missing information
     *
     * @param ctx the interpolation context
     */
    void Refine_(m_ptr<const interp_ctx_t> ctx) const override {
        //-------------------------------------------------------------------------
        // the gs_lim is know @ compile time
        constexpr sid_t gs_lim = len_gs_<TN, TNT> / 2 - 1;

        const real_t one      = 1.0;
        const real_t alpha    = ctx->alpha;
        const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        const real_t* gs = gs_<TN, TNT> + gs_lim;

        // restrict the pointers
        real_t*       tdata = ctx->tdata.Write();
        const real_t* cdata = ctx->cdata.Read();
        const real_t* sdata = ctx->sdata.Read();

        auto op = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const sid_t ix = m_sign(i0) * (i0 % 2);
            const sid_t iy = m_sign(i1) * (i1 % 2);
            const sid_t iz = m_sign(i2) * (i2 % 2);
            m_assert(ix == 0 || ix == 1, "this are the two possible values");
            m_assert(iy == 0 || iy == 1, "this are the two possible values");
            m_assert(iz == 0 || iz == 1, "this are the two possible values");

            // get the target location
            real_t* ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
            // const data_ptr lcdata = cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
            m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");

            //get the local adress of the source, a bit more complicated to handle the negative numbers
            const lid_t   i0_s   = (i0 - ix) / 2;
            const lid_t   i1_s   = (i1 - iy) / 2;
            const lid_t   i2_s   = (i2 - iz) / 2;
            const real_t* lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);
            m_assert((i0_s * 2) <= i0, "if not, we made something wrong...: source = %d, target = %d", i0_s, i0);
            m_assert((i1_s * 2) <= i1, "if not, we made something wrong...: source = %d, target = %d", i1_s, i1);
            m_assert((i2_s * 2) <= i2, "if not, we made something wrong...: source = %d, target = %d", i2_s, i2);

            // get the filter, depending on if I am odd or even
            const real_t* gs_x         = (ix == 1) ? (gs) : (&one);
            const real_t* gs_y         = (iy == 1) ? (gs) : (&one);
            const real_t* gs_z         = (iz == 1) ? (gs) : (&one);
            const bidx_t  lim_start[3] = {gs_lim * ix, gs_lim * iy, gs_lim * iz};
            const bidx_t  lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

            m_assert(((i0 / 2 - lim_start[0]) >= ctx->srcstart[0]) && ((i0 / 2 + lim_end[0]) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", i0 - gs_lim, ctx->srcstart[0], i0 + gs_lim, ctx->srcend[0]);
            m_assert(((i1 / 2 - lim_start[1]) >= ctx->srcstart[1]) && ((i1 / 2 + lim_end[1]) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", i1 - gs_lim, ctx->srcstart[1], i1 + gs_lim, ctx->srcend[1]);
            m_assert(((i2 / 2 - lim_start[2]) >= ctx->srcstart[2]) && ((i2 / 2 + lim_end[2]) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", i2 - gs_lim, ctx->srcstart[2], i2 + gs_lim, ctx->srcend[2]);

            // add the constant array
            m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
            ltdata[0] = 0.0;  //alpha * lcdata[0];

            // if one dim is even, id = 0, -> gs[0] = 1 and that's it
            // if one dim is odd, id = 1, -> we loop on gs, business as usual
            for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                    for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                        const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                        ltdata[0] += fact * lsdata[m_idx(id0, id1, id2, 0, ctx->srcstr)];
                    }
                }
            }
            m_assert(ltdata[0] == ltdata[0], "the value cannot be nan: block @ %d %d %d: %f", i0, i1, i2, ltdata[0]);
        };

        for_loop(&op, start, end);

        // // for each of the data for the needed target
        // for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
        //     for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
        //         for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
        //             // get 0 if odd, 1 if even (even if negative!!)
        //             const sid_t iy = m_sign(ik1) * (ik1 % 2);
        //             const sid_t ix = m_sign(ik0) * (ik0 % 2);
        //             const sid_t iz = m_sign(ik2) * (ik2 % 2);
        //             m_assert(ix == 0 || ix == 1, "this are the two possible values");
        //             m_assert(iy == 0 || iy == 1, "this are the two possible values");
        //             m_assert(iz == 0 || iz == 1, "this are the two possible values");

        //             // get the target location
        //             data_ptr ltdata = tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
        //             m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
        //             // const data_ptr lcdata = cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);

        //             //get the local adress of the source, a bit more complicated to handle the negative numbers
        //             const lid_t    ik0_s  = (ik0 - ix) / 2;
        //             const lid_t    ik1_s  = (ik1 - iy) / 2;
        //             const lid_t    ik2_s  = (ik2 - iz) / 2;
        //             const data_ptr lsdata = sdata + m_sidx(ik0_s, ik1_s, ik2_s, 0, ctx->srcstr);
        //             m_assert((ik0_s * 2) <= ik0, "if not, we made something wrong...: source = %d, target = %d", ik0_s, ik0);
        //             m_assert((ik1_s * 2) <= ik1, "if not, we made something wrong...: source = %d, target = %d", ik1_s, ik1);
        //             m_assert((ik2_s * 2) <= ik2, "if not, we made something wrong...: source = %d, target = %d", ik2_s, ik2);

        //             // get the filter, depending on if I am odd or even
        //             const_mem_ptr gs_x         = (ix == 1) ? (gs) : (&one);
        //             const_mem_ptr gs_y         = (iy == 1) ? (gs) : (&one);
        //             const_mem_ptr gs_z         = (iz == 1) ? (gs) : (&one);
        //             const sid_t   lim_start[3] = {gs_lim * ix, gs_lim * iy, gs_lim * iz};
        //             const sid_t   lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

        //             m_assert(((ik0 / 2 - lim_start[0]) >= ctx->srcstart[0]) && ((ik0 / 2 + lim_end[0]) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
        //             m_assert(((ik1 / 2 - lim_start[1]) >= ctx->srcstart[1]) && ((ik1 / 2 + lim_end[1]) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
        //             m_assert(((ik2 / 2 - lim_start[2]) >= ctx->srcstart[2]) && ((ik2 / 2 + lim_end[2]) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);

        //             // add the constant array
        //             m_assert(cdata == nullptr, "the constant data must be nullptr for the moment");
        //             ltdata[0] = 0.0;  //alpha * lcdata[0];

        //             // if one dim is even, id = 0, -> gs[0] = 1 and that's it
        //             // if one dim is odd, id = 1, -> we loop on gs, business as usual
        //             for (sid_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
        //                 for (sid_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
        //                     for (sid_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
        //                         const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
        //                         ltdata[0] += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
        //                     }
        //                 }
        //             }
        //             m_assert(ltdata[0] == ltdata[0], "the value cannot be nan: block @ %d %d %d: %f", ik0, ik1, ik2, ltdata[0]);
        //         }
        //     }
        // }
        //-------------------------------------------------------------------------
    };

    /**
     * @brief gets the maximum infinite norm of the detail coefficients over the block
     * 
     * we use the gs filter, which is nothing but the dual lifting step. We reconstruct the value that we would be able to get
     * and compute the difference with the actual value.
     * 
     * @param ctx only the trgdata information are used, the source is considered empty
     * @param details_max the maximum of the local detail coefficients
     */
    void Detail_(m_ptr<const interp_ctx_t> ctx, m_ptr<real_t> details_max) const override {
        m_assert(*details_max == 0.0, "the value must be 0.0");
        //-------------------------------------------------------------------------
        // the size is know @ compiler time
        constexpr sid_t gs_lim = (len_gs_<TN, TNT> / 2 - 1);

        const real_t one      = 1.0;
        const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
        const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

        const real_t*       sdata = ctx->sdata.Read();
        const real_t* const gs    = gs_<TN, TNT> + gs_lim;

        auto op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const lda_t iy = m_sign(i1) * (i1 % 2);
            const lda_t ix = m_sign(i0) * (i0 % 2);
            const lda_t iz = m_sign(i2) * (i2 % 2);
            m_assert(ix == 0 || ix == 1, "this are the two possible values");
            m_assert(iy == 0 || iy == 1, "this are the two possible values");
            m_assert(iz == 0 || iz == 1, "this are the two possible values");

            const lid_t   i0_s   = (i0 - ix);
            const lid_t   i1_s   = (i1 - iy);
            const lid_t   i2_s   = (i2 - iz);
            const real_t* lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);

            // get the filter, depending on if I am odd or even
            const real_t* const gs_x         = (ix == 1) ? (gs) : (&one);
            const real_t* const gs_y         = (iy == 1) ? (gs) : (&one);
            const real_t* const gs_z         = (iz == 1) ? (gs) : (&one);
            const bidx_t        lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
            const bidx_t        lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

            // if one dim is even, id = 0, -> gs[0] = 1 and that's it
            // if one dim is odd, id = 1, -> we loop on gs, business as usual
            real_t interp = 0.0;
            for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                    for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                        const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                        interp += fact * lsdata[m_idx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->srcstr)];
                    }
                }
            }
            real_t detail = sdata[m_idx(i0, i1, i2, 0, ctx->srcstr)] - interp;

            // check the maximum
            (*details_max) = m_max(fabs(detail), (*details_max));
        };

        // run that
        *details_max = 0.0;
        for_loop(&op, start, end);
    };
    // /**
    //  * @brief gets the norm-2 of the detail coefficients over the block
    //  * 
    //  * we use the gs filter, which is nothing but the dual lifting step. We reconstruct the value that we would be able to get
    //  * and compute the difference with the actual value.
    //  * 
    //  * @param ctx only the trgdata information are used, the source is considered empty
    //  * @param details_max the maximum of the local detail coefficients
    //  */
    // void Detail_2_(m_ptr<const interp_ctx_t> ctx, m_ptr<real_t> details_max) const override {
    //     m_assert(*details_max == 0.0, "the value must be 0.0");
    //     //-------------------------------------------------------------------------
    //     // the size is know @ compiler time
    //     constexpr sid_t gs_lim = (len_gs_<TN, TNT> / 2 - 1);

    //     const real_t one      = 1.0;
    //     const real_t vol      = ctx->alpha;
    //     const lid_t  start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    //     const lid_t  end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

    //     const real_t*       sdata = ctx->sdata.Read();
    //     const real_t* const gs    = gs_<TN, TNT> + gs_lim;

    //     auto op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
    //         // get 0 if odd, 1 if even (even if negative!!)
    //         const lda_t iy = m_sign(i1) * (i1 % 2);
    //         const lda_t ix = m_sign(i0) * (i0 % 2);
    //         const lda_t iz = m_sign(i2) * (i2 % 2);
    //         m_assert(ix == 0 || ix == 1, "this are the two possible values");
    //         m_assert(iy == 0 || iy == 1, "this are the two possible values");
    //         m_assert(iz == 0 || iz == 1, "this are the two possible values");

    //         const lid_t   i0_s   = (i0 - ix);
    //         const lid_t   i1_s   = (i1 - iy);
    //         const lid_t   i2_s   = (i2 - iz);
    //         const real_t* lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);

    //         // get the filter, depending on if I am odd or even
    //         const real_t* const gs_x         = (ix == 1) ? (gs) : (&one);
    //         const real_t* const gs_y         = (iy == 1) ? (gs) : (&one);
    //         const real_t* const gs_z         = (iz == 1) ? (gs) : (&one);
    //         const bidx_t        lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
    //         const bidx_t        lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

    //         // if one dim is even, id = 0, -> gs[0] = 1 and that's it
    //         // if one dim is odd, id = 1, -> we loop on gs, business as usual
    //         real_t interp = 0.0;
    //         for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
    //             for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
    //                 for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
    //                     const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
    //                     interp += fact * lsdata[m_idx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->srcstr)];
    //                 }
    //             }
    //         }
    //         real_t detail = sdata[m_idx(i0, i1, i2, 0, ctx->srcstr)] - interp;

    //         // check the maximum
    //         *details_max += detail * detail * vol;
    //     };

    //     // run that
    //     *details_max = 0.0;
    //     for_loop(&op, start, end);

    //     // take the sqrt of the result
    //     *details_max = sqrt(*details_max);
    //     //-------------------------------------------------------------------------
    // };

    /**
     * @brief Compute the details using the source field and write them in the target field
     * 
     * @param ctx the interpolation context
     */
    void WriteDetail_(m_ptr<const interp_ctx_t> ctx) const override {
        m_assert(ctx->srcstart[0] == 0 && ctx->srcend[0] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->srcstart[0], ctx->srcend[0]);
        m_assert(ctx->srcstart[1] == 0 && ctx->srcend[1] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->srcstart[1], ctx->srcend[1]);
        m_assert(ctx->srcstart[2] == 0 && ctx->srcend[2] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->srcstart[2], ctx->srcend[2]);
        m_assert(ctx->trgstart[0] == 0 && ctx->trgend[0] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->trgstart[0], ctx->trgend[0]);
        m_assert(ctx->trgstart[1] == 0 && ctx->trgend[1] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->trgstart[1], ctx->trgend[1]);
        m_assert(ctx->trgstart[2] == 0 && ctx->trgend[2] == M_N, "the start index = %d and the end one = %d must be 0 and M_N", ctx->trgstart[2], ctx->trgend[2]);
        //-------------------------------------------------------------------------
        // the size is know @ compiler time
        constexpr sid_t gs_lim = (len_gs_<TN, TNT> / 2 - 1);

        const real_t        one = 1.0;
        const real_t* const gs  = gs_<TN, TNT> + gs_lim;

        // get the target + criterion field
        const real_t* sdata = ctx->sdata.Read();
        real_t*       tdata = ctx->tdata.Write();

        // declare the lambda to run
        auto lambda = [=, &tdata](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get 0 if odd, 1 if even (even if negative!!)
            const lda_t iy = m_sign(i1) * (i1 % 2);
            const lda_t ix = m_sign(i0) * (i0 % 2);
            const lda_t iz = m_sign(i2) * (i2 % 2);
            m_assert(ix == 0 || ix == 1, "this are the two possible values");
            m_assert(iy == 0 || iy == 1, "this are the two possible values");
            m_assert(iz == 0 || iz == 1, "this are the two possible values");

            const lid_t   i0_s   = (i0 - ix);
            const lid_t   i1_s   = (i1 - iy);
            const lid_t   i2_s   = (i2 - iz);
            const real_t* lsdata = sdata + m_idx(i0_s, i1_s, i2_s, 0, ctx->srcstr);

            // get the filter, depending on if I am odd or even
            const real_t* const gs_x         = (ix == 1) ? (gs) : (&one);
            const real_t* const gs_y         = (iy == 1) ? (gs) : (&one);
            const real_t* const gs_z         = (iz == 1) ? (gs) : (&one);
            const bidx_t        lim_start[3] = {(gs_lim)*ix, (gs_lim)*iy, (gs_lim)*iz};
            const bidx_t        lim_end[3]   = {(gs_lim + 1) * ix, (gs_lim + 1) * iy, (gs_lim + 1) * iz};

            // if one dim is even, id = 0, -> gs[0] = 1 and that's it
            // if one dim is odd, id = 1, -> we loop on gs, business as usual
            real_t interp = 0.0;
            for (bidx_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                for (bidx_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                    for (bidx_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                        const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                        interp += fact * lsdata[m_idx(id0 * 2, id1 * 2, id2 * 2, 0, ctx->srcstr)];
                    }
                }
            }
            real_t* ltdata = tdata + m_idx(i0, i1, i2, 0, ctx->trgstr);
            ltdata[0]      = sdata[m_idx(i0, i1, i2, 0, ctx->srcstr)] - interp;

            // check that we retrieve the original value if we are a scaling coef
            m_assert(!(ix == 0 && iy == 0 && iz == 0 && tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)] != 0.0), "the target value should be 0.0 instead of %e", tdata[m_idx(i0, i1, i2, 0, ctx->trgstr)]);
        };

        // do the loop
        for_loop<0, M_N>(&lambda);
        //-------------------------------------------------------------------------
    };
};

#endif  // SRC_INTERPOLATING_WAVELET_HPP_