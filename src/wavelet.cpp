#ifndef SRC_WAVELET_IPP_
#define SRC_WAVELET_IPP_

#include "wavelet.hpp"

#include <cmath>

#include "murphy.hpp"

// this is needed, see https://stackoverflow.com/questions/8016780/undefined-reference-to-static-constexpr-char
constexpr lid_t  Wavelet::len_ha_;
constexpr lid_t  Wavelet::len_ga_;
constexpr lid_t  Wavelet::len_gs_;
constexpr real_t Wavelet::ha_[];
constexpr real_t Wavelet::ga_[];
constexpr real_t Wavelet::gs_[];

/**
 * @brief return the biggest detail coefficient (infinite norm) as a refinement/coarsening criterion
 * 
 * The max/min is computed on an extension of the memory Layout (by lifting_len/2).
 * The reason for that comes from the coarsening of a block. 
 * While coarsening a block, we implicitly assume that all detail coefficients involved in the refinement are 0.0.
 * By doing so, we can coarsen the block, forget the detail coefficients and still retreive a perfect information.
 * Practically, it means that the lifting step, i.e. the contribution of the detail coefficients to my scaling coeff, is useless.
 * 
 * The lifting step is driven by the lifting coefficient, whose total length is (2*Nt-1).
 * In front of the block, I need to ensure that the details from my neighbor are 0.0, which means (2*Nt-1)/2 detail are zero
 * At the back of the block, I need to ensure (2*Nt-1)/2 -1 details are 0 (as the last point is a detail)
 *
 * @param block the block to analyze
 * @param data the data
 * @return real_t the infinite norm of the max detail coefficient in the extended region
 */
real_t Wavelet::Criterion(MemLayout* block, data_ptr data) {
    //-------------------------------------------------------------------------
    // get the extended memory layout
    constexpr lid_t lift_len = (2 * M_WAVELET_NT - 1);
    lid_t           start[3];
    lid_t           end[3];
    for (lda_t id = 0; id < 3; id++) {
        start[id] = block->start(id) - m_max(lift_len / 2, 0);
        end[id]   = block->end(id) + m_max((lift_len / 2) - 1, 0);
    }
    SubBlock extended_block(block->gs(), block->stride(), start, end);

    // get the detail coefficients
    real_t details_max = 0.0;
    Details(&extended_block, data, &details_max);

    return details_max;
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the max detail coefficients of a given MemLayout
 * 
 * @param block the block on which we computed
 * @param data the memory pointer to the point (0,0,0) of that block
 * @param details_max an array of size 8 that will contain the detail coefficients: dx, dy, dz, dxy, dyz, dxz, dxyz, mean
 */
void Wavelet::Details(MemLayout* block, data_ptr data, real_t* details_max) {
    //-------------------------------------------------------------------------
    interp_ctx_t ctx;
    // get memory details
    for (int id = 0; id < 3; id++) {
#ifndef NDEBUG
        ctx.srcstart[id] = -1;
        ctx.srcend[id]   = -2;
#endif
        ctx.trgstart[id] = block->start(id);
        ctx.trgend[id]   = block->end(id);
    }
    ctx.srcstr = -1;
    ctx.sdata  = nullptr;
    ctx.trgstr = block->stride();
    ctx.tdata  = data;
    Detail_(&ctx, details_max);
    //-------------------------------------------------------------------------
}

/**
 * @brief coarsen the values of the source memory to gather them in the target memory.
 * 
 * @tparam N the number of vanishing moment
 * @tparam Nt the order of interpolation
 * @param ctx the interpolation context
 */
void Wavelet::Coarsen_(const interp_ctx_t* ctx) {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    // m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    // m_assume_aligned(ctx->cdata);

    const real_t  alpha  = ctx->alpha;
    const lid_t   ha_lim = len_ha_ / 2;
    const_mem_ptr ha     = Wavelet::ha_ + ha_lim;

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
                            ltdata[0] += lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)] * ha[id0] * ha[id1] * ha[id2];
                        }
                    }
                }
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
void Wavelet::Refine_(const interp_ctx_t* ctx) {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    const real_t  alpha    = ctx->alpha;
    const sid_t   gs_lim   = (len_gs_ / 2) - 1;
    const_mem_ptr gs       = Wavelet::gs_ + gs_lim;
    const lid_t   start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    const lid_t   end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};
    const real_t  one      = 1.0;

    // for each of the data for the needed target
    for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
        for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
            for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
                // get the current data
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx((ik0 / 2), (ik1 / 2), (ik2 / 2), 0, ctx->srcstr);

                // get the needed filters and the loop lim
                const lda_t iy = std::fabs(ik1 % 2);
                const lda_t ix = std::fabs(ik0 % 2);
                const lda_t iz = std::fabs(ik2 % 2);
                m_assert(ix == 0 || ix == 1, "this are the two possible values");
                m_assert(iy == 0 || iy == 1, "this are the two possible values");
                m_assert(iz == 0 || iz == 1, "this are the two possible values");

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
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];

                // if one dim is even, id = 0, -> gs[0] = 1 and that's it
                // if one dim is odd, id = 1, -> we loop on gs, business as usual
                for (sid_t id2 = -lim_start[2]; id2 <= lim_end[2]; ++id2) {
                    for (sid_t id1 = -lim_start[1]; id1 <= lim_end[1]; ++id1) {
                        for (sid_t id0 = -lim_start[0]; id0 <= lim_end[0]; ++id0) {
                            const real_t fact = gs_x[id0] * gs_y[id1] * gs_z[id2];
                            ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += fact * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                        }
                    }
                }
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
void Wavelet::Detail_(const interp_ctx_t* ctx, real_t* details_max) {
    //-------------------------------------------------------------------------
    const sid_t   ga_lim   = len_ga_ / 2;
    const_mem_ptr ga       = Wavelet::ga_ + ga_lim;
    const real_t  zero     = 0.0;
    const lid_t   start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    const lid_t   end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

    // do a small sanity check
    m_assert(ga[0] == 1.0, "we need that to ensure a correct result bellow");

    // for each of the data for the considered children
    (*details_max) = 0.0;
    for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
        for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
            for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
                // get the current data
                data_ptr ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                // m_assume_aligned(ltdata);

                const lda_t iy = std::fabs(ik1 % 2);
                const lda_t ix = std::fabs(ik0 % 2);
                const lda_t iz = std::fabs(ik2 % 2);
                m_assert(ix == 0 || ix == 1, "this are the two possible values");
                m_assert(iy == 0 || iy == 1, "this are the two possible values");
                m_assert(iz == 0 || iz == 1, "this are the two possible values");

                // get the filter, depending on if I am odd or even
                const sid_t lim[3] = {ga_lim * ix, ga_lim * iy, ga_lim * iz};

                // if one dim is even, id = 0, -> ga = 1 and that's it
                // if one dim is odd, id=1, -> we loop on ga, business as usual
                real_t detail = 0.0;
                for (sid_t id2 = -lim[2]; id2 <= lim[2]; ++id2) {
                    for (sid_t id1 = -lim[1]; id1 <= lim[1]; ++id1) {
                        for (sid_t id0 = -lim[0]; id0 <= lim[0]; ++id0) {
                            detail += ga[id0] * ga[id1] * ga[id2] * ltdata[m_sidx(id0, id1, id2, 0, ctx->trgstr)];
                        }
                    }
                }
                // update the max detail if needed and set the detail to 0.0 if we are even everywhere
                detail         = ((ix + iy + iz) > 0) ? (detail) : (0.0);
                (*details_max) = m_max(std::fabs(detail), (*details_max));
            }
        }
    }
    //-------------------------------------------------------------------------
}

#endif