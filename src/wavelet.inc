#ifndef SRC_WAVELET_IPP_
#define SRC_WAVELET_IPP_

#include <cmath>

#include "murphy.hpp"
#include "wavelet.hpp"

using std::pow;

/**
 * @brief return the biggest detail coefficient as a refinement/coarsening criterion
 * 
 * The max/min is computed on an extension of the Memory Layout (by gs_len/2).
 * The reason for that comes from the refinement of a block.
 * To refine a block, we have to assume that all the detail coefficients involved in the refinement are 0.0.
 * Taking this way of thinking in the opposite direction, we may discard a block, only if the detail coefficients involved in its refinement are 0.0.
 * Practically, it means that the lifting step, i.e. the contribution of the detail coefficients to my scaling coeff are 0.0.
 * 
 * The lifting step is driven by the lifting coefficient, whose total length is (2*Nt-1).
 * in front, I need to ensure that the details from my neighbor are 0.0, which means (2*Nt-1)/2 detail are zero
 * at the back, I need to ensure (2*Nt-1)/2 -1 details are 0 (as the last point is a detail)
 * 
 * 
 * @tparam order 
 * @param block 
 * @param data 
 * @return real_t 
 */
template <int N, int Nt>
real_t Wavelet<N, Nt>::Criterion(MemLayout* block, data_ptr data) {
    //-------------------------------------------------------------------------
    // get the extended memory layout
    const lid_t lift_len = (2 * Nt - 1);
    lid_t       start[3];
    lid_t       end[3];
    for (lda_t id = 0; id < 3; id++) {
        start[id] = block->start(id) - (lift_len / 2);
        end[id]   = block->end(id) + (lift_len / 2) - 1;
    }
    SubBlock extended_block(block->gs(), block->stride(), start, end);

    // get the detail coefficients
    real_t details_max = 0.0;
    Details(&extended_block, data, &details_max);

    // get the max out of all the details
    // real_t criterion = 0.0;
    // for (lda_t id = 0; id < 7; id++) {
    //     criterion = m_max(criterion, details_max[id]);
    // }
    // m_log("hello the details: %f %f %f %f %f %f %f", details_max[0], details_max[1], details_max[2], details_max[3], details_max[4], details_max[5], details_max[6]);
    return details_max;
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the detail coefficients of a given MemLayout
 * 
 * @tparam order the order of the wavelet
 * @param block the block on which we computed
 * @param data the memory pointer to the point (0,0,0) of that block
 * @param details_max an array of size 8 that will contain the detail coefficients: dx, dy, dz, dxy, dyz, dxz, dxyz, mean
 */
template <int N, int Nt>
void Wavelet<N, Nt>::Details(MemLayout* block, data_ptr data, real_t* details_max) {
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
    // get the details
    // m_assert((block->gs() + block->start(0)) >= (2 * NGhostCoarse()), "the detail computation requires at least %d ghost points", 2 * NGhostCoarse());
    // m_assert((block->gs() + block->start(1)) >= (2 * NGhostCoarse()), "the detail computation requires at least %d ghost points", 2 * NGhostCoarse());
    // m_assert((block->gs() + block->start(2)) >= (2 * NGhostCoarse()), "the detail computation requires at least %d ghost points", 2 * NGhostCoarse());
    Detail_(&ctx, details_max);
    //-------------------------------------------------------------------------
}

/**
 * @brief coarsen the values of the source memory to gather them in the target memory.
 * It can refine on up to 2 successive levels
 * 
 * @tparam N the number of vanishing moment
 * @tparam Nt the order of interpolation
 * @param ctx the interpolation context
 */
template <int N, int Nt>
void Wavelet<N, Nt>::Coarsen_(const interp_ctx_t* ctx) {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    // m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    // m_assume_aligned(ctx->cdata);

    const real_t  alpha  = ctx->alpha;
    const sid_t   ha_lim = len_ha_ / 2;
    const mem_ptr ha     = ha_ + ha_lim;

    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
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
 * @tparam N the number of vanishing moment
 * @tparam Nt the order of interpolation
 * @param ctx the interpolation context
 */
template <int N, int Nt>
void Wavelet<N, Nt>::Refine_(const interp_ctx_t* ctx) {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    // m_assume_aligned(ctx->tdata);
    // m_assume_aligned(ctx->sdata);
    // m_assume_aligned(ctx->cdata);

    const real_t  alpha  = ctx->alpha;
    const sid_t   gs_lim = (len_gs_ / 2) - 1;
    const mem_ptr gs     = gs_ + gs_lim;

    // m_log("gslim = %d", gs_lim);

    // m_log("trg = %d %d %d to %d %d %d ", ctx->trgstart[0] / 2, ctx->trgstart[1] / 2, ctx->trgstart[2] / 2, ctx->trgend[0] / 2, ctx->trgend[1] / 2, ctx->trgend[2] / 2);

    // m_assert((ctx->trgstart[0]%2==0) && (ctx->trgend[0]%2==0),"the starting %d and ending %d location must be even in dir %d",ctx->trgstart[0],ctx->trgend[0],0);
    // m_assert((ctx->trgstart[1]%2==0) && (ctx->trgend[1]%2==0),"the starting %d and ending %d location must be even in dir %d",ctx->trgstart[1],ctx->trgend[1],1);
    // m_assert((ctx->trgstart[2]%2==0) && (ctx->trgend[2]%2==0),"the starting %d and ending %d location must be even in dir %d",ctx->trgstart[2],ctx->trgend[2],2);
    m_assert((ctx->trgstart[0] % 2 == 0), "the starting location = %d must be even in dir %d", ctx->trgstart[0], 0);
    m_assert((ctx->trgstart[1] % 2 == 0), "the starting location = %d must be even in dir %d", ctx->trgstart[1], 1);
    m_assert((ctx->trgstart[2] % 2 == 0), "the starting location = %d must be even in dir %d", ctx->trgstart[2], 2);

    for (lid_t ik2 = ctx->trgstart[2] / 2; ik2 < ctx->trgend[2] / 2; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                m_assert(((ik0 - gs_lim) >= ctx->srcstart[0]) && ((ik0 + gs_lim + 1) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 1, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(id0, 0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // point (1,1,0) and (0,1,1) and (1,0,1)
                for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                    for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                        ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, id1, 0, 0, ctx->srcstr)];
                        ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                        ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, 0, id1, 0, ctx->srcstr)];
                    }
                }
                // point (1,1,1)
                for (sid_t id2 = (-gs_lim); id2 <= (gs_lim + 1); id2++) {
                    for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                        for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                            ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * gs[id2] * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                        }
                    }
                }
            }
            // I cannot go in x=1
            for (lid_t ik0 = ctx->trgend[0] / 2; ik0 < (ctx->trgend[0] + 1) / 2; ik0++) {
                m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // point (1,1,0) and (0,1,1) and (1,0,1)
                for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                    for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                        ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                    }
                }
            }
        }
        // I cannot go in y=1
        for (lid_t ik1 = ctx->trgstart[1] / 2; ik1 < (ctx->trgend[1] + 1) / 2; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                m_assert(((ik0 - gs_lim) >= ctx->srcstart[0]) && ((ik0 + gs_lim + 1) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                // m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(id0, 0, 0, 0, ctx->srcstr)];
                    // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // point (1,1,0) and (0,1,1) and (1,0,1)
                for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                    for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                        // ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, id1, 0, 0, ctx->srcstr)];
                        // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                        ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, 0, id1, 0, ctx->srcstr)];
                    }
                }
                // // point (1,1,1)
                // for (sid_t id2 = (-gs_lim); id2 <= (gs_lim + 1); id2++) {
                //     for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //         for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //             ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * gs[id2] * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                //         }
                //     }
                // }
            }
            // I cannot go in x=1
            for (lid_t ik0 = ctx->trgend[0] / 2; ik0 < (ctx->trgend[0] + 1) / 2; ik0++) {
                m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // // point (1,1,0) and (0,1,1) and (1,0,1)
                // for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //     for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //         ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                //     }
                // }
            }
        }
    }
    // I cannot go to z=1
    for (lid_t ik2 = ctx->trgstart[2] / 2; ik2 < (ctx->trgend[2] + 1) / 2; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                m_assert(((ik0 - gs_lim) >= ctx->srcstart[0]) && ((ik0 + gs_lim + 1) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                // m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(id0, 0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // point (1,1,0) and (0,1,1) and (1,0,1)
                for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                    for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                        ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, id1, 0, 0, ctx->srcstr)];
                        // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                        // ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, 0, id1, 0, ctx->srcstr)];
                    }
                }
                // // point (1,1,1)
                // for (sid_t id2 = (-gs_lim); id2 <= (gs_lim + 1); id2++) {
                //     for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //         for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //             ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * gs[id2] * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                //         }
                //     }
                // }
            }
            // I cannot go in x=1
            for (lid_t ik0 = ctx->trgend[0] / 2; ik0 < (ctx->trgend[0] + 1) / 2; ik0++) {
                m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                // m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // // point (1,1,0) and (0,1,1) and (1,0,1)
                // for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //     for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //         ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                //     }
                // }
            }
        }
        // I cannot go in y=1
        for (lid_t ik1 = ctx->trgstart[1] / 2; ik1 < (ctx->trgend[1] + 1) / 2; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                m_assert(((ik0 - gs_lim) >= ctx->srcstart[0]) && ((ik0 + gs_lim + 1) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                // m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                // m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                data_ptr       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const data_ptr lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                    ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(id0, 0, 0, 0, ctx->srcstr)];
                    // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                    // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                }
                // // point (1,1,0) and (0,1,1) and (1,0,1)
                // for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //     for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //         // ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, id1, 0, 0, ctx->srcstr)];
                //         // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                //         ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, 0, id1, 0, ctx->srcstr)];
                //     }
                // }
                // // point (1,1,1)
                // for (sid_t id2 = (-gs_lim); id2 <= (gs_lim + 1); id2++) {
                //     for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //         for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //             ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * gs[id2] * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
                //         }
                //     }
                // }
            }
            // I cannot go in x=1
            for (lid_t ik0 = ctx->trgend[0] / 2; ik0 < (ctx->trgend[0] + 1) / 2; ik0++) {
                m_assert(((ik1 - gs_lim) >= ctx->srcstart[1]) && ((ik1 + gs_lim + 1) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert(((ik2 - gs_lim) >= ctx->srcstart[2]) && ((ik2 + gs_lim + 1) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                real_p       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const real_p lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const real_p lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to the constant
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
                // ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];

                // point (0,0,0)
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];

                // // point (1,0,0) and (0,1,0) and (0,0,1)
                // for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //     // ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, id0, 0, 0, ctx->srcstr)];
                //     ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0, id0, 0, ctx->srcstr)];
                // }
                // // point (1,1,0) and (0,1,1) and (1,0,1)
                // for (sid_t id1 = (-gs_lim); id1 <= (gs_lim + 1); id1++) {
                //     for (sid_t id0 = (-gs_lim); id0 <= (gs_lim + 1); id0++) {
                //         ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                //     }
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
template <int N, int Nt>
void Wavelet<N, Nt>::Detail_(const interp_ctx_t* ctx, real_t* details_max) {
    //-------------------------------------------------------------------------
    m_assume_aligned(ctx->tdata);

    const sid_t   ga_lim = len_ga_ / 2;
    const mem_ptr ga     = ga_ + ga_lim;
    m_assert(ga[0] == 1.0, "this is annoying as the looping procedure is based on that assumption");

    const lid_t start[3] = {ctx->trgstart[0], ctx->trgstart[1], ctx->trgstart[2]};
    const lid_t end[3]   = {ctx->trgend[0], ctx->trgend[1], ctx->trgend[2]};

    (*details_max) = 0.0;

    // for each of the data for the considered children
    for (lid_t ik2 = start[2]; ik2 < end[2]; ++ik2) {
        for (lid_t ik1 = start[1]; ik1 < end[1]; ++ik1) {
            for (lid_t ik0 = start[0]; ik0 < end[0]; ++ik0) {
                // get the current data
                data_ptr ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                m_assume_aligned(ltdata);

                const lda_t ix = ik0 % 2;
                const lda_t iy = ik0 % 2;
                const lda_t iz = ik0 % 2;

                // if one dim is even, id = 0, -> ga = 1 and that's it
                // if one dim is odd, id=1, -> we loop on ga, business as usual
                real_t detail = 0.0;
                for (sid_t id2 = -(ga_lim * iz); id2 <= (ga_lim * iz); ++id2) {
                    for (sid_t id1 = -(ga_lim * iy); id1 <= (ga_lim * iy); ++id1) {
                        for (sid_t id0 = -(ga_lim * ix); id0 <= (ga_lim * ix); ++id0) {
                            detail += ga[id0] * ga[id1] * ga[id2] * ltdata[m_sidx(1 + id0, 1 + id1, 1 + id2, 0, ctx->trgstr)];
                        }
                    }
                }
                // update the max detail if needed
                (*details_max) = m_max(std::fabs(detail), (*details_max));

                // // init the details
                // real_t detail[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                // // detail(0,0,0) = scaling!
                // detail[0] = ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];

                // // point (1,0,0) and (0,1,0) and (0,0,1)
                // for (sid_t id0 = -ga_lim; id0 <= ga_lim; ++id0) {
                //     detail[1] += ga[id0] * ltdata[m_sidx(1 + id0, 0, 0, 0, ctx->trgstr)];
                //     detail[2] += ga[id0] * ltdata[m_sidx(0, 1 + id0, 0, 0, ctx->trgstr)];
                //     detail[3] += ga[id0] * ltdata[m_sidx(0, 0, 1 + id0, 0, ctx->trgstr)];
                // }
                // // point (1,1,0) and (0,1,1) and (1,0,1)
                // for (sid_t id1 = -ga_lim; id1 <= ga_lim; ++id1) {
                //     for (sid_t id0 = -ga_lim; id0 <= ga_lim; ++id0) {
                //         detail[4] += ga[id0] * ga[id1] * ltdata[m_sidx(1 + id0, 1 + id1, 0, 0, ctx->trgstr)];
                //         detail[5] += ga[id0] * ga[id1] * ltdata[m_sidx(0, 1 + id0, 1 + id1, 0, ctx->trgstr)];
                //         detail[6] += ga[id0] * ga[id1] * ltdata[m_sidx(1 + id0, 0, 1 + id1, 0, ctx->trgstr)];
                //     }
                // }
                // // point (1,1,1)
                // for (sid_t id2 = -ga_lim; id2 <= ga_lim; ++id2) {
                //     for (sid_t id1 = -ga_lim; id1 <= ga_lim; ++id1) {
                //         for (sid_t id0 = -ga_lim; id0 <= ga_lim; ++id0) {
                //             detail[7] += ga[id0] * ga[id1] * ga[id2] * ltdata[m_sidx(1 + id0, 1 + id1, 1 + id2, 0, ctx->trgstr)];
                //         }
                //     }
                // }
                // // store the max, relatively to the current value
                // for (int id = 1; id < 8; id++) {
                //     details_inf_norm[id] = m_max(std::fabs(detail[id]), details_inf_norm[id]);
                // }
            }
        }
    }
    //-------------------------------------------------------------------------
}

#endif