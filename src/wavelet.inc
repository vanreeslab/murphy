#ifndef SRC_WAVELET_IPP_
#define SRC_WAVELET_IPP_

#include <cmath>

#include "murphy.hpp"
#include "wavelet.hpp"

using std::pow;

/**
 * @brief return the biggest detail coefficient as a refinement/coarsening criterion
 * 
 * @tparam order 
 * @param block 
 * @param data 
 * @return real_t 
 */
template <int N, int Nt>
real_t Wavelet<N,Nt>::Criterion(MemLayout* block, real_p data) {
    //-------------------------------------------------------------------------
    real_t details_max[8] = {0};
    // get memory details
    Details(block, data, details_max);

    // get the max out of all the details
    real_t criterion = 0.0;
    for (int id = 0; id < 7; id++) {
        criterion = m_max(criterion, details_max[id]);
    }
    return criterion;
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
void Wavelet<N,Nt>::Details(MemLayout* block, real_p data, real_t* details_max) {
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
    m_assert((block->gs() + block->start(0)) >= (2 * NGhostCoarse()), "the detail computation requires at least %d ghost points", 2 * NGhostCoarse());
    m_assert((block->gs() + block->start(1)) >= (2 * NGhostCoarse()), "the detail computation requires at least %d ghost points", 2 * NGhostCoarse());
    m_assert((block->gs() + block->start(2)) >= (2 * NGhostCoarse()), "the detail computation requires at least %d ghost points", 2 * NGhostCoarse());
    Detail_(&ctx, details_max);
    //-------------------------------------------------------------------------
}

/**
 * @brief copy the value of the source memory to the target memory
 * 
 * @tparam order 
 * @param ctx 
 */
template <int N,int Nt>
void Wavelet<N,Nt>::Copy_(const interp_ctx_t* ctx) const {
    //-------------------------------------------------------------------------
    // ensure alignment for target, constant and source
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    m_assume_aligned(ctx->cdata);

    const real_t alpha = ctx->alpha;
    // do the copy
    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                m_assert((ik0 >= ctx->srcstart[0]) && (ik0 < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d<%d", ik0, ctx->srcstart[0], ik0, ctx->srcend[0]);
                m_assert((ik1 >= ctx->srcstart[1]) && (ik1 < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d<%d", ik1, ctx->srcstart[1], ik1, ctx->srcend[1]);
                m_assert((ik2 >= ctx->srcstart[2]) && (ik2 < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d<%d", ik2, ctx->srcstart[2], ik2, ctx->srcend[2]);
                // get the current parent's data
                real_p       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);
                // do the simple copy
                ltdata[0] = alpha * lcdata[0] + lsdata[0];
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief coarsen the values of the source memory to gather them in the target memory.
 * It can refine on up to 2 successive levels
 * 
 * @tparam order 
 * @param ctx 
 * @param dlvl the number of levels we have to refine
 */
template <int N, int Nt>
void Wavelet<N, Nt>::Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    m_assume_aligned(ctx->cdata);

    const real_t  alpha  = ctx->alpha;
    const sid_t   ha_lim = (dlvl == 1) ? (len_ha_ / 2) : (dlvl == 2) ? (len_ha_2_ / 2) : 0;
    const real_t* ha     = (dlvl == 1) ? ha_ : (dlvl == 2) ? ha_2_ : nullptr;

    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                // do some checks
                m_assert(((ik0 - ha_lim) >= ctx->srcstart[0]) && ((ik0 + ha_lim) <= ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - ha_lim, ctx->srcstart[0], ik0 + ha_lim, ctx->srcend[0]);
                m_assert(((ik1 - ha_lim) >= ctx->srcstart[1]) && ((ik1 + ha_lim) <= ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - ha_lim, ctx->srcstart[1], ik1 + ha_lim, ctx->srcend[1]);
                m_assert(((ik2 - ha_lim) >= ctx->srcstart[2]) && ((ik2 + ha_lim) <= ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - ha_lim, ctx->srcstart[2], ik2 + ha_lim, ctx->srcend[2]);
                //get the local adress of the source, the target and the constant
                real_p       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lsdata = ctx->sdata + m_sidx(2 * dlvl * ik0, 2 * dlvl * ik1, 2 * dlvl * ik2, 0, ctx->srcstr);

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
 * @tparam order 
 * @param ctx 
 */
template <int N, int Nt>
void Wavelet<N,Nt>::Refine_(const interp_ctx_t* ctx) const {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    m_assume_aligned(ctx->cdata);

    const real_t  alpha  = ctx->alpha;
    const sid_t   gs_lim = (len_gs_ / 2);
    const real_t* gs     = gs_;

    m_verb("gslim = %d",gs_lim);

    for (lid_t ik2 = ctx->trgstart[2] / 2; ik2 < ctx->trgend[2] / 2; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                m_assert((ik0 - gs_lim + 1 >= ctx->srcstart[0]) && (ik0 + gs_lim < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - gs_lim, ctx->srcstart[0], ik0 + gs_lim, ctx->srcend[0]);
                m_assert((ik1 - gs_lim + 1 >= ctx->srcstart[1]) && (ik1 + gs_lim < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - gs_lim, ctx->srcstart[1], ik1 + gs_lim, ctx->srcend[1]);
                m_assert((ik2 - gs_lim + 1 >= ctx->srcstart[2]) && (ik2 + gs_lim < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - gs_lim, ctx->srcstart[2], ik2 + gs_lim, ctx->srcend[2]);
                //get the local adress of the arrays
                real_p       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const real_p lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                const real_p lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

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
                for (sid_t id0 = -gs_lim + 1; id0 <= gs_lim; id0++) {
                    ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(id0, 0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0,id0, 0, 0, ctx->srcstr)];
                    ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += gs[id0] * lsdata[m_sidx(0, 0,id0, 0, ctx->srcstr)];
                }
                // point (1,1,0) and (0,1,1) and (1,0,1)
                for (sid_t id1 = -gs_lim + 1; id1 <= gs_lim; id1++) {
                    for (sid_t id0 = -gs_lim + 1; id0 <= gs_lim; id0++) {
                        ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, id1, 0, 0, ctx->srcstr)];
                        ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(0, id0, id1, 0, ctx->srcstr)];
                        ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * lsdata[m_sidx(id0, 0, id1, 0, ctx->srcstr)];
                    }
                }
                // point (1,1,1)
                for (sid_t id2 = -gs_lim + 1; id2 <= gs_lim; id2++) {
                    for (sid_t id1 = -gs_lim + 1; id1 <= gs_lim; id1++) {
                        for (sid_t id0 = -gs_lim + 1; id0 <= gs_lim; id0++) {
                            ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += gs[id0] * gs[id1] * gs[id2] * lsdata[m_sidx(id0, id1, id2, 0, ctx->srcstr)];
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
 * @param details_inf_norm the maximum of the local detail coefficients: (d_x,d_y,d_z,d_xy,d_yz,d_xz,d_xyz,mean)
 */
template <int N, int Nt>
void Wavelet<N,Nt>::Detail_(const interp_ctx_t* ctx, real_t* details_inf_norm) const {
    //-------------------------------------------------------------------------
    m_assume_aligned(ctx->tdata);

    const real_t  alpha  = ctx->alpha;
    const sid_t   ga_lim = len_ga_/2;
    const real_t* ga     = ga_;

    // for each of the data for the considered children
    for (int ik2 = ctx->trgstart[2] / 2; ik2 < ctx->trgend[2] / 2; ik2++) {
        for (int ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
            for (int ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                // get the current data
                real_p ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                m_assume_aligned(ltdata);

                // init the details
                real_t detail[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

                // detail(0,0,0) = scaling!
                detail[0] = ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];

                // point (1,0,0) and (0,1,0) and (0,0,1)
                for (sid_t id0 = -ga_lim; id0 <= ga_lim; id0++) {
                    detail[1] += ga[id0] * ltdata[m_sidx(1 + id0, 0, 0, 0, ctx->trgstr)];
                    detail[2] += ga[id0] * ltdata[m_sidx(0, 1 + id0, 0, 0, ctx->trgstr)];
                    detail[3] += ga[id0] * ltdata[m_sidx(0, 0, 1 + id0, 0, ctx->trgstr)];
                }
                // point (1,1,0) and (0,1,1) and (1,0,1)
                for (sid_t id1 = -ga_lim; id1 <= ga_lim; id1++) {
                    for (sid_t id0 = -ga_lim; id0 <= ga_lim; id0++) {
                        detail[4] += ga[id0] * ga[id1] * ltdata[m_sidx(1 + id0, 1 + id1, 0, 0, ctx->trgstr)];
                        detail[5] += ga[id0] * ga[id1] * ltdata[m_sidx(0, 1 + id0, 1 + id1, 0, ctx->trgstr)];
                        detail[6] += ga[id0] * ga[id1] * ltdata[m_sidx(1 + id0, 0, 1 + id1, 0, ctx->trgstr)];
                    }
                }
                // point (1,1,1)
                for (sid_t id2 = -ga_lim; id2 <= ga_lim; id2++) {
                    for (sid_t id1 = -ga_lim; id1 <= ga_lim; id1++) {
                        for (sid_t id0 = -ga_lim; id0 <= ga_lim; id0++) {
                            detail[7] += ga[id0] * ga[id1] * ga[id2] * ltdata[m_sidx(1 + id0, 1 + id1, 1 + id2, 0, ctx->trgstr)];
                        }
                    }
                }
                // store the max, relatively to the current value
                for (int id = 1; id < 8; id++) {
                    details_inf_norm[id] = m_max(std::fabs(detail[id]), details_inf_norm[id]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

#endif