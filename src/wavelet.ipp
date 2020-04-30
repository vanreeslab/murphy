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
template <int order>
real_t Wavelet<order>::Criterion(MemLayout* block, real_p data) {
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
template <int order>
void Wavelet<order>::Details(MemLayout* block, real_p data, real_t* details_max) {
    //-------------------------------------------------------------------------
    interp_ctx_t ctx;
    // get memory details
    for (int id = 0; id < 3; id++) {
        ctx.srcstart[id] = block->start(id);
        ctx.srcend[id]   = block->end(id);
        ctx.trgstart[id] = -1;
        ctx.trgend[id]   = -2;
    }
    ctx.srcstr = block->stride();
    ctx.sdata  = data;
    ctx.trgstr = -1;
    ctx.tdata  = nullptr;
    // get the details
    m_assert((block->gs()+block->start(0)) >= (2*NGhostCoarse()), "the detail computation requires at least %d ghost points",2*NGhostCoarse());
    m_assert((block->gs()+block->start(1)) >= (2*NGhostCoarse()), "the detail computation requires at least %d ghost points",2*NGhostCoarse());
    m_assert((block->gs()+block->start(2)) >= (2*NGhostCoarse()), "the detail computation requires at least %d ghost points",2*NGhostCoarse());
    Detail_(&ctx, details_max);
    //-------------------------------------------------------------------------
}

/**
 * @brief copy the value of the source memory to the target memory
 * 
 * @tparam order 
 * @param ctx 
 */
template <int order>
void Wavelet<order>::Copy_(const interp_ctx_t* ctx) const {
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);

    // do the copy
    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                // for every block of 8 child data, get how much information I can get from my parent
                m_assert((ik0 >= ctx->srcstart[0]) && (ik0 < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d<%d", ik0, ctx->srcstart[0], ik0, ctx->srcend[0]);
                m_assert((ik1 >= ctx->srcstart[1]) && (ik1 < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d<%d", ik1, ctx->srcstart[1], ik1, ctx->srcend[1]);
                m_assert((ik2 >= ctx->srcstart[2]) && (ik2 < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d<%d", ik2, ctx->srcstart[2], ik2, ctx->srcend[2]);
                // get the current parent's data
                ctx->tdata[m_sidx(ik0, ik1, ik2, 0, ctx->trgstr)] = ctx->sdata[m_sidx(ik0, ik1, ik2, 0, ctx->srcstr)];
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
template <int order>
void Wavelet<order>::Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const {
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    // get the factor = 0.125^(#d level)
    const real_t fact = pow(0.125, dlvl);
    // for each of the data for the considered children
    for (int ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (int ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (int ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                m_assert(((2 * dlvl) * ik0 >= ctx->srcstart[0]) && ((2 * dlvl) * ik0 + (2 * dlvl - 1) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", (2 * dlvl) * ik0, ctx->srcstart[0], (2 * dlvl) * ik0 + (2 * dlvl - 1), ctx->srcend[0]);
                m_assert(((2 * dlvl) * ik1 >= ctx->srcstart[1]) && ((2 * dlvl) * ik1 + (2 * dlvl - 1) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", (2 * dlvl) * ik1, ctx->srcstart[1], (2 * dlvl) * ik1 + (2 * dlvl - 1), ctx->srcend[1]);
                m_assert(((2 * dlvl) * ik2 >= ctx->srcstart[2]) && ((2 * dlvl) * ik2 + (2 * dlvl - 1) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", (2 * dlvl) * ik2, ctx->srcstart[2], (2 * dlvl) * ik2 + (2 * dlvl - 1), ctx->srcend[2]);
                //get the local adress of the dady
                real_p lsdata = ctx->sdata + m_sidx((2 * dlvl) * ik0, (2 * dlvl) * ik1, (2 * dlvl) * ik2, 0, ctx->srcstr);
                m_assume_aligned(lsdata);
                const real_p ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);

                // do the mean on the current patch
                ltdata[0] = 0.0;
                for (int it2 = 0; it2 < 2 * dlvl; it2++) {
                    for (int it1 = 0; it1 < 2 * dlvl; it1++) {
                        for (int it0 = 0; it0 < 2 * dlvl; it0++) {
                            ltdata[0] += fact * lsdata[m_sidx(it0, it1, it2, 0, ctx->srcstr)];
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
template <int order>
void Wavelet<order>::Refine_(const interp_ctx_t* ctx) const {
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);

    const lid_t   hslen = order / 2;
    const real_t* hs    = hs_ + hslen;
    const real_t* sign  = sgn_hs_ + hslen;

    // for each of the data for the considered children
    for (int ik2 = ctx->trgstart[2] / 2; ik2 < ctx->trgend[2] / 2; ik2++) {
        for (int ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
            for (int ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
                m_assert((ik0 - hslen >= ctx->srcstart[0]) && (ik0 + hslen < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - hslen, ctx->srcstart[0], ik0 + hslen, ctx->srcend[0]);
                m_assert((ik1 - hslen >= ctx->srcstart[1]) && (ik1 + hslen < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - hslen, ctx->srcstart[1], ik1 + hslen, ctx->srcend[1]);
                m_assert((ik2 - hslen >= ctx->srcstart[2]) && (ik2 + hslen < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - hslen, ctx->srcstart[2], ik2 + hslen, ctx->srcend[2]);
                //get the local adress of the dady
                real_p ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
                m_assume_aligned(ltdata);
                const real_p lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

                // reset to 0.0
                ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] = 0.0;
                ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] = 0.0;

                for (int dm2 = -(order / 2); dm2 <= (order / 2); dm2++) {
                    for (int dm1 = -(order / 2); dm1 <= (order / 2); dm1++) {
                        for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
                            // get the current parent's data
                            const real_t ldata = lsdata[m_sidx(dm0, dm1, dm2, 0, ctx->srcstr)];
                            const real_t fact  = hs[dm2] * hs[dm1] * hs[dm0];
                            // we give the information to every wavelet which is inside my "block"
                            ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += fact * ldata;
                            ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += fact * ldata * sign[dm0];
                            ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += fact * ldata * sign[dm1];
                            ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm2];
                            ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm1] * sign[dm2];
                            ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += fact * ldata * sign[dm0] * sign[dm1];
                            ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm0] * sign[dm2];
                            ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm0] * sign[dm1] * sign[dm2];
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
 * @param ctx 
 * @param details_inf_norm the maximum of the local detail coefficients: (d_x,d_y,d_z,d_xy,d_yz,d_xz,d_xyz,mean)
 */
template <int order>
void Wavelet<order>::Detail_(const interp_ctx_t* ctx, real_t* details_inf_norm) const {
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(ctx->sdata);

    const lid_t   hslen = order / 2;
    const real_t* ga    = ga_ + hslen;
    const real_t* sign  = sgn_ga_ + hslen;

    real_t ha_tmp[order] = {0};
    ha_tmp[order / 2]    = 0.5;
    const real_t* ha     = ha_tmp + hslen;

    // for each of the data for the considered children
    for (int ik2 = ctx->srcstart[2] / 2; ik2 < ctx->srcend[2] / 2; ik2++) {
        for (int ik1 = ctx->srcstart[1] / 2; ik1 < ctx->srcend[1] / 2; ik1++) {
            for (int ik0 = ctx->srcstart[0] / 2; ik0 < ctx->srcend[0] / 2; ik0++) {
                //get the local adress of the dady
                real_p lsdata = ctx->sdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);
                m_assume_aligned(lsdata);

                // set the detail coefficient to the
                real_t detail[8];

                // compute all the single coefficients
                detail[0] = 0.0;
                detail[1] = 0.0;
                detail[2] = 0.0;
                for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
                    const real_t fact = ha[0] * ha[0] * ga[dm0];
                    // this is dx
                    detail[0] += lsdata[m_sidx(2 * dm0 + 0, 0, 0, 0, ctx->srcstr)] * fact;
                    detail[0] += lsdata[m_sidx(2 * dm0 + 1, 0, 0, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[0] += lsdata[m_sidx(2 * dm0 + 0, 1, 0, 0, ctx->srcstr)] * fact;
                    detail[0] += lsdata[m_sidx(2 * dm0 + 0, 0, 1, 0, ctx->srcstr)] * fact;
                    detail[0] += lsdata[m_sidx(2 * dm0 + 0, 1, 1, 0, ctx->srcstr)] * fact;
                    detail[0] += lsdata[m_sidx(2 * dm0 + 1, 1, 0, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[0] += lsdata[m_sidx(2 * dm0 + 1, 0, 1, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[0] += lsdata[m_sidx(2 * dm0 + 1, 1, 1, 0, ctx->srcstr)] * fact * sign[dm0];

                    // this is dy
                    detail[1] += lsdata[m_sidx(0, 2 * dm0 + 0, 0, 0, ctx->srcstr)] * fact;
                    detail[1] += lsdata[m_sidx(1, 2 * dm0 + 0, 0, 0, ctx->srcstr)] * fact;
                    detail[1] += lsdata[m_sidx(0, 2 * dm0 + 1, 0, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[1] += lsdata[m_sidx(0, 2 * dm0 + 0, 1, 0, ctx->srcstr)] * fact;
                    detail[1] += lsdata[m_sidx(0, 2 * dm0 + 1, 1, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[1] += lsdata[m_sidx(1, 2 * dm0 + 1, 0, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[1] += lsdata[m_sidx(1, 2 * dm0 + 0, 1, 0, ctx->srcstr)] * fact;
                    detail[1] += lsdata[m_sidx(1, 2 * dm0 + 1, 1, 0, ctx->srcstr)] * fact * sign[dm0];

                    // this is dz
                    detail[2] += lsdata[m_sidx(0, 0, 2 * dm0 + 0, 0, ctx->srcstr)] * fact;
                    detail[2] += lsdata[m_sidx(1, 0, 2 * dm0 + 0, 0, ctx->srcstr)] * fact;
                    detail[2] += lsdata[m_sidx(0, 1, 2 * dm0 + 0, 0, ctx->srcstr)] * fact;
                    detail[2] += lsdata[m_sidx(0, 0, 2 * dm0 + 1, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[2] += lsdata[m_sidx(0, 1, 2 * dm0 + 1, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[2] += lsdata[m_sidx(1, 1, 2 * dm0 + 0, 0, ctx->srcstr)] * fact;
                    detail[2] += lsdata[m_sidx(1, 0, 2 * dm0 + 1, 0, ctx->srcstr)] * fact * sign[dm0];
                    detail[2] += lsdata[m_sidx(1, 1, 2 * dm0 + 1, 0, ctx->srcstr)] * fact * sign[dm0];

                    // m_verb("order = %d, dm0 = %d: fact = %f, sign = %f -> dx = %e, dy = %e, dz = %e",order,dm0,fact,sign[dm0],detail[0],detail[1],detail[2]);
                }
                // m_verb("order = %d, =================> dx = %e, dy = %e, dz = %e",order,detail[0],detail[1],detail[2]);

                detail[3] = 0.0;
                detail[4] = 0.0;
                detail[5] = 0.0;
                for (int dm1 = -(order / 2); dm1 <= (order / 2); dm1++) {
                    for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
                        const real_t fact = ha[0] * ga[dm1] * ga[dm0];
                        // this is dxy
                        detail[3] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 0, 0, ctx->srcstr)] * fact;
                        detail[3] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 0, 0, ctx->srcstr)] * fact * sign[dm0];
                        detail[3] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 0, 0, ctx->srcstr)] * fact * sign[dm1];
                        detail[3] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 1, 0, ctx->srcstr)] * fact;
                        detail[3] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 1, 0, ctx->srcstr)] * fact * sign[dm1];
                        detail[3] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 0, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];
                        detail[3] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 1, 0, ctx->srcstr)] * fact * sign[dm0];
                        detail[3] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];

                        // this is dyz
                        detail[4] += lsdata[m_sidx(0, 2 * dm0 + 0, 2 * dm1 + 0, 0, ctx->srcstr)] * fact;
                        detail[4] += lsdata[m_sidx(1, 2 * dm0 + 0, 2 * dm1 + 0, 0, ctx->srcstr)] * fact;
                        detail[4] += lsdata[m_sidx(0, 2 * dm0 + 1, 2 * dm1 + 0, 0, ctx->srcstr)] * fact * sign[dm0];
                        detail[4] += lsdata[m_sidx(0, 2 * dm0 + 0, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm1];
                        detail[4] += lsdata[m_sidx(0, 2 * dm0 + 1, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];
                        detail[4] += lsdata[m_sidx(1, 2 * dm0 + 1, 2 * dm1 + 0, 0, ctx->srcstr)] * fact * sign[dm0];
                        detail[4] += lsdata[m_sidx(1, 2 * dm0 + 0, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm1];
                        detail[4] += lsdata[m_sidx(1, 2 * dm0 + 1, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];

                        // this is dxz
                        detail[5] += lsdata[m_sidx(2 * dm0 + 0, 0, 2 * dm1 + 0, 0, ctx->srcstr)] * fact;
                        detail[5] += lsdata[m_sidx(2 * dm0 + 1, 0, 2 * dm1 + 0, 0, ctx->srcstr)] * fact * sign[dm0];
                        detail[5] += lsdata[m_sidx(2 * dm0 + 0, 1, 2 * dm1 + 0, 0, ctx->srcstr)] * fact;
                        detail[5] += lsdata[m_sidx(2 * dm0 + 0, 0, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm1];
                        detail[5] += lsdata[m_sidx(2 * dm0 + 0, 1, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm1];
                        detail[5] += lsdata[m_sidx(2 * dm0 + 1, 1, 2 * dm1 + 0, 0, ctx->srcstr)] * fact * sign[dm0];
                        detail[5] += lsdata[m_sidx(2 * dm0 + 1, 0, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];
                        detail[5] += lsdata[m_sidx(2 * dm0 + 1, 1, 2 * dm1 + 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];
                    }
                }

                detail[6] = 0.0;
                for (int dm2 = -(order / 2); dm2 <= (order / 2); dm2++) {
                    for (int dm1 = -(order / 2); dm1 <= (order / 2); dm1++) {
                        for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
                            const real_t fact = ha[0] * ga[dm1] * ga[dm0];
                            // this is dxy
                            detail[6] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 2 * dm2 + 0, 0, ctx->srcstr)] * fact;
                            detail[6] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 2 * dm2 + 0, 0, ctx->srcstr)] * fact * sign[dm0];
                            detail[6] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 2 * dm2 + 0, 0, ctx->srcstr)] * fact * sign[dm1];
                            detail[6] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 2 * dm2 + 1, 0, ctx->srcstr)] * fact * sign[dm2];
                            detail[6] += lsdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 2 * dm2 + 1, 0, ctx->srcstr)] * fact * sign[dm1] * sign[dm2];
                            detail[6] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 2 * dm2 + 0, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1];
                            detail[6] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 2 * dm2 + 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm2];
                            detail[6] += lsdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 2 * dm2 + 1, 0, ctx->srcstr)] * fact * sign[dm0] * sign[dm1] * sign[dm2];
                        }
                    }
                }

                // get the mean value, the last missing ingredient
                detail[7] = 0.0;
                detail[7] += lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(1, 0, 0, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(0, 1, 0, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(0, 0, 1, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(0, 1, 1, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(1, 1, 0, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(1, 0, 1, 0, ctx->srcstr)] * 0.125;
                detail[7] += lsdata[m_sidx(1, 1, 1, 0, ctx->srcstr)] * 0.125;

                // m_verb("order = %d my details = %f %f %f %f %f %f %f %f",order,detail[0],detail[1],detail[2],detail[3],detail[4],detail[5],detail[6],detail[7]);

                // store the max, relatively to the current value
                for (int id = 0; id < 8; id++) {
                    details_inf_norm[id] = m_max(std::fabs(detail[id]), details_inf_norm[id]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

#endif