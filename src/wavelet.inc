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
real_t Wavelet<N,Nt>::Criterion(MemLayout* block, real_p data, MemPool* mem_pool) {
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
                // for every block of 8 child data, get how much information I can get from my parent
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
void Wavelet<N,Nt>::Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const {
    //-------------------------------------------------------------------------
    // assure alignment for the target, the source, the constant and the temp data
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    m_assume_aligned(ctx->cdata);
    m_assume_aligned(ctx->wdata);

    const lid_t  n_ghost = this->NGhostFine();
    const real_t alpha   = ctx->alpha;

    //copy the source data to the working memory (so we can change it's values)
    for (lid_t ik2 = ctx->trgstart[2] - n_ghost; ik2 < ctx->trgend[2] + n_ghost; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1] - n_ghost; ik1 < ctx->trgend[1] + n_ghost; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0] - n_ghost; ik0 < ctx->trgend[0] + n_ghost; ik0++) {
                // check memory access
                m_assert(((2 * dlvl) * ik0 >= ctx->srcstart[0]) && ((2 * dlvl) * ik0 + (2 * dlvl - 1) < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", (2 * dlvl) * ik0, ctx->srcstart[0], (2 * dlvl) * ik0 + (2 * dlvl - 1), ctx->srcend[0]);
                m_assert(((2 * dlvl) * ik1 >= ctx->srcstart[1]) && ((2 * dlvl) * ik1 + (2 * dlvl - 1) < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", (2 * dlvl) * ik1, ctx->srcstart[1], (2 * dlvl) * ik1 + (2 * dlvl - 1), ctx->srcend[1]);
                m_assert(((2 * dlvl) * ik2 >= ctx->srcstart[2]) && ((2 * dlvl) * ik2 + (2 * dlvl - 1) < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", (2 * dlvl) * ik2, ctx->srcstart[2], (2 * dlvl) * ik2 + (2 * dlvl - 1), ctx->srcend[2]);
                // do the copy of the 8 associated elements
                real_p       lwdata = ctx->wdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);
                const real_p lsdata = ctx->sdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);
                m_assume_aligned(lwdata);
                m_assume_aligned(lsdata);
                // do the copy
                lwdata[m_sidx(0, 0, 0, 0, ctx->srcstr)] = lsdata[m_sidx(0, 0, 0, 0, ctx->srcstr)];
                lwdata[m_sidx(1, 0, 0, 0, ctx->srcstr)] = lsdata[m_sidx(1, 0, 0, 0, ctx->srcstr)];
                lwdata[m_sidx(0, 1, 0, 0, ctx->srcstr)] = lsdata[m_sidx(0, 1, 0, 0, ctx->srcstr)];
                lwdata[m_sidx(0, 0, 1, 0, ctx->srcstr)] = lsdata[m_sidx(0, 0, 1, 0, ctx->srcstr)];
                lwdata[m_sidx(1, 1, 0, 0, ctx->srcstr)] = lsdata[m_sidx(1, 1, 0, 0, ctx->srcstr)];
                lwdata[m_sidx(1, 0, 1, 0, ctx->srcstr)] = lsdata[m_sidx(1, 0, 1, 0, ctx->srcstr)];
                lwdata[m_sidx(0, 1, 1, 0, ctx->srcstr)] = lsdata[m_sidx(0, 1, 1, 0, ctx->srcstr)];
                lwdata[m_sidx(1, 1, 1, 0, ctx->srcstr)] = lsdata[m_sidx(1, 1, 1, 0, ctx->srcstr)];
            }
        }
    }

    for (sid_t il = 0; il < dlvl; il++) {
        //-----------------
        // DUAL-LIFTING
        // loop over the scaling coef and give the contribution to every detail
        for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
            for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
                for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                    real_p lwdata = ctx->wdata + m_sidx(2 * il * ik0, 2 * il * ik1, 2 * il * ik2, 0, ctx->srcstr);
                    // save the current scaling value
                    real_t scaling = lwdata[0];
                    // update the details associated with the current scaling
                    for (sid_t id2 = -(N / 2); id2 < (N / 2); id2++) {
                        for (sid_t id1 = -(N / 2); id1 < (N / 2); id1++) {
                            for (sid_t id0 = -(N / 2); id0 < (N / 2); id0++) {
                                lwdata[m_sidx(id0 * il, id1 * il, id2 * il, 0, ctx->srcstr)] -= scaling * s_dual_[id0] * s_dual_[id1] * s_dual_[id2];
                            }
                        }
                    }
                    // the scaling has been updated, so we reset its value
                    lwdata[0] = scaling;
                }
            }
        }
        //-----------------
        // LIFTING
        // loop over the scaling and get the contribution from every detail
        for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
            for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
                for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                    real_p lwdata = ctx->wdata + m_sidx(2 * il * ik0, 2 * il * ik1, 2 * il * ik2, 0, ctx->srcstr);
                    // save the current scaling value
                    real_t scaling = lwdata[0];
                    // update the scaling with the neighboring details
                    for (sid_t id2 = -(Nt / 2); id2 < (Nt / 2); id2++) {
                        for (sid_t id1 = -(Nt / 2); id1 < (Nt / 2); id1++) {
                            for (sid_t id0 = -(Nt / 2); id0 < (Nt / 2); id0++) {
                                lwdata[0] += lwdata[m_sidx(id0 * il, id1 * il, id2 * il, 0, ctx->srcstr)] * s_lift_[id0] * s_lift_[id1] * s_lift_[id2];
                            }
                        }
                    }
                    // we added the scaling twice, so removing it
                    lwdata[0] -= lwdata[0] * s_lift_[0] * s_lift_[0] * s_lift_[0];
                }
            }
        }
    }

    //-----------------
    // SPLITTING - we just take one data out of 2
    // note: this cannot be merged with the lifting since we apply the splitting recursivelly to coarsen multiple levels at a time
    for (lid_t ik2 = ctx->trgstart[2]; ik2 < ctx->trgend[2]; ik2++) {
        for (lid_t ik1 = ctx->trgstart[1]; ik1 < ctx->trgend[1]; ik1++) {
            for (lid_t ik0 = ctx->trgstart[0]; ik0 < ctx->trgend[0]; ik0++) {
                //get the local adress of the source, the target and the constant
                real_p       ltdata = ctx->tdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lcdata = ctx->cdata + m_sidx(ik0, ik1, ik2, 0, ctx->trgstr);
                const real_p lwdata = ctx->wdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->srcstr);
                // fill the new scaling value
                ltdata[0] = alpha * lcdata[0] + lwdata[0];
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
    // assure alignment for both target and source
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);
    m_assume_aligned(ctx->cdata);

    // const lid_t   hslen = order / 2;
    // const real_t* hs    = hs_ + hslen;
    // const real_t* sign  = sgn_hs_ + hslen;
    // const real_t  alpha = ctx->alpha;

    // // for each of the data for the considered children
    // for (int ik2 = ctx->trgstart[2] / 2; ik2 < ctx->trgend[2] / 2; ik2++) {
    //     for (int ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
    //         for (int ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
    //             m_assert((ik0 - hslen >= ctx->srcstart[0]) && (ik0 + hslen < ctx->srcend[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - hslen, ctx->srcstart[0], ik0 + hslen, ctx->srcend[0]);
    //             m_assert((ik1 - hslen >= ctx->srcstart[1]) && (ik1 + hslen < ctx->srcend[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - hslen, ctx->srcstart[1], ik1 + hslen, ctx->srcend[1]);
    //             m_assert((ik2 - hslen >= ctx->srcstart[2]) && (ik2 + hslen < ctx->srcend[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - hslen, ctx->srcstart[2], ik2 + hslen, ctx->srcend[2]);
    //             //get the local adress of the dady
    //             real_p       ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
    //             const real_p lcdata = ctx->cdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
    //             const real_p lsdata = ctx->sdata + m_sidx(ik0, ik1, ik2, 0, ctx->srcstr);

    //             // reset to 0.0
    //             ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 0, 0, ctx->trgstr)];
    //             ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 0, 0, ctx->trgstr)];
    //             ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 0, 0, ctx->trgstr)];
    //             ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 0, 1, 0, ctx->trgstr)];
    //             ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(0, 1, 1, 0, ctx->trgstr)];
    //             ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 0, 0, ctx->trgstr)];
    //             ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 0, 1, 0, ctx->trgstr)];
    //             ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] = alpha * lcdata[m_sidx(1, 1, 1, 0, ctx->trgstr)];

    //             for (int dm2 = -(order / 2); dm2 <= (order / 2); dm2++) {
    //                 for (int dm1 = -(order / 2); dm1 <= (order / 2); dm1++) {
    //                     for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
    //                         // get the current parent's data
    //                         const real_t ldata = lsdata[m_sidx(dm0, dm1, dm2, 0, ctx->srcstr)];
    //                         const real_t fact  = hs[dm2] * hs[dm1] * hs[dm0];
    //                         // we give the information to every wavelet which is inside my "block"
    //                         ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] += fact * ldata;
    //                         ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] += fact * ldata * sign[dm0];
    //                         ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] += fact * ldata * sign[dm1];
    //                         ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm2];
    //                         ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm1] * sign[dm2];
    //                         ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] += fact * ldata * sign[dm0] * sign[dm1];
    //                         ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm0] * sign[dm2];
    //                         ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] += fact * ldata * sign[dm0] * sign[dm1] * sign[dm2];
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
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
    // assure alignment for both target and source
    m_assume_aligned(ctx->tdata);

    // const lid_t   hslen = order / 2;
    // const real_t* ga    = ga_ + hslen;
    // const real_t* sign  = sgn_ga_ + hslen;

    // real_t ha_tmp[order] = {0};
    // ha_tmp[order / 2]    = 0.5;
    // const real_t* ha     = ha_tmp + hslen;

    // // for each of the data for the considered children
    // for (int ik2 = ctx->trgstart[2] / 2; ik2 < ctx->trgend[2] / 2; ik2++) {
    //     for (int ik1 = ctx->trgstart[1] / 2; ik1 < ctx->trgend[1] / 2; ik1++) {
    //         for (int ik0 = ctx->trgstart[0] / 2; ik0 < ctx->trgend[0] / 2; ik0++) {
    //             //get the local adress of the dady
    //             real_p ltdata = ctx->tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, ctx->trgstr);
    //             m_assume_aligned(ltdata);

    //             // set the detail coefficient to the
    //             real_t detail[8];

    //             // compute all the single coefficients
    //             detail[0] = 0.0;
    //             detail[1] = 0.0;
    //             detail[2] = 0.0;
    //             for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
    //                 const real_t fact = ha[0] * ha[0] * ga[dm0];
    //                 // this is dx
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 0, 0, 0, 0, ctx->trgstr)] * fact;
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 1, 0, 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 0, 1, 0, 0, ctx->trgstr)] * fact;
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 0, 0, 1, 0, ctx->trgstr)] * fact;
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 0, 1, 1, 0, ctx->trgstr)] * fact;
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 1, 1, 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 1, 0, 1, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[0] += ltdata[m_sidx(2 * dm0 + 1, 1, 1, 0, ctx->trgstr)] * fact * sign[dm0];

    //                 // this is dy
    //                 detail[1] += ltdata[m_sidx(0, 2 * dm0 + 0, 0, 0, ctx->trgstr)] * fact;
    //                 detail[1] += ltdata[m_sidx(1, 2 * dm0 + 0, 0, 0, ctx->trgstr)] * fact;
    //                 detail[1] += ltdata[m_sidx(0, 2 * dm0 + 1, 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[1] += ltdata[m_sidx(0, 2 * dm0 + 0, 1, 0, ctx->trgstr)] * fact;
    //                 detail[1] += ltdata[m_sidx(0, 2 * dm0 + 1, 1, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[1] += ltdata[m_sidx(1, 2 * dm0 + 1, 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[1] += ltdata[m_sidx(1, 2 * dm0 + 0, 1, 0, ctx->trgstr)] * fact;
    //                 detail[1] += ltdata[m_sidx(1, 2 * dm0 + 1, 1, 0, ctx->trgstr)] * fact * sign[dm0];

    //                 // this is dz
    //                 detail[2] += ltdata[m_sidx(0, 0, 2 * dm0 + 0, 0, ctx->trgstr)] * fact;
    //                 detail[2] += ltdata[m_sidx(1, 0, 2 * dm0 + 0, 0, ctx->trgstr)] * fact;
    //                 detail[2] += ltdata[m_sidx(0, 1, 2 * dm0 + 0, 0, ctx->trgstr)] * fact;
    //                 detail[2] += ltdata[m_sidx(0, 0, 2 * dm0 + 1, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[2] += ltdata[m_sidx(0, 1, 2 * dm0 + 1, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[2] += ltdata[m_sidx(1, 1, 2 * dm0 + 0, 0, ctx->trgstr)] * fact;
    //                 detail[2] += ltdata[m_sidx(1, 0, 2 * dm0 + 1, 0, ctx->trgstr)] * fact * sign[dm0];
    //                 detail[2] += ltdata[m_sidx(1, 1, 2 * dm0 + 1, 0, ctx->trgstr)] * fact * sign[dm0];

    //                 // m_verb("order = %d, dm0 = %d: fact = %f, sign = %f -> dx = %e, dy = %e, dz = %e",order,dm0,fact,sign[dm0],detail[0],detail[1],detail[2]);
    //             }
    //             // m_verb("order = %d, =================> dx = %e, dy = %e, dz = %e",order,detail[0],detail[1],detail[2]);

    //             detail[3] = 0.0;
    //             detail[4] = 0.0;
    //             detail[5] = 0.0;
    //             for (int dm1 = -(order / 2); dm1 <= (order / 2); dm1++) {
    //                 for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
    //                     const real_t fact = ha[0] * ga[dm1] * ga[dm0];
    //                     // this is dxy
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 0, 0, ctx->trgstr)] * fact;
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 0, 0, ctx->trgstr)] * fact * sign[dm1];
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 1, 0, ctx->trgstr)] * fact;
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 1, 0, ctx->trgstr)] * fact * sign[dm1];
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 0, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 1, 0, ctx->trgstr)] * fact * sign[dm0];
    //                     detail[3] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];

    //                     // this is dyz
    //                     detail[4] += ltdata[m_sidx(0, 2 * dm0 + 0, 2 * dm1 + 0, 0, ctx->trgstr)] * fact;
    //                     detail[4] += ltdata[m_sidx(1, 2 * dm0 + 0, 2 * dm1 + 0, 0, ctx->trgstr)] * fact;
    //                     detail[4] += ltdata[m_sidx(0, 2 * dm0 + 1, 2 * dm1 + 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                     detail[4] += ltdata[m_sidx(0, 2 * dm0 + 0, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm1];
    //                     detail[4] += ltdata[m_sidx(0, 2 * dm0 + 1, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];
    //                     detail[4] += ltdata[m_sidx(1, 2 * dm0 + 1, 2 * dm1 + 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                     detail[4] += ltdata[m_sidx(1, 2 * dm0 + 0, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm1];
    //                     detail[4] += ltdata[m_sidx(1, 2 * dm0 + 1, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];

    //                     // this is dxz
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 0, 0, 2 * dm1 + 0, 0, ctx->trgstr)] * fact;
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 1, 0, 2 * dm1 + 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 0, 1, 2 * dm1 + 0, 0, ctx->trgstr)] * fact;
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 0, 0, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm1];
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 0, 1, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm1];
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 1, 1, 2 * dm1 + 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 1, 0, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];
    //                     detail[5] += ltdata[m_sidx(2 * dm0 + 1, 1, 2 * dm1 + 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];
    //                 }
    //             }

    //             detail[6] = 0.0;
    //             for (int dm2 = -(order / 2); dm2 <= (order / 2); dm2++) {
    //                 for (int dm1 = -(order / 2); dm1 <= (order / 2); dm1++) {
    //                     for (int dm0 = -(order / 2); dm0 <= (order / 2); dm0++) {
    //                         const real_t fact = ha[0] * ga[dm1] * ga[dm0];
    //                         // this is dxy
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 2 * dm2 + 0, 0, ctx->trgstr)] * fact;
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 2 * dm2 + 0, 0, ctx->trgstr)] * fact * sign[dm0];
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 2 * dm2 + 0, 0, ctx->trgstr)] * fact * sign[dm1];
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 0, 2 * dm2 + 1, 0, ctx->trgstr)] * fact * sign[dm2];
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 0, 2 * dm1 + 1, 2 * dm2 + 1, 0, ctx->trgstr)] * fact * sign[dm1] * sign[dm2];
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 2 * dm2 + 0, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1];
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 0, 2 * dm2 + 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm2];
    //                         detail[6] += ltdata[m_sidx(2 * dm0 + 1, 2 * dm1 + 1, 2 * dm2 + 1, 0, ctx->trgstr)] * fact * sign[dm0] * sign[dm1] * sign[dm2];
    //                     }
    //                 }
    //             }

    //             // get the mean value, the last missing ingredient
    //             detail[7] = 0.0;
    //             detail[7] += ltdata[m_sidx(0, 0, 0, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(1, 0, 0, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(0, 1, 0, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(0, 0, 1, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(0, 1, 1, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(1, 1, 0, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(1, 0, 1, 0, ctx->trgstr)] * 0.125;
    //             detail[7] += ltdata[m_sidx(1, 1, 1, 0, ctx->trgstr)] * 0.125;

    //             // m_verb("order = %d my details = %f %f %f %f %f %f %f %f",order,detail[0],detail[1],detail[2],detail[3],detail[4],detail[5],detail[6],detail[7]);

    //             // store the max, relatively to the current value
    //             for (int id = 0; id < 8; id++) {
    //                 details_inf_norm[id] = m_max(std::fabs(detail[id]), details_inf_norm[id]);
    //             }
    //         }
    //     }
    // }
    //-------------------------------------------------------------------------
}

#endif