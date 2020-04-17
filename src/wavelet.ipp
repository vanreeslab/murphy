#ifndef SRC_WAVELET_IPP_
#define SRC_WAVELET_IPP_

#include <cmath>

#include "murphy.hpp"
#include "wavelet.hpp"

using std::pow;

template <int order>
void Wavelet<order>::Criterion(MemLayout* block, real_p data, real_t* criterion) {
    //-------------------------------------------------------------------------
    interp_ctx_t ctx;
    real_t details_max[7] = {0};
    // get memory details
    for (int id = 0; id < 3; id++) {
        ctx.srcstart[id] = block->start(id);
        ctx.srcend[id]   = block->end(id);
        ctx.trgstart[id] = -1;
        ctx.trgend[id] = -2;
    }
    ctx.srcstr = block->stride();
    ctx.sdata = data;
    ctx.trgstr = -1000;
    ctx.tdata = nullptr;


    Detail_(&ctx,details_max);

    // get the max out of all the details
    for (int id = 0; id < 7; id++) {
        (*criterion) = m_max(*criterion, details_max[id]);
    }
    //-------------------------------------------------------------------------
}

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

template <int order>
void Wavelet<order>::Detail_(const interp_ctx_t* ctx, real_t* details_inf_norm) const {
    m_assert(!(order==5 && M_GS<4), "the detail computation requires at least 4 ghost points on the current block");
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
                real_t detail[7] = {0};

                // compute all the single coefficients
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
                }

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

                // store the max
                for (int id = 0; id < 7; id++) {
                    details_inf_norm[id] = m_max(std::fabs(detail[id]), details_inf_norm[id]);
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

#endif