#include <cmath>

using std::pow;

template <int order>
void Wavelet<order>::Copy_(const interp_ctx_t* ctx) const {
    m_begin;
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
    m_end;
}
template <int order>
void Wavelet<order>::Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const {
    m_begin;
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
    m_end;
}
template <int order>
void Wavelet<order>::Refine_(const interp_ctx_t* ctx) const {
    m_begin;
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(ctx->tdata);
    m_assume_aligned(ctx->sdata);

    const lid_t   hslen = order / 2;
    const real_t* hs    = hs_ + hslen;
    const real_t* sign  = sgn_ + hslen;

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
    m_end;
}
