#include "wavelet.hpp"

#include <cmath>

using std::pow;

void Wavelet::Copy_(const real_p sdata, real_p tdata) const {
    m_begin;
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(tdata);
    m_assume_aligned(sdata);

    // for each of the data for the considered children
    for (int ik2 = trgstart_[2]; ik2 < trgend_[2]; ik2++) {
        for (int ik1 = trgstart_[1]; ik1 < trgend_[1]; ik1++) {
            for (int ik0 = trgstart_[0]; ik0 < trgend_[0]; ik0++) {
                // for every block of 8 child data, get how much information I can get from my parent
                m_assert((ik0 >= srcstart_[0]) && (ik0 < srcend_[0]), "the source domain is too small in dir 0: %d >= %d and %d<%d", ik0, srcstart_[0], ik0, srcend_[0]);
                m_assert((ik1 >= srcstart_[1]) && (ik1 < srcend_[1]), "the source domain is too small in dir 1: %d >= %d and %d<%d", ik1, srcstart_[1], ik1, srcend_[1]);
                m_assert((ik2 >= srcstart_[2]) && (ik2 < srcend_[2]), "the source domain is too small in dir 2: %d >= %d and %d<%d", ik2, srcstart_[2], ik2, srcend_[2]);
                // get the current parent's data
                tdata[m_sidx(ik0, ik1, ik2, 0, trgstr_, 0)] = sdata[m_sidx(ik0, ik1, ik2, 0, srcstr_, 0)];
            }
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Wavelet::Coarsen_(const lid_t dlvl, const real_p sdata, real_p tdata) const {
    m_begin;
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(tdata);
    m_assume_aligned(sdata);

    // get the factor = 0.125^(#d level)
    const real_t fact = pow(0.125, dlvl);
    // for each of the data for the considered children
    for (int ik2 = trgstart_[2]; ik2 < trgend_[2]; ik2++) {
        for (int ik1 = trgstart_[1]; ik1 < trgend_[1]; ik1++) {
            for (int ik0 = trgstart_[0]; ik0 < trgend_[0]; ik0++) {
                m_assert(((2 * dlvl) * ik0 >= srcstart_[0]) && ((2 * dlvl) * ik0 + (2 * dlvl - 1) < srcend_[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", (2 * dlvl) * ik0, srcstart_[0], (2 * dlvl) * ik0 + (2 * dlvl - 1), srcend_[0]);
                m_assert(((2 * dlvl) * ik1 >= srcstart_[1]) && ((2 * dlvl) * ik1 + (2 * dlvl - 1) < srcend_[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", (2 * dlvl) * ik1, srcstart_[1], (2 * dlvl) * ik1 + (2 * dlvl - 1), srcend_[1]);
                m_assert(((2 * dlvl) * ik2 >= srcstart_[2]) && ((2 * dlvl) * ik2 + (2 * dlvl - 1) < srcend_[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", (2 * dlvl) * ik2, srcstart_[2], (2 * dlvl) * ik2 + (2 * dlvl - 1), srcend_[2]);
                //get the local adress of the dady
                real_p lsdata = sdata + m_sidx((2 * dlvl) * ik0, (2 * dlvl) * ik1, (2 * dlvl) * ik2, 0, srcstr_, 0);
                m_assume_aligned(lsdata);
                const real_p ltdata = tdata + m_sidx(ik0, ik1, ik2, 0, trgstr_, 0);

                // do the mean on the current patch
                ltdata[0] = 0.0;
                for (int it2 = 0; it2 < 2 * dlvl; it2++) {
                    for (int it1 = 0; it1 < 2 * dlvl; it1++) {
                        for (int it0 = 0; it0 < 2 * dlvl; it0++) {
                            ltdata[0] += fact * lsdata[m_sidx(it0, it1, it2, 0, srcstr_, 0)];
                        }
                    }
                }
            }
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}

void Wavelet::Refine_(const real_p sdata, real_p tdata) const {
    m_begin;
    //-------------------------------------------------------------------------
    // assure alignment for both target and source
    m_assume_aligned(tdata);
    m_assume_aligned(sdata);

    const lid_t   hslen = hslen_ / 2;
    const real_t* hs    = hs_ + hslen;
    const real_t* sign  = sgn_ + hslen;

    // for each of the data for the considered children
    for (int ik2 = trgstart_[2] / 2; ik2 < trgend_[2] / 2; ik2++) {
        for (int ik1 = trgstart_[1] / 2; ik1 < trgend_[1] / 2; ik1++) {
            for (int ik0 = trgstart_[0] / 2; ik0 < trgend_[0] / 2; ik0++) {
                m_assert((ik0 - hslen >= srcstart_[0]) && (ik0 + hslen < srcend_[0]), "the source domain is too small in dir 0: %d >= %d and %d < %d", ik0 - hslen, srcstart_[0], ik0 + hslen, srcend_[0]);
                m_assert((ik1 - hslen >= srcstart_[1]) && (ik1 + hslen < srcend_[1]), "the source domain is too small in dir 1: %d >= %d and %d < %d", ik1 - hslen, srcstart_[1], ik1 + hslen, srcend_[1]);
                m_assert((ik2 - hslen >= srcstart_[2]) && (ik2 + hslen < srcend_[2]), "the source domain is too small in dir 2: %d >= %d and %d < %d", ik2 - hslen, srcstart_[2], ik2 + hslen, srcend_[2]);
                //get the local adress of the dady
                real_p ltdata = tdata + m_sidx(2 * ik0, 2 * ik1, 2 * ik2, 0, trgstr_, 0);
                m_assume_aligned(ltdata);
                const real_p lsdata = sdata + m_sidx(ik0, ik1, ik2, 0, srcstr_, 0);

                for (int dm2 = -hslen; dm2 <= hslen; dm2++) {
                    for (int dm1 = -hslen; dm1 <= hslen; dm1++) {
                        for (int dm0 = -hslen; dm0 <= hslen; dm0++) {
                            // get the current parent's data
                            const real_t ldata = lsdata[m_sidx(dm0, dm1, dm2, 0, srcstr_, 0)];
                            const real_t fact  = hs[dm2] * hs[dm1] * hs[dm0];
                            // we give the information to every wavelet which is inside my "block"
                            ltdata[m_sidx(0, 0, 0, 0, trgstr_, 0)] += fact * ldata;
                            ltdata[m_sidx(1, 0, 0, 0, trgstr_, 0)] += fact * ldata * sign[dm0];
                            ltdata[m_sidx(0, 1, 0, 0, trgstr_, 0)] += fact * ldata * sign[dm1];
                            ltdata[m_sidx(0, 0, 1, 0, trgstr_, 0)] += fact * ldata * sign[dm2];
                            ltdata[m_sidx(0, 1, 1, 0, trgstr_, 0)] += fact * ldata * sign[dm1] * sign[dm2];
                            ltdata[m_sidx(1, 1, 0, 0, trgstr_, 0)] += fact * ldata * sign[dm0] * sign[dm1];
                            ltdata[m_sidx(1, 0, 1, 0, trgstr_, 0)] += fact * ldata * sign[dm0] * sign[dm2];
                            ltdata[m_sidx(1, 1, 1, 0, trgstr_, 0)] += fact * ldata * sign[dm0] * sign[dm1] * sign[dm2];
                        }
                    }
                }
            }
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
