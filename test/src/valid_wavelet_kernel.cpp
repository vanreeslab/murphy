
#include "error.hpp"
#include "gtest/gtest.h"
#include "murphy.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-13

// cf Abramowitz p887, formula 25.4.20
static real_t Integral(const lid_t N, const real_t* f, const real_t h) {
    m_assert((N % 10) == 1, "the size MUST be 11 ");

    real_t sum = 0.0;
    for (lid_t id = 0; id < (N % 10); ++id) {
        lid_t cid = id * 10;
        sum += (5 * h) / (299376) *
               (16067 * (f[cid + 0] + f[cid + 10]) +
                106300 * (f[cid + 1] + f[cid + 9]) -
                48525 * (f[cid + 2] + f[cid + 8]) +
                272400 * (f[cid + 3] + f[cid + 7]) -
                260550 * (f[cid + 4] + f[cid + 6]) +
                427368 * f[cid + 5]);
    }
    return sum;
};

static double poly(const double x) {
    if constexpr (M_WAVELET_N == 2) {
        return M_PI * x + M_PI_2;
    } else if (M_WAVELET_N == 4) {
        return M_PI * pow(x, 3) + M_PI_2 * pow(x, 2) + M_SQRT2 * x + M_PI_4;
    } else {
        return M_PI * pow(x, M_WAVELET_N - 1) + M_PI_2 * pow(x, M_WAVELET_N - 2) + M_SQRT2 * pow(x, M_WAVELET_N - 3) + M_PI_4;
    }
}

static constexpr lid_t m_gs     = M_GS;
static constexpr lid_t m_n      = 48;
static constexpr lid_t m_hn     = 24;
static constexpr lid_t m_stride = m_n + 2 * m_gs;

/**
 * @brief test the wavelets on blocks of length 12 + 6x2 gp = 24
 * 
 */
class valid_Wavelet_Kernel : public ::testing::Test {
   protected:
    SubBlock* block_coarse_;
    SubBlock* block_fine_;
    real_p    data_coarse_;
    real_p    data_fine_;
    real_t    hcoarse_;
    real_t    hfine_;

    lid_t coarse_start_[3];
    lid_t coarse_end_[3];
    lid_t fine_start_[3];
    lid_t fine_end_[3];

    void SetUp() override {
        data_fine_   = (real_t*)m_calloc(m_stride * m_stride * m_stride * sizeof(real_t));
        data_coarse_ = (real_t*)m_calloc(m_stride * m_stride * m_stride * sizeof(real_t));

        for (int id = 0; id < 3; id++) {
            coarse_start_[id] = -m_gs;
            coarse_end_[id]   = m_hn + m_gs;
            fine_start_[id]   = 0;
            fine_end_[id]     = m_n;
        }
        block_coarse_ = new SubBlock(m_gs, m_hn + 2 * m_gs, coarse_start_, coarse_end_);
        block_fine_   = new SubBlock(m_gs, m_stride, fine_start_, fine_end_);

        hcoarse_ = 1.0 / (m_hn);
        hfine_   = 1.0 / (m_n);
    };
    void TearDown() override {
        m_free(data_fine_);
        m_free(data_coarse_);

        delete (block_coarse_);
        delete (block_fine_);
    };
};

TEST_F(valid_Wavelet_Kernel, filter_lentgh) {
    InterpolatingWavelet interp;

    if (interp.N() == 2 && interp.Nt() == 2) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 2);
        ASSERT_EQ(interp.ncoarsen_back(), 1);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 0);
        ASSERT_EQ(interp.nrefine_back(), 1);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 2);
        ASSERT_EQ(interp.ncriterion_back(), 1);
        // shift
        ASSERT_EQ(interp.shift_front(), 1);
        ASSERT_EQ(interp.shift_back(), 0);
    } else if (interp.N() == 4 && interp.Nt() == 0) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 0);
        ASSERT_EQ(interp.ncoarsen_back(), 0);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 1);
        ASSERT_EQ(interp.nrefine_back(), 2);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 2);
        ASSERT_EQ(interp.ncriterion_back(), 3);
        // shift
        ASSERT_EQ(interp.shift_front(), 0);
        ASSERT_EQ(interp.shift_back(), 0);
    } else if (interp.N() == 4 && interp.Nt() == 2) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 4);
        ASSERT_EQ(interp.ncoarsen_back(), 3);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 1);
        ASSERT_EQ(interp.nrefine_back(), 2);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 4);
        ASSERT_EQ(interp.ncriterion_back(), 3);
        // shift
        ASSERT_EQ(interp.shift_front(), 1);
        ASSERT_EQ(interp.shift_back(), 0);
    } else if (interp.N() == 4 && interp.Nt() == 4) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 6);
        ASSERT_EQ(interp.ncoarsen_back(), 5);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 1);
        ASSERT_EQ(interp.nrefine_back(), 2);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 6);
        ASSERT_EQ(interp.ncriterion_back(), 5);
        // shift
        ASSERT_EQ(interp.shift_front(), 3);
        ASSERT_EQ(interp.shift_back(), 2);
    } else if (interp.N() == 6 && interp.Nt() == 0) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 0);
        ASSERT_EQ(interp.ncoarsen_back(), 0);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 2);
        ASSERT_EQ(interp.nrefine_back(), 3);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 4);
        ASSERT_EQ(interp.ncriterion_back(), 5);
        // shift
        ASSERT_EQ(interp.shift_front(), 0);
        ASSERT_EQ(interp.shift_back(), 0);
    } else if (interp.N() == 6 && interp.Nt() == 2) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 6);
        ASSERT_EQ(interp.ncoarsen_back(), 5);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 2);
        ASSERT_EQ(interp.nrefine_back(), 3);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 6);
        ASSERT_EQ(interp.ncriterion_back(), 5);
        // shift
        ASSERT_EQ(interp.shift_front(), 1);
        ASSERT_EQ(interp.shift_back(), 0);
    } else if (interp.N() == 6 && interp.Nt() == 4) {
        // coarsen
        ASSERT_EQ(interp.ncoarsen_front(), 8);
        ASSERT_EQ(interp.ncoarsen_back(), 7);
        // refine
        ASSERT_EQ(interp.nrefine_front(), 2);
        ASSERT_EQ(interp.nrefine_back(), 3);
        // criterion
        ASSERT_EQ(interp.ncriterion_front(), 8);
        ASSERT_EQ(interp.ncriterion_back(), 7);
        // shift
        ASSERT_EQ(interp.shift_front(), 3);
        ASSERT_EQ(interp.shift_back(), 2);
    }
}

TEST_F(valid_Wavelet_Kernel, coarsen_detail) {
    for (int id = 0; id < 3; id++) {
        for (int id = 0; id < 3; id++) {
            coarse_start_[id] = 0;
            coarse_end_[id]   = m_hn;
            fine_start_[id]   = 0;
            fine_end_[id]     = m_n;
        }
        block_coarse_->Reset(m_gs, m_stride, coarse_start_, coarse_end_);
        block_fine_->Reset(m_gs, m_stride, fine_start_, fine_end_);

        InterpolatingWavelet interp;

        real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
        real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

        real_p ptr_tmp = (real_t*)m_calloc(m_stride * m_stride * m_stride * sizeof(real_t));

        // fill the fine block!
        for (int i2 = -m_gs; i2 < (m_n + m_gs); i2++) {
            for (int i1 = -m_gs; i1 < (m_n + m_gs); i1++) {
                for (int i0 = -m_gs; i0 < (m_n + m_gs); i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly(pos[id]);
                }
            }
        }

        // do the interpolation (-1 is refinement)
        //................................................
        // get the coarse version of life
        const lid_t m_hgs = m_gs / 2;
        SubBlock    tmp_block(m_hgs, m_hn + 2 * m_hgs, -m_hgs, m_hn + m_hgs);
        data_ptr    data_tmp = ptr_tmp + m_zeroidx(0, &tmp_block);
        for (lid_t i2 = tmp_block.start(2); i2 < tmp_block.end(2); i2++) {
            for (lid_t i1 = tmp_block.start(1); i1 < tmp_block.end(1); i1++) {
                for (lid_t i0 = tmp_block.start(0); i0 < tmp_block.end(0); i0++) {
                    data_tmp[m_midx(i0, i1, i2, 0, &tmp_block)] = data_fine[m_midx(i0 * 2, i1 * 2, i2 * 2, 0, block_fine_)];
                }
            }
        }
        real_t detail_max;
        interp.Details(block_fine_, data_fine, &tmp_block, data_tmp, &detail_max);
        ASSERT_NEAR(detail_max, 0.0, DOUBLE_TOL) << "during test: dir = " << id << " detail max = " << detail_max;

        m_free(ptr_tmp);

        // reset the indexes for the fine
        for (int id = 0; id < 3; id++) {
            fine_start_[id] = -m_gs;
            fine_end_[id]   = m_n + m_gs;
        }
        block_fine_->Reset(m_gs, m_stride, fine_start_, fine_end_);

        // do the coarsening
        lid_t shift[3] = {0};
        interp.Interpolate(1, shift, block_fine_, data_fine, block_coarse_, data_coarse);

        // get the moment, must be the same along every dimension
        const lid_t dir0 = id;
        const lid_t dir1 = (id + 1) % 3;
        const lid_t dir2 = (id + 2) % 3;

        real_t mom_coarse[M_WAVELET_NT];
        real_t mom_fine[M_WAVELET_NT];
        real_t f_coarse[M_WAVELET_NT][21];
        real_t f_fine[M_WAVELET_NT][41];

        // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)

        for (lid_t i2 = 0; i2 < 21; ++i2) {
            for (lid_t i1 = 0; i1 < 21; ++i1) {
                // fill the f array with the correct dimension
                real_t* datastick_fine;
                real_t* datastick_coarse;
                //------------------------
                if (dir0 == 0) {
                    datastick_fine   = data_fine + m_midx(0, 2 * i1, 2 * i2, 0, block_fine_);
                    datastick_coarse = data_coarse + m_midx(0, i1, i2, 0, block_fine_);
                    for (lid_t i0 = 0; i0 < 21; ++i0) {
                        const real_t x0 = (2 * i0) * hfine_;
                        const real_t x1 = (2 * i0 + 1) * hfine_;
                        //check the results
                        real_t val_coarse = datastick_coarse[m_midx(i0, 0, 0, 0, block_coarse_)];
                        real_t val_fine0  = datastick_fine[m_midx(2 * i0, 0, 0, 0, block_fine_)];
                        real_t val_fine1  = datastick_fine[m_midx(2 * i0 + 1, 0, 0, 0, block_fine_)];
                        // m_log("checking value: %f vs %f", val_fine0, val_coarse);
                        ASSERT_NEAR(poly(x0), val_coarse, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;

                        // fill the f array
                        for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                            f_coarse[im][i0]   = pow(x0, im) * val_coarse;
                            f_fine[im][i0]     = pow(x0, im) * val_fine0;
                            f_fine[im][i0 + 1] = pow(x1, im) * val_fine1;
                        }
                    }
                }
                //------------------------
                else if (dir0 == 1) {
                    datastick_fine   = data_fine + m_midx(2 * i1, 0, 2 * i2, 0, block_fine_);
                    datastick_coarse = data_coarse + m_midx(i1, 0, i2, 0, block_fine_);
                    for (lid_t i0 = 0; i0 < 21; ++i0) {
                        const real_t x0 = (2 * i0) * hfine_;
                        const real_t x1 = (2 * i0 + 1) * hfine_;
                        //check the results
                        real_t val_coarse = datastick_coarse[m_midx(0, i0, 0, 0, block_coarse_)];
                        real_t val_fine0  = datastick_fine[m_midx(0, 2 * i0, 0, 0, block_fine_)];
                        real_t val_fine1  = datastick_fine[m_midx(0, 2 * i0 + 1, 0, 0, block_fine_)];
                        // m_log("checking value: %f vs %f", val_fine0, val_coarse);
                        ASSERT_NEAR(poly(x0), val_coarse, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;

                        // fill the f array
                        for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                            f_coarse[im][i0]   = pow(x0, im) * val_coarse;
                            f_fine[im][i0]     = pow(x0, im) * val_fine0;
                            f_fine[im][i0 + 1] = pow(x1, im) * val_fine1;
                        }
                    }
                }
                //------------------------
                else if (dir0 == 2) {
                    datastick_fine   = data_fine + m_midx(2 * i1, 2 * i2, 0, 0, block_fine_);
                    datastick_coarse = data_coarse + m_midx(i1, i2, 0, 0, block_fine_);
                    for (lid_t i0 = 0; i0 < 21; ++i0) {
                        const real_t x0 = (2 * i0) * hfine_;
                        const real_t x1 = (2 * i0 + 1) * hfine_;
                        //check the results
                        real_t val_coarse = datastick_coarse[m_midx(0, 0, i0, 0, block_coarse_)];
                        real_t val_fine0  = datastick_fine[m_midx(0, 0, 2 * i0, 0, block_fine_)];
                        real_t val_fine1  = datastick_fine[m_midx(0, 0, 2 * i0 + 1, 0, block_fine_)];
                        // m_log("checking value: %f vs %f", val_fine0, val_coarse);
                        ASSERT_NEAR(poly(x0), val_coarse, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;

                        // fill the f array
                        for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                            f_coarse[im][i0]   = pow(x0, im) * val_coarse;
                            f_fine[im][i0]     = pow(x0, im) * val_fine0;
                            f_fine[im][i0 + 1] = pow(x1, im) * val_fine1;
                        }
                    }
                }

                //------------------------
                // check the moment
                for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                    real_t mom_coarse = Integral(21, f_coarse[im], hcoarse_);
                    real_t mom_fine   = Integral(41, f_fine[im], hcoarse_);
                    // m_log("moment = %f vs %f",mom_coarse,mom_fine);
                    ASSERT_NEAR(mom_fine, mom_coarse, DOUBLE_TOL) << "two moment do not match: fine = " << mom_fine << "and coarse = " << mom_coarse;
                }
            }
        }
    }
}

TEST_F(valid_Wavelet_Kernel, refine) {
    for (int id = 0; id < 3; id++) {
        for (int id = 0; id < 3; id++) {
            coarse_start_[id] = -m_gs;
            coarse_end_[id]   = m_hn + m_gs;
            fine_start_[id]   = 0;
            fine_end_[id]     = m_n;
        }
        block_coarse_->Reset(m_gs, m_stride, coarse_start_, coarse_end_);
        block_fine_->Reset(m_gs, m_stride, fine_start_, fine_end_);

        real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
        real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

        // fill the fine block!
        for (int i2 = -m_gs; i2 < (m_hn + m_gs); i2++) {
            for (int i1 = -m_gs; i1 < (m_hn + m_gs); i1++) {
                for (int i0 = -m_gs; i0 < (m_n + m_gs); i0++) {
                    real_t x      = i0 * hcoarse_;
                    real_t y      = i1 * hcoarse_;
                    real_t z      = i2 * hcoarse_;
                    real_t pos[3] = {x, y, z};

                    data_coarse[m_midx(i0, i1, i2, 0, block_fine_)] = poly(pos[id]);
                }
            }
        }

        // do the coarsening
        InterpolatingWavelet interp;
        lid_t                shift[3] = {0};
        interp.Interpolate(-1, shift, block_coarse_, data_coarse, block_fine_, data_fine);

        // get the moment, must be the same along every dimension
        const lid_t dir0 = id;
        const lid_t dir1 = (id + 1) % 3;
        const lid_t dir2 = (id + 2) % 3;

        real_t mom_coarse[M_WAVELET_NT];
        real_t mom_fine[M_WAVELET_NT];
        real_t f_coarse[M_WAVELET_NT][21];
        real_t f_fine[M_WAVELET_NT][41];

        // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)

        for (lid_t i2 = 0; i2 < 21; ++i2) {
            for (lid_t i1 = 0; i1 < 21; ++i1) {
                // fill the f array with the correct dimension
                real_t* datastick_fine;
                real_t* datastick_coarse;
                //------------------------
                if (dir0 == 0) {
                    datastick_fine   = data_fine + m_midx(0, 2 * i1, 2 * i2, 0, block_fine_);
                    datastick_coarse = data_coarse + m_midx(0, i1, i2, 0, block_fine_);
                    for (lid_t i0 = 0; i0 < 21; ++i0) {
                        const real_t x0 = (2 * i0) * hfine_;
                        const real_t x1 = (2 * i0 + 1) * hfine_;
                        //check the results
                        real_t val_coarse = datastick_coarse[m_midx(i0, 0, 0, 0, block_coarse_)];
                        real_t val_fine0  = datastick_fine[m_midx(2 * i0, 0, 0, 0, block_fine_)];
                        real_t val_fine1  = datastick_fine[m_midx(2 * i0 + 1, 0, 0, 0, block_fine_)];
                        // m_log("checking value: %f vs %f", val_fine0, val_coarse);
                        ASSERT_NEAR(poly(x0), val_fine0, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;
                        ASSERT_NEAR(poly(x1), val_fine1, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;

                        // fill the f array
                        for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                            f_coarse[im][i0]   = pow(x0, im) * val_coarse;
                            f_fine[im][i0]     = pow(x0, im) * val_fine0;
                            f_fine[im][i0 + 1] = pow(x1, im) * val_fine1;
                        }
                    }
                }
                //------------------------
                else if (dir0 == 1) {
                    datastick_fine   = data_fine + m_midx(2 * i1, 0, 2 * i2, 0, block_fine_);
                    datastick_coarse = data_coarse + m_midx(i1, 0, i2, 0, block_fine_);
                    for (lid_t i0 = 0; i0 < 21; ++i0) {
                        const real_t x0 = (2 * i0) * hfine_;
                        const real_t x1 = (2 * i0 + 1) * hfine_;
                        //check the results
                        real_t val_coarse = datastick_coarse[m_midx(0, i0, 0, 0, block_coarse_)];
                        real_t val_fine0  = datastick_fine[m_midx(0, 2 * i0, 0, 0, block_fine_)];
                        real_t val_fine1  = datastick_fine[m_midx(0, 2 * i0 + 1, 0, 0, block_fine_)];
                        // m_log("checking value: %f vs %f", val_fine0, val_coarse);
                        ASSERT_NEAR(poly(x0), val_fine0, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;
                        ASSERT_NEAR(poly(x1), val_fine1, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;

                        // fill the f array
                        for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                            f_coarse[im][i0]   = pow(x0, im) * val_coarse;
                            f_fine[im][i0]     = pow(x0, im) * val_fine0;
                            f_fine[im][i0 + 1] = pow(x1, im) * val_fine1;
                        }
                    }
                }
                //------------------------
                else if (dir0 == 2) {
                    datastick_fine   = data_fine + m_midx(2 * i1, 2 * i2, 0, 0, block_fine_);
                    datastick_coarse = data_coarse + m_midx(i1, i2, 0, 0, block_fine_);
                    for (lid_t i0 = 0; i0 < 21; ++i0) {
                        const real_t x0 = (2 * i0) * hfine_;
                        const real_t x1 = (2 * i0 + 1) * hfine_;
                        //check the results
                        real_t val_coarse = datastick_coarse[m_midx(0, 0, i0, 0, block_coarse_)];
                        real_t val_fine0  = datastick_fine[m_midx(0, 0, 2 * i0, 0, block_fine_)];
                        real_t val_fine1  = datastick_fine[m_midx(0, 0, 2 * i0 + 1, 0, block_fine_)];
                        // m_log("checking value: %f vs %f", val_fine0, val_coarse);
                        ASSERT_NEAR(poly(x0), val_fine0, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;
                        ASSERT_NEAR(poly(x1), val_fine1, DOUBLE_TOL) << "during test: dir = " << id << " @ " << i0 << " " << i1 << " " << i2;

                        // fill the f array
                        for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                            f_coarse[im][i0]   = pow(x0, im) * val_coarse;
                            f_fine[im][i0]     = pow(x0, im) * val_fine0;
                            f_fine[im][i0 + 1] = pow(x1, im) * val_fine1;
                        }
                    }
                }

                //------------------------
                // check the moment
                for (lda_t im = 0; im < M_WAVELET_NT; im++) {
                    real_t mom_coarse = Integral(21, f_coarse[im], hcoarse_);
                    real_t mom_fine   = Integral(41, f_fine[im], hcoarse_);
                    // m_log("moment = %f vs %f",mom_coarse,mom_fine);
                    ASSERT_NEAR(mom_fine, mom_coarse, DOUBLE_TOL) << "two moment do not match: fine = " << mom_fine << "and coarse = " << mom_coarse;
                }
            }
        }
    }
}
