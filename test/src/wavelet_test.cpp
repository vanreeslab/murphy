
#include "wavelet.hpp"

#include "error.hpp"
#include "gtest/gtest.h"
#include "murphy.hpp"
#include "subblock.hpp"

#define DOUBLE_TOL 1e-13

class valid_Wavelet : public ::testing::Test {
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
        data_fine_   = (real_t*)m_calloc(M_STRIDE * M_STRIDE * M_STRIDE * sizeof(real_t));
        data_coarse_ = (real_t*)m_calloc(M_STRIDE * M_STRIDE * M_STRIDE * sizeof(real_t));

        for (int id = 0; id < 3; id++) {
            coarse_start_[id] = -M_GS;
            coarse_end_[id]   = M_HN + M_GS;
            fine_start_[id]   = 0;
            fine_end_[id]     = M_N;
        }
        block_coarse_ = new SubBlock(M_GS, M_HN + 2 * M_GS, coarse_start_, coarse_end_);
        block_fine_   = new SubBlock(M_GS, M_STRIDE, fine_start_, fine_end_);

        hcoarse_ = 1.0 / (M_HN);
        hfine_   = 1.0 / (M_N);
    };
    void TearDown() override {
        m_free(data_fine_);
        m_free(data_coarse_);

        delete (block_coarse_);
        delete (block_fine_);
    };
};

//==============================================================================================================================
#if (WAVELET_N == 2 && WAVELET_NT == 2)
static double poly_1(const double x) {
    return M_PI_2 * x + M_SQRT2;
}
// test the wavelet + moments
TEST_F(valid_Wavelet, coarsen_order_2_2) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = 0;
        coarse_end_[id]   = 8;
        fine_start_[id]   = -2;
        fine_end_[id]     = 18;
    }
    block_coarse_->Reset(2, 20, coarse_start_, coarse_end_);
    block_fine_->Reset(2, 20, fine_start_, fine_end_);

    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);
        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_1(pos[id]);
                }
            }
        }
        // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)
        for (int i2 = 2 * coarse_start_[2]; i2 < 2 * coarse_end_[2] - 3; i2++) {
            for (int i1 = 2 * coarse_start_[1]; i1 < 2 * coarse_end_[1] - 3; i1++) {
                for (int i0 = 2 * coarse_start_[0]; i0 < 2 * coarse_end_[0] - 3; i0++) {
                    real_t x = i0 * hfine_;
                    real_t y = i1 * hfine_;
                    real_t z = i2 * hfine_;
                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (2 * coarse_start_[0]) || i0 == (2 * coarse_end_[0] - 4)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (2 * coarse_start_[1]) || i1 == (2 * coarse_end_[1] - 4)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (2 * coarse_start_[2]) || i2 == (2 * coarse_end_[2] - 4)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    mom_coarse[0] += val * vol_fine * corr0 * corr1 * corr2;
                    mom_coarse[1] += x * val * vol_fine * corr0 * corr1 * corr2;
                    mom_coarse[2] += y * val * vol_fine * corr0 * corr1 * corr2;
                    mom_coarse[3] += z * val * vol_fine * corr0 * corr1 * corr2;
                }
            }
        }

        // do the coarsening
        Wavelet interp;
        lid_t   shift[3] = {0};
        interp.Interpolate(1, shift, block_fine_, data_fine, block_coarse_, data_coarse);

        //check the result
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x      = i0 * hcoarse_;
                    real_t y      = i1 * hcoarse_;
                    real_t z      = i2 * hcoarse_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];
                    ASSERT_NEAR(val, poly_1(pos[id]), DOUBLE_TOL) << "while testing direction " << id;
                }
            }
        }
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2] - 1; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1] - 1; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0] - 1; i0++) {
                    real_t x = i0 * hcoarse_;
                    real_t y = i1 * hcoarse_;
                    real_t z = i2 * hcoarse_;

                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];

                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (coarse_start_[0]) || i0 == (coarse_end_[0] - 2)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (coarse_start_[1]) || i1 == (coarse_end_[1] - 2)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (coarse_start_[2]) || i2 == (coarse_end_[2] - 2)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    mom_fine[0] += val * vol_coarse * corr0 * corr1 * corr2;
                    mom_fine[1] += x * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_fine[2] += y * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_fine[3] += z * val * vol_coarse * corr0 * corr1 * corr2;
                }
                // m_verb("-------------\n");
            }
        }
        printf("moments 0: %e vs %e\n", mom_coarse[0], mom_fine[0]);
        printf("moments 1 in x: %e vs %e\n", mom_coarse[1], mom_fine[1]);
        printf("moments 1 in y: %e vs %e\n", mom_coarse[2], mom_fine[2]);
        printf("moments 1 in z: %e vs %e\n", mom_coarse[3], mom_fine[3]);
        ASSERT_NEAR(mom_coarse[0], mom_fine[0], DOUBLE_TOL);
        ASSERT_NEAR(mom_coarse[1], mom_fine[1], DOUBLE_TOL);
        ASSERT_NEAR(mom_coarse[2], mom_fine[2], DOUBLE_TOL);
        ASSERT_NEAR(mom_coarse[3], mom_fine[3], DOUBLE_TOL);
    }
}

// test the wavelets + moments
TEST_F(valid_Wavelet, refine_order_2_2) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = -2;
        coarse_end_[id]   = 10;
        fine_start_[id]   = 0;
        fine_end_[id]     = 16;
    }

    block_coarse_->Reset(2, 20, coarse_start_, coarse_end_);
    block_fine_->Reset(2, 20, fine_start_, fine_end_);

    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);
        // fill the function
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x                                          = i0 * hcoarse_;
                    real_t y                                          = i1 * hcoarse_;
                    real_t z                                          = i2 * hcoarse_;
                    real_t pos[3]                                     = {x, y, z};
                    data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)] = poly_1(pos[id]);
                }
            }
        }
        // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)
        for (int i2 = fine_start_[2] / 2; i2 < fine_end_[2] / 2 - 1; i2++) {
            for (int i1 = fine_start_[1] / 2; i1 < fine_end_[1] / 2 - 1; i1++) {
                for (int i0 = fine_start_[0] / 2; i0 < fine_end_[0] / 2 - 1; i0++) {
                    real_t x = i0 * hcoarse_;
                    real_t y = i1 * hcoarse_;
                    real_t z = i2 * hcoarse_;
                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (fine_start_[0] / 2) || i0 == (fine_end_[0] / 2 - 2)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (fine_start_[1] / 2) || i1 == (fine_end_[1] / 2 - 2)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (fine_start_[2] / 2) || i2 == (fine_end_[2] / 2 - 2)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];
                    // real_t val =  1.0;//data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    mom_coarse[0] += val * vol_coarse * corr0 * corr1 * corr2;
                    mom_coarse[1] += x * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_coarse[2] += y * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_coarse[3] += z * val * vol_coarse * corr0 * corr1 * corr2;
                }
            }
        }

        // do the coarsening
        Wavelet interp;
        lid_t   shift[3] = {0};
        interp.Interpolate(-1, shift, block_coarse_, data_coarse, block_fine_, data_fine);

        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    ASSERT_NEAR(val, poly_1(pos[id]), DOUBLE_TOL) << "while testing direction " << id << "@ " << i0 << " " << i1 << " " << i2;
                }
            }
        }
        for (int i2 = fine_start_[2]; i2 < fine_end_[2] - 3; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1] - 3; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0] - 3; i0++) {
                    real_t x = i0 * hfine_;
                    real_t y = i1 * hfine_;
                    real_t z = i2 * hfine_;

                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    // real_t val = 1.0;//data_fine[m_midx(i0, i1, i2, 0, block_fine_)];

                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (fine_start_[0]) || i0 == (fine_end_[0] - 4)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (fine_start_[1]) || i1 == (fine_end_[1] - 4)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (fine_start_[2]) || i2 == (fine_end_[2] - 4)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    mom_fine[0] += val * vol_fine * corr0 * corr1 * corr2;
                    mom_fine[1] += x * val * vol_fine * corr0 * corr1 * corr2;
                    mom_fine[2] += y * val * vol_fine * corr0 * corr1 * corr2;
                    mom_fine[3] += z * val * vol_fine * corr0 * corr1 * corr2;
                }
            }
        }
        printf("moments 0: %e vs %e\n", mom_coarse[0], mom_fine[0]);
        printf("moments 1 in x: %e vs %e\n", mom_coarse[1], mom_fine[1]);
        printf("moments 1 in y: %e vs %e\n", mom_coarse[2], mom_fine[2]);
        printf("moments 1 in z: %e vs %e\n", mom_coarse[3], mom_fine[3]);
        ASSERT_NEAR(mom_coarse[0], mom_fine[0], DOUBLE_TOL);
        ASSERT_NEAR(mom_coarse[1], mom_fine[1], DOUBLE_TOL);
        ASSERT_NEAR(mom_coarse[2], mom_fine[2], DOUBLE_TOL);
        ASSERT_NEAR(mom_coarse[3], mom_fine[3], DOUBLE_TOL);
    }
}

// test the details
TEST_F(valid_Wavelet, detail_order_2_2) {
    real_t* data_fine = data_fine_ + m_zeroidx(0, block_fine_);
    // fill the source
    for (int i2 = fine_start_[2] - block_fine_->gs(); i2 < fine_end_[2] + block_fine_->gs(); i2++) {
        for (int i1 = fine_start_[1] - block_fine_->gs(); i1 < fine_end_[1] + block_fine_->gs(); i1++) {
            for (int i0 = fine_start_[0] - block_fine_->gs(); i0 < fine_end_[0] + block_fine_->gs(); i0++) {
                real_t x = i0 * hfine_;
                real_t y = i1 * hfine_;
                real_t z = i2 * hfine_;

                data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_1(x) + poly_1(y) + poly_1(z);
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet interp;
    real_t  detail_max;
    interp.Details(block_fine_, data_fine_, &detail_max);
    ASSERT_NEAR(detail_max, 0.0, DOUBLE_TOL);
}
//==============================================================================================================================
#elif (WAVELET_N == 4 && WAVELET_NT == 0)

static double poly_3(const double x) {
    return M_PI * pow(x, 3) + M_PI_2 * pow(x, 2) + M_SQRT2 * x + M_PI_4;
}
// test the wavelet + moments
TEST_F(valid_Wavelet, coarsen_order_4_0) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = 0;
        coarse_end_[id]   = 8;
        fine_start_[id]   = -2;
        fine_end_[id]     = 18;
    }
    block_coarse_->Reset(2, 20, coarse_start_, coarse_end_);
    block_fine_->Reset(2, 20, fine_start_, fine_end_);

    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);
        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_3(pos[id]);
                }
            }
        }
        // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)
        for (int i2 = 2 * coarse_start_[2]; i2 < 2 * coarse_end_[2] - 3; i2++) {
            for (int i1 = 2 * coarse_start_[1]; i1 < 2 * coarse_end_[1] - 3; i1++) {
                for (int i0 = 2 * coarse_start_[0]; i0 < 2 * coarse_end_[0] - 3; i0++) {
                    real_t x = i0 * hfine_;
                    real_t y = i1 * hfine_;
                    real_t z = i2 * hfine_;
                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (2 * coarse_start_[0]) || i0 == (2 * coarse_end_[0] - 4)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (2 * coarse_start_[1]) || i1 == (2 * coarse_end_[1] - 4)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (2 * coarse_start_[2]) || i2 == (2 * coarse_end_[2] - 4)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    mom_coarse[0] += val * vol_fine * corr0 * corr1 * corr2;
                    mom_coarse[1] += x * val * vol_fine * corr0 * corr1 * corr2;
                    mom_coarse[2] += y * val * vol_fine * corr0 * corr1 * corr2;
                    mom_coarse[3] += z * val * vol_fine * corr0 * corr1 * corr2;
                }
            }
        }

        // do the coarsening
        Wavelet interp;
        lid_t   shift[3] = {0};
        interp.Interpolate(1, shift, block_fine_, data_fine, block_coarse_, data_coarse);

        //check the result
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x      = i0 * hcoarse_;
                    real_t y      = i1 * hcoarse_;
                    real_t z      = i2 * hcoarse_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];
                    ASSERT_NEAR(val, poly_3(pos[id]), DOUBLE_TOL) << "while testing direction " << id;
                }
            }
        }
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2] - 1; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1] - 1; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0] - 1; i0++) {
                    real_t x = i0 * hcoarse_;
                    real_t y = i1 * hcoarse_;
                    real_t z = i2 * hcoarse_;

                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];

                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (coarse_start_[0]) || i0 == (coarse_end_[0] - 2)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (coarse_start_[1]) || i1 == (coarse_end_[1] - 2)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (coarse_start_[2]) || i2 == (coarse_end_[2] - 2)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    mom_fine[0] += val * vol_coarse * corr0 * corr1 * corr2;
                    mom_fine[1] += x * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_fine[2] += y * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_fine[3] += z * val * vol_coarse * corr0 * corr1 * corr2;
                }
                // m_verb("-------------\n");
            }
        }
        printf("moments 0: %e vs %e\n", mom_coarse[0], mom_fine[0]);
        printf("moments 1 in x: %e vs %e\n", mom_coarse[1], mom_fine[1]);
        printf("moments 1 in y: %e vs %e\n", mom_coarse[2], mom_fine[2]);
        printf("moments 1 in z: %e vs %e\n", mom_coarse[3], mom_fine[3]);
        // ASSERT_NEAR(mom_coarse[0], mom_fine[0], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[1], mom_fine[1], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[2], mom_fine[2], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[3], mom_fine[3], DOUBLE_TOL);
    }
}
// test the wavelets + moments
TEST_F(valid_Wavelet, refine_order_4_0) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = -2;
        coarse_end_[id]   = 10;
        fine_start_[id]   = 0;
        fine_end_[id]     = 16;
    }

    block_coarse_->Reset(2, 20, coarse_start_, coarse_end_);
    block_fine_->Reset(2, 20, fine_start_, fine_end_);

    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);
        // fill the function
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x                                          = i0 * hcoarse_;
                    real_t y                                          = i1 * hcoarse_;
                    real_t z                                          = i2 * hcoarse_;
                    real_t pos[3]                                     = {x, y, z};
                    data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)] = poly_3(pos[id]);
                }
            }
        }
        // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)
        for (int i2 = fine_start_[2] / 2; i2 < fine_end_[2] / 2 - 1; i2++) {
            for (int i1 = fine_start_[1] / 2; i1 < fine_end_[1] / 2 - 1; i1++) {
                for (int i0 = fine_start_[0] / 2; i0 < fine_end_[0] / 2 - 1; i0++) {
                    real_t x = i0 * hcoarse_;
                    real_t y = i1 * hcoarse_;
                    real_t z = i2 * hcoarse_;
                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (fine_start_[0] / 2) || i0 == (fine_end_[0] / 2 - 2)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (fine_start_[1] / 2) || i1 == (fine_end_[1] / 2 - 2)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (fine_start_[2] / 2) || i2 == (fine_end_[2] / 2 - 2)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];
                    // real_t val =  1.0;//data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    mom_coarse[0] += val * vol_coarse * corr0 * corr1 * corr2;
                    mom_coarse[1] += x * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_coarse[2] += y * val * vol_coarse * corr0 * corr1 * corr2;
                    mom_coarse[3] += z * val * vol_coarse * corr0 * corr1 * corr2;
                }
            }
        }

        // do the coarsening
        Wavelet interp;
        lid_t   shift[3] = {0};
        interp.Interpolate(-1, shift, block_coarse_, data_coarse, block_fine_, data_fine);

        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    ASSERT_NEAR(val, poly_3(pos[id]), DOUBLE_TOL) << "while testing direction " << id;
                }
            }
        }
        for (int i2 = fine_start_[2]; i2 < fine_end_[2] - 3; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1] - 3; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0] - 3; i0++) {
                    real_t x = i0 * hfine_;
                    real_t y = i1 * hfine_;
                    real_t z = i2 * hfine_;

                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    // real_t val = 1.0;//data_fine[m_midx(i0, i1, i2, 0, block_fine_)];

                    // trapezoidal moment computation
                    real_t corr0 = (i0 == (fine_start_[0]) || i0 == (fine_end_[0] - 4)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr1 = (i1 == (fine_start_[1]) || i1 == (fine_end_[1] - 4)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    real_t corr2 = (i2 == (fine_start_[2]) || i2 == (fine_end_[2] - 4)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
                    // get the moments
                    mom_fine[0] += val * vol_fine * corr0 * corr1 * corr2;
                    mom_fine[1] += x * val * vol_fine * corr0 * corr1 * corr2;
                    mom_fine[2] += y * val * vol_fine * corr0 * corr1 * corr2;
                    mom_fine[3] += z * val * vol_fine * corr0 * corr1 * corr2;
                }
            }
        }

        printf("moments 0: %e vs %e\n", mom_coarse[0], mom_fine[0]);
        printf("moments 1 in x: %e vs %e\n", mom_coarse[1], mom_fine[1]);
        printf("moments 1 in y: %e vs %e\n", mom_coarse[2], mom_fine[2]);
        printf("moments 1 in z: %e vs %e\n", mom_coarse[3], mom_fine[3]);
        // ASSERT_NEAR(mom_coarse[0], mom_fine[0], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[1], mom_fine[1], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[2], mom_fine[2], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[3], mom_fine[3], DOUBLE_TOL);
    }
}

// test the details
TEST_F(valid_Wavelet, detail_order_4_0) {
    for (int id = 0; id < 3; id++) {
        fine_start_[id] = 0;
        fine_end_[id]   = 12;
    }
    block_fine_->Reset(4, 20, fine_start_, fine_end_);
    real_t* data_fine = data_fine_ + m_zeroidx(0, block_fine_);
    // fill the source
    m_verb("filling between %d %d %d", fine_start_[0] - block_fine_->gs(), fine_start_[1] - block_fine_->gs(), fine_start_[2] - block_fine_->gs());
    for (int i2 = fine_start_[2] - block_fine_->gs(); i2 < fine_end_[2] + block_fine_->gs(); i2++) {
        for (int i1 = fine_start_[1] - block_fine_->gs(); i1 < fine_end_[1] + block_fine_->gs(); i1++) {
            for (int i0 = fine_start_[0] - block_fine_->gs(); i0 < fine_end_[0] + block_fine_->gs(); i0++) {
                real_t x = i0 * hfine_;
                real_t y = i1 * hfine_;
                real_t z = i2 * hfine_;

                data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_3(x) + poly_3(y) + poly_3(z);
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet interp;
    real_t  detail_max;
    interp.Details(block_fine_, data_fine, &detail_max);

    ASSERT_NEAR(detail_max, 0.0, DOUBLE_TOL);  // d_x = 0.0
}

//==============================================================================================================================
// test the wavelet + moments
#elif (WAVELET_N == 4 && WAVELET_NT == 2)


static double poly_3(const double x) {
    return M_PI * pow(x, 3) + M_PI_2 * pow(x, 2) + M_SQRT2 * x + M_PI_4;
}

TEST_F(valid_Wavelet, coarsen_order_4_2) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = 0;
        coarse_end_[id]   = 6;
        fine_start_[id]   = -4;
        fine_end_[id]     = 16;
    }
    block_coarse_->Reset(4, 20, coarse_start_, coarse_end_);
    block_fine_->Reset(4, 20, fine_start_, fine_end_);

    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);
        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_3(pos[id]);
                }
            }
        }
        // // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)
        // for (int i2 = 2 * coarse_start_[2]; i2 < 2 * coarse_end_[2] - 3; i2++) {
        //     for (int i1 = 2 * coarse_start_[1]; i1 < 2 * coarse_end_[1] - 3; i1++) {
        //         for (int i0 = 2 * coarse_start_[0]; i0 < 2 * coarse_end_[0] - 3; i0++) {
        //             real_t x = i0 * hfine_;
        //             real_t y = i1 * hfine_;
        //             real_t z = i2 * hfine_;
        //             // trapezoidal moment computation
        //             real_t corr0 = (i0 == (2 * coarse_start_[0]) || i0 == (2 * coarse_end_[0] - 4)) ? (1.0 / 3.0) : (i0%2) == 0 ? (2.0 / 3.0) : 4.0/3.0;
        //             real_t corr1 = (i1 == (2 * coarse_start_[1]) || i1 == (2 * coarse_end_[1] - 4)) ? (1.0 / 3.0) : (i1%2) == 0 ? (2.0 / 3.0) : 4.0/3.0;
        //             real_t corr2 = (i2 == (2 * coarse_start_[2]) || i2 == (2 * coarse_end_[2] - 4)) ? (1.0 / 3.0) : (i2%2) == 0 ? (2.0 / 3.0) : 4.0/3.0;
        //             // get the moments
        //             real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
        //             mom_coarse[0] += val * vol_fine * corr0 * corr1 * corr2;
        //             mom_coarse[1] += x * val * vol_fine * corr0 * corr1 * corr2;
        //             mom_coarse[2] += y * val * vol_fine * corr0 * corr1 * corr2;
        //             mom_coarse[3] += z * val * vol_fine * corr0 * corr1 * corr2;
        //         }
        //     }
        // }

        // do the coarsening
        Wavelet interp;
        lid_t   shift[3] = {0};
        interp.Interpolate(1, shift, block_fine_, data_fine, block_coarse_, data_coarse);

        //check the result
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x      = i0 * hcoarse_;
                    real_t y      = i1 * hcoarse_;
                    real_t z      = i2 * hcoarse_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];
                    ASSERT_NEAR(val, poly_3(pos[id]), DOUBLE_TOL) << "while testing direction " << id;
                }
            }
        }
        // for (int i2 = coarse_start_[2]; i2 < coarse_end_[2] - 1; i2++) {
        //     for (int i1 = coarse_start_[1]; i1 < coarse_end_[1] - 1; i1++) {
        //         for (int i0 = coarse_start_[0]; i0 < coarse_end_[0] - 1; i0++) {
        //             real_t x = i0 * hcoarse_;
        //             real_t y = i1 * hcoarse_;
        //             real_t z = i2 * hcoarse_;

        //             real_t val = data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)];

        //             // trapezoidal moment computation
        //             real_t corr0 = (i0 == (coarse_start_[0]) || i0 == (coarse_end_[0] - 2)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
        //             real_t corr1 = (i1 == (coarse_start_[1]) || i1 == (coarse_end_[1] - 2)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
        //             real_t corr2 = (i2 == (coarse_start_[2]) || i2 == (coarse_end_[2] - 2)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
        //             // get the moments
        //             mom_fine[0] += val * vol_coarse * corr0 * corr1 * corr2;
        //             mom_fine[1] += x * val * vol_coarse * corr0 * corr1 * corr2;
        //             mom_fine[2] += y * val * vol_coarse * corr0 * corr1 * corr2;
        //             mom_fine[3] += z * val * vol_coarse * corr0 * corr1 * corr2;
        //         }
        //         // m_verb("-------------\n");
        //     }
        // }
        // printf("moments 0: %e vs %e\n", mom_coarse[0], mom_fine[0]);
        // printf("moments 1 in x: %e vs %e\n", mom_coarse[1], mom_fine[1]);
        // printf("moments 1 in y: %e vs %e\n", mom_coarse[2], mom_fine[2]);
        // printf("moments 1 in z: %e vs %e\n", mom_coarse[3], mom_fine[3]);
        // ASSERT_NEAR(mom_coarse[0], mom_fine[0], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[1], mom_fine[1], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[2], mom_fine[2], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[3], mom_fine[3], DOUBLE_TOL);
    }
}

// test the wavelets + moments
TEST_F(valid_Wavelet, refine_order_4_2) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = -4;
        coarse_end_[id]   = 6 + 2;
        fine_start_[id]   = 0;
        fine_end_[id]     = 12;
    }

    block_coarse_->Reset(4, 20, coarse_start_, coarse_end_);
    block_fine_->Reset(4, 20, fine_start_, fine_end_);

    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);
        // fill the function
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x                                          = i0 * hcoarse_;
                    real_t y                                          = i1 * hcoarse_;
                    real_t z                                          = i2 * hcoarse_;
                    real_t pos[3]                                     = {x, y, z};
                    data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)] = poly_3(pos[id]);
                }
            }
        }
        // // we need to have an odd number of coarse points -> coarse -1 points (coarse = 8)
        // for (int i2 = fine_start_[2]/2; i2 < fine_end_[2]/2 - 1; i2++) {
        //     for (int i1 = fine_start_[1]/2; i1 < fine_end_[1]/2 - 1; i1++) {
        //         for (int i0 = fine_start_[0]/2; i0 < fine_end_[0]/2 - 1; i0++) {
        //             real_t x      = i0 * hcoarse_;
        //             real_t y      = i1 * hcoarse_;
        //             real_t z      = i2 * hcoarse_;
        //             // trapezoidal moment computation
        //             real_t corr0 = (i0 == (fine_start_[0]/2) || i0 == (fine_end_[0]/2 - 2)) ? (1.0 / 3.0) : (i0%2) == 0 ? (2.0 / 3.0) : 4.0/3.0;
        //             real_t corr1 = (i1 == (fine_start_[1]/2) || i1 == (fine_end_[1]/2 - 2)) ? (1.0 / 3.0) : (i1%2) == 0 ? (2.0 / 3.0) : 4.0/3.0;
        //             real_t corr2 = (i2 == (fine_start_[2]/2) || i2 == (fine_end_[2]/2 - 2)) ? (1.0 / 3.0) : (i2%2) == 0 ? (2.0 / 3.0) : 4.0/3.0;
        //             // get the moments
        //             real_t val =  data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)] ;
        //             // real_t val =  1.0;//data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
        //             mom_coarse[0] += val * vol_coarse * corr0 * corr1 * corr2;
        //             mom_coarse[1] += x * val * vol_coarse * corr0 * corr1 * corr2;
        //             mom_coarse[2] += y * val * vol_coarse * corr0 * corr1 * corr2;
        //             mom_coarse[3] += z * val * vol_coarse * corr0 * corr1 * corr2;
        //         }
        //     }
        // }

        // do the coarsening
        Wavelet interp;
        lid_t   shift[3] = {0};
        interp.Interpolate(-1, shift, block_coarse_, data_coarse, block_fine_, data_fine);

        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    ASSERT_NEAR(val, poly_3(pos[id]), DOUBLE_TOL) << "while testing direction " << id;
                }
            }
        }
        // for (int i2 = fine_start_[2]; i2 < fine_end_[2] - 3; i2++) {
        //     for (int i1 = fine_start_[1]; i1 < fine_end_[1] - 3; i1++) {
        //         for (int i0 = fine_start_[0]; i0 < fine_end_[0] - 3; i0++) {
        //             real_t x = i0 * hfine_;
        //             real_t y = i1 * hfine_;
        //             real_t z = i2 * hfine_;

        //             real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
        //             // real_t val = 1.0;//data_fine[m_midx(i0, i1, i2, 0, block_fine_)];

        //             // trapezoidal moment computation
        //             real_t corr0 = (i0 == (fine_start_[0]) || i0 == (fine_end_[0] - 4)) ? (1.0 / 3.0) : (i0 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
        //             real_t corr1 = (i1 == (fine_start_[1]) || i1 == (fine_end_[1] - 4)) ? (1.0 / 3.0) : (i1 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
        //             real_t corr2 = (i2 == (fine_start_[2]) || i2 == (fine_end_[2] - 4)) ? (1.0 / 3.0) : (i2 % 2) == 0 ? (2.0 / 3.0) : 4.0 / 3.0;
        //             // get the moments
        //             mom_fine[0] += val * vol_fine * corr0 * corr1 * corr2;
        //             mom_fine[1] += x * val * vol_fine * corr0 * corr1 * corr2;
        //             mom_fine[2] += y * val * vol_fine * corr0 * corr1 * corr2;
        //             mom_fine[3] += z * val * vol_fine * corr0 * corr1 * corr2;
        //         }
        //     }
        // }
        // printf("moments 0: %e vs %e\n", mom_coarse[0], mom_fine[0]);
        // printf("moments 1 in x: %e vs %e\n", mom_coarse[1], mom_fine[1]);
        // printf("moments 1 in y: %e vs %e\n", mom_coarse[2], mom_fine[2]);
        // printf("moments 1 in z: %e vs %e\n", mom_coarse[3], mom_fine[3]);
        // ASSERT_NEAR(mom_coarse[0], mom_fine[0], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[1], mom_fine[1], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[2], mom_fine[2], DOUBLE_TOL);
        // ASSERT_NEAR(mom_coarse[3], mom_fine[3], DOUBLE_TOL);
    }
}

// test the details
TEST_F(valid_Wavelet, detail_order_4_2) {
    for (int id = 0; id < 3; id++) {
        fine_start_[id] = 0;
        fine_end_[id]   = 12;
    }
    block_fine_->Reset(4, 20, fine_start_, fine_end_);
    real_t* data_fine = data_fine_ + m_zeroidx(0, block_fine_);
    // fill the source
    m_verb("filling between %d %d %d", fine_start_[0] - block_fine_->gs(), fine_start_[1] - block_fine_->gs(), fine_start_[2] - block_fine_->gs());
    for (int i2 = fine_start_[2] - block_fine_->gs(); i2 < fine_end_[2] + block_fine_->gs(); i2++) {
        for (int i1 = fine_start_[1] - block_fine_->gs(); i1 < fine_end_[1] + block_fine_->gs(); i1++) {
            for (int i0 = fine_start_[0] - block_fine_->gs(); i0 < fine_end_[0] + block_fine_->gs(); i0++) {
                real_t x = i0 * hfine_;
                real_t y = i1 * hfine_;
                real_t z = i2 * hfine_;

                data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_3(x) + poly_3(y) + poly_3(z);
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet interp;
    real_t  detail_max;
    interp.Details(block_fine_, data_fine, &detail_max);

    ASSERT_NEAR(detail_max, 0.0, DOUBLE_TOL);
}
#endif