
#include "wavelet.hpp"

#include "error.hpp"
#include "gtest/gtest.h"
#include "murphy.hpp"
#include "subblock.hpp"

#define DOUBLE_TOL 1e-13

class valid_Wavelet_Ghost : public ::testing::Test {
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
        data_fine_   = (real_t*)m_calloc(20 * 20 * 20 * sizeof(real_t));
        data_coarse_ = (real_t*)m_calloc(20 * 20 * 20 * sizeof(real_t));

        for (int id = 0; id < 3; id++) {
            coarse_start_[id] = 0;
            coarse_end_[id]   = M_N;
            fine_start_[id]   = -M_GS;
            fine_end_[id]     = 0;
        }
        block_coarse_ = new SubBlock(M_GS, M_STRIDE, coarse_start_, coarse_end_);
        block_fine_   = new SubBlock(M_GS, M_STRIDE, fine_start_, fine_end_);

        hcoarse_ = 1.0 / (M_N);
        hfine_   = 0.5 / (M_N);
    };
    void TearDown() override {
        m_free(data_fine_);
        m_free(data_coarse_);

        delete (block_coarse_);
        delete (block_fine_);
    };
};

//==============================================================================================================================
static double poly_1(const double x) {
    return x;// M_PI_2 * x + M_SQRT2;
}
// test the wavelet + moments
TEST_F(valid_Wavelet_Ghost, ghost_order_2_2) {
    real_t* data_coarse = data_coarse_ + m_zeroidx(0, block_coarse_);
    real_t* data_fine   = data_fine_ + m_zeroidx(0, block_fine_);
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = 0;
        coarse_end_[id]   = M_N;
        fine_start_[id]   = 0;
        fine_end_[id]     = M_N;
    }
    fine_start_[0] = -M_GS;
    fine_end_[0]   = 0;
    fine_end_[1]   = M_N - 2;
    fine_end_[2]   = M_N - 2;

    block_coarse_->Reset(M_GS, M_STRIDE, coarse_start_, coarse_end_);
    block_fine_->Reset(M_GS, M_STRIDE, fine_start_, fine_end_);

    for (int id = 0; id < 3; id++) {
        real_t mom_coarse[4] = {0.0, 0.0, 0.0, 0.0};
        real_t mom_fine[4]   = {0.0, 0.0, 0.0, 0.0};
        real_t vol_fine      = (hfine_ * hfine_ * hfine_);
        real_t vol_coarse    = (hcoarse_ * hcoarse_ * hcoarse_);

        real_t shift_fine[3]   = {1.0, 0.0, 0.0};
        real_t shift_coarse[3] = {0.0, 0.0, 0.0};

        // fill the fine block
        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    real_t x      = shift_fine[0] + i0 * hfine_;
                    real_t y      = shift_fine[1] + i1 * hfine_;
                    real_t z      = shift_fine[2] + i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    data_fine[m_midx(i0, i1, i2, 0, block_fine_)] = poly_1(pos[id]);
                }
            }
        }


        // fill the coarse block
        for (int i2 = 0; i2 < M_N; i2++) {
            for (int i1 = 0; i1 < M_N; i1++) {
                for (int i0 = 0; i0 < M_N; i0++) {
                    real_t x      = shift_coarse[0] + i0 * hcoarse_;
                    real_t y      = shift_coarse[1] + i1 * hcoarse_;
                    real_t z      = shift_coarse[2] + i2 * hcoarse_;
                    real_t pos[3] = {x, y, z};

                    data_coarse[m_midx(i0, i1, i2, 0, block_coarse_)] = poly_1(pos[id]);
                }
            }
        }
        // do the refinement
        Wavelet<2, 2>* interp    = new Wavelet<2, 2>();
        lid_t          shift[3]  = {M_N, 0, 0};
        sid_t          normal[3] = {-1, 0, 0};
        interp->Interpolate(-1, shift, block_coarse_, data_coarse, block_fine_, data_fine, normal);

        //check the result
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = shift_fine[0] + i0 * hfine_;
                    real_t y      = shift_fine[1] + i1 * hfine_;
                    real_t z      = shift_fine[2] + i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    real_t val = data_fine[m_midx(i0, i1, i2, 0, block_fine_)];
                    ASSERT_NEAR(val, poly_1(pos[id]), DOUBLE_TOL) << "in" << i0 << " " << i1 << " " << i2 << " while testing direction " << id;
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
        delete (interp);
    }
}
