
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

    void
    SetUp() override {
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
        data_coarse_  = data_coarse_ + m_zeroidx(0, block_coarse_);
        data_fine_    = data_fine_ + m_zeroidx(0, block_fine_);

        hcoarse_ = 1.0 / (M_HN);
        hfine_   = 1.0 / (M_N);
    };
    void TearDown() override {
        m_free(data_fine_ - m_zeroidx(0, block_fine_));
        m_free(data_coarse_ - m_zeroidx(0, block_coarse_));

        delete (block_coarse_);
        delete (block_fine_);
    };
};

//==============================================================================================================================
TEST_F(valid_Wavelet, refine_order_3_linear) {
    // fill the source
    for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
        for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
            for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hcoarse_;
                real_t y = ((real_t)i1 + 0.5) * hcoarse_;
                real_t z = ((real_t)i2 + 0.5) * hcoarse_;

                data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = x + y + z;
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet<3>* interp   = new Wavelet<3>();
    lid_t       shift[3] = {0};
    interp->Interpolate(-1, shift, block_coarse_, data_coarse_, block_fine_, data_fine_);

    // check the result

    for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
        for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
            for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hfine_;
                real_t y = ((real_t)i1 + 0.5) * hfine_;
                real_t z = ((real_t)i2 + 0.5) * hfine_;
                EXPECT_NEAR(data_fine_[m_midx(i0, i1, i2, 0, block_fine_)], x + y + z, DOUBLE_TOL) << "fails for " << i0 << " " << i1 << " " << i2 << " hgrid = " << hfine_;
            }
        }
    }
}
//==============================================================================================================================
// TEST_F(valid_Wavelet, order_3_quad) {
//     // fill the source
//     for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
//         for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
//             for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
//                 real_t x = ((real_t)i0 + 0.5) * hcoarse_;
//                 real_t y = ((real_t)i1 + 0.5) * hcoarse_;
//                 real_t z = ((real_t)i2 + 0.5) * hcoarse_;

//                 data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = x * x + y * y + z * z;
//             }
//         }
//     }

//     // do the interpolation (-1 is refinement)
//     Wavelet<3>* interp   = new Wavelet<3>();
//     lid_t       shift[3] = {0};
//     interp->Interpolate(-1, shift, block_coarse_, data_coarse_, block_fine_, data_fine_);

//     // check the result

//     for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
//         for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
//             for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
//                 real_t x = ((real_t)i0 + 0.5) * hfine_;
//                 real_t y = ((real_t)i1 + 0.5) * hfine_;
//                 real_t z = ((real_t)i2 + 0.5) * hfine_;
//                 EXPECT_NEAR(data_fine_[m_midx(i0, i1, i2, 0, block_fine_)], x * x + y * y + z * z, DOUBLE_TOL) << "fails for " << i0 << " " << i1 << " " << i2 << " hgrid = " << hfine_;
//             }
//         }
//     }
// }

//==============================================================================================================================
TEST_F(valid_Wavelet, refine_order_5_linear) {
    // fill the source
    for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
        for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
            for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hcoarse_;
                real_t y = ((real_t)i1 + 0.5) * hcoarse_;
                real_t z = ((real_t)i2 + 0.5) * hcoarse_;

                data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = x + y + z;
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet<5>* interp   = new Wavelet<5>();
    lid_t       shift[3] = {0};
    interp->Interpolate(-1, shift, block_coarse_, data_coarse_, block_fine_, data_fine_);

    // check the result

    for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
        for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
            for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hfine_;
                real_t y = ((real_t)i1 + 0.5) * hfine_;
                real_t z = ((real_t)i2 + 0.5) * hfine_;
                EXPECT_NEAR(data_fine_[m_midx(i0, i1, i2, 0, block_fine_)], x + y + z, DOUBLE_TOL) << "fails for " << i0 << " " << i1 << " " << i2 << " hgrid = " << hfine_;
            }
        }
    }
}

//==============================================================================================================================
// TEST_F(valid_Wavelet, order_5_quad) {
//     // fill the source
//     for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
//         for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
//             for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
//                 real_t x = ((real_t)i0 + 0.5) * hcoarse_;
//                 real_t y = ((real_t)i1 + 0.5) * hcoarse_;
//                 real_t z = ((real_t)i2 + 0.5) * hcoarse_;

//                 data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = x * x + y * y + z * z;
//             }
//         }
//     }

//     // do the interpolation (-1 is refinement)
//     Wavelet<5>* interp   = new Wavelet<5>();
//     lid_t       shift[3] = {0};
//     interp->Interpolate(-1, shift, block_coarse_, data_coarse_, block_fine_, data_fine_);

//     // check the result

//     for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
//         for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
//             for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
//                 real_t x = ((real_t)i0 + 0.5) * hfine_;
//                 real_t y = ((real_t)i1 + 0.5) * hfine_;
//                 real_t z = ((real_t)i2 + 0.5) * hfine_;
//                 EXPECT_NEAR(data_fine_[m_midx(i0, i1, i2, 0, block_fine_)], x * x + y * y + z * z, DOUBLE_TOL) << "fails for " << i0 << " " << i1 << " " << i2 << " hgrid = " << hfine_;
//             }
//         }
//     }
// }

//==============================================================================================================================
TEST_F(valid_Wavelet, detail_order_3) {
    // fill the source
    for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
        for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
            for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hcoarse_;
                real_t y = ((real_t)i1 + 0.5) * hcoarse_;
                real_t z = ((real_t)i2 + 0.5) * hcoarse_;

                data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = x + y + z;
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet<3>* interp     = new Wavelet<3>();
    real_t      detail_max = 0;
    interp->Criterion(block_coarse_, data_coarse_, &detail_max);
    ASSERT_NEAR(detail_max, 0.0, DOUBLE_TOL);
}
//==============================================================================================================================
TEST_F(valid_Wavelet, detail_order_5) {
    // fill the source
    for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
        for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
            for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                real_t x = ((real_t)i0 + 0.5) * hcoarse_;
                real_t y = ((real_t)i1 + 0.5) * hcoarse_;
                real_t z = ((real_t)i2 + 0.5) * hcoarse_;

                data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = x + y + z;
            }
        }
    }

    // do the interpolation (-1 is refinement)
    Wavelet<5>* interp     = new Wavelet<5>();
    real_t      detail_max = 0;
    interp->Criterion(block_coarse_, data_coarse_, &detail_max);

    ASSERT_NEAR(detail_max, 0.0, DOUBLE_TOL);
}