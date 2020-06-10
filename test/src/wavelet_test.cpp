
#include "wavelet.hpp"

#include "error.hpp"
#include "gtest/gtest.h"
#include "murphy.hpp"
#include "subblock.hpp"
#include "mempool.hpp"

#define DOUBLE_TOL 1e-13

class valid_Wavelet : public ::testing::Test {
   protected:
    SubBlock* block_coarse_;
    SubBlock* block_fine_;
    real_p    data_coarse_;
    real_p    data_fine_;
    real_t    hcoarse_;
    real_t    hfine_;

    MemPool* mem_pool_;

    lid_t coarse_start_[3];
    lid_t coarse_end_[3];
    lid_t fine_start_[3];
    lid_t fine_end_[3];



    void SetUp() override {
        data_fine_   = (real_t*)m_calloc(M_STRIDE * M_STRIDE * M_STRIDE * sizeof(real_t));
        data_coarse_ = (real_t*)m_calloc(M_STRIDE * M_STRIDE * M_STRIDE * sizeof(real_t));
        
        mem_pool_ = new MemPool();
        
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
        real_p data_fine   = data_fine_ - m_zeroidx(0, block_fine_);
        real_p data_coarse = data_coarse_ - m_zeroidx(0, block_coarse_);
        m_free(data_fine);
        m_free(data_coarse);

        delete (block_coarse_);
        delete (block_fine_);
        delete (mem_pool_);
    };
};

//==============================================================================================================================
static double order3_int(const double x) {
    return M_PI * (x * x * x / 3.0) + M_PI_2 * (x * x / 2.0) - M_SQRT2 * x;
}
static double poly_1(const double x) {
    return M_PI_2 * x - M_SQRT2;
}
TEST_F(valid_Wavelet, refine_order_3) {
    // real_t coarse_mom_1 = 0.0;
    // real_t coarse_mom_2 = 0.0;

    // for (int id = 0; id < 3; id++) {
    //     // fill the source
    //     for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
    //         for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
    //             for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
    //                 real_t x       = ((real_t)i0 + 0.5) * hcoarse_;
    //                 real_t y       = ((real_t)i1 + 0.5) * hcoarse_;
    //                 real_t z       = ((real_t)i2 + 0.5) * hcoarse_;
    //                 real_t pos[3]  = {x, y, z};
    //                 real_t x_left  = pos[id] - hcoarse_ / 2.0;
    //                 real_t x_right = pos[id] + hcoarse_ / 2.0;

    //                 data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = (order3_int(x_right) - order3_int(x_left)) / hcoarse_;
    //             }
    //         }
    //     }

    //     // do the interpolation (-1 is refinement)
    //     Wavelet<3>* interp   = new Wavelet<3>();
    //     lid_t       shift[3] = {0};
    //     interp->Interpolate(-1, shift, block_coarse_, data_coarse_, block_fine_, data_fine_);

    //     // check the result
    //     real_t fine_mom_1 = 0.0;
    //     real_t fine_mom_2 = 0.0;
    //     for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
    //         for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
    //             for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
    //                 real_t x       = ((real_t)i0 + 0.5) * hfine_;
    //                 real_t y       = ((real_t)i1 + 0.5) * hfine_;
    //                 real_t z       = ((real_t)i2 + 0.5) * hfine_;
    //                 real_t pos[3]  = {x, y, z};
    //                 real_t x_left  = pos[id] - hfine_ / 2.0;
    //                 real_t x_right = pos[id] + hfine_ / 2.0;

    //                 ASSERT_NEAR(data_fine_[m_midx(i0, i1, i2, 0, block_fine_)], (order3_int(x_right) - order3_int(x_left)) / hfine_, DOUBLE_TOL) << "while testing direction " << id;
    //             }
    //         }
    //     }
    // }
}
TEST_F(valid_Wavelet, coarsen_1_order_2_2) {
    for (int id = 0; id < 3; id++) {
        coarse_start_[id] = 0;
        coarse_end_[id]   = M_HN;
        fine_start_[id]   = -M_GS;
        fine_end_[id]     = M_N + M_GS;
    }
    block_coarse_->Reset(M_GS, M_HN + 2 * M_GS, coarse_start_, coarse_end_);
    block_fine_->Reset(M_GS, M_STRIDE, fine_start_, fine_end_);

    for (int id = 0; id < 3; id++) {
        // fill the source
        for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
            for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
                for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
                    real_t x      = i0 * hfine_;
                    real_t y      = i1 * hfine_;
                    real_t z      = i2 * hfine_;
                    real_t pos[3] = {x, y, z};

                    data_fine_[m_midx(i0, i1, i2, 0, block_fine_)] = poly_2(pos[id]);
                }
            }
        }
        
        // do the coarsening
        Wavelet<2, 2>* interp   = new Wavelet<2, 2>();
        lid_t          shift[3] = {0};
        interp->Interpolate(1, shift, block_fine_, data_fine_, block_coarse_, data_coarse_, mem_pool_);

        for (sid_t it = 0; it < omp_get_max_threads(); it++) {
            ASSERT_EQ(mem_pool_->n_locked(it), 0);
        }

        //check the result
        for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
            for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
                for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
                    real_t x      = i0 * hcoarse_;
                    real_t y      = i1 * hcoarse_;
                    real_t z      = i2 * hcoarse_;
                    real_t pos[3] = {x, y, z};

                    ASSERT_NEAR(data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)], poly_1(pos[id]), DOUBLE_TOL) << "while testing direction " << id;
                }
            }
        }
        delete(interp);

    }
}
// //==============================================================================================================================
// static double order5_int(const double x) {
//     return M_2_SQRTPI*(x*x*x*x*x/5.0) - M_E*(x*x*x*x/4.0) + M_PI * (x * x * x / 3.0) + M_PI_2 * (x * x / 2.0) - M_SQRT2 * x;
// }
// TEST_F(valid_Wavelet, refine_order_5) {
//     for (int id = 0; id < 3; id++) {
//         // fill the source
//         for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
//             for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
//                 for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
//                     real_t x       = ((real_t)i0 + 0.5) * hcoarse_;
//                     real_t y       = ((real_t)i1 + 0.5) * hcoarse_;
//                     real_t z       = ((real_t)i2 + 0.5) * hcoarse_;
//                     real_t pos[3]  = {x, y, z};
//                     real_t x_left  = pos[id] - hcoarse_ / 2.0;
//                     real_t x_right = pos[id] + hcoarse_ / 2.0;

//                     data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)] = (order5_int(x_right) - order5_int(x_left)) / hcoarse_;
//                 }
//             }
//         }

//         // do the refinement (-1 is refinement)
//         Wavelet<5>* interp   = new Wavelet<5>();
//         lid_t       shift[3] = {0};
//         interp->Interpolate(-1, shift, block_coarse_, data_coarse_, block_fine_, data_fine_);

//         // check the result
//         for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
//             for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
//                 for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
//                     real_t x       = ((real_t)i0 + 0.5) * hfine_;
//                     real_t y       = ((real_t)i1 + 0.5) * hfine_;
//                     real_t z       = ((real_t)i2 + 0.5) * hfine_;
//                     real_t pos[3]  = {x, y, z};
//                     real_t x_left  = pos[id] - hfine_ / 2.0;
//                     real_t x_right = pos[id] + hfine_ / 2.0;

//                     ASSERT_NEAR(data_fine_[m_midx(i0, i1, i2, 0, block_fine_)], (order5_int(x_right) - order5_int(x_left)) / hfine_, DOUBLE_TOL) << "while testing direction " << id;
//                 }
//             }
//         }
//     }
// }
// TEST_F(valid_Wavelet, coarsen_order_5) {
//     for (int id = 0; id < 3; id++) {
//         coarse_start_[id] = 0;
//         coarse_end_[id]   = M_HN;
//         fine_start_[id]   = -M_GS;
//         fine_end_[id]     = M_N + M_GS;
//     }
//     block_coarse_->Reset(M_GS, M_HN + 2 * M_GS, coarse_start_, coarse_end_);
//     block_fine_->Reset(M_GS, M_STRIDE, fine_start_, fine_end_);

//     for (int id = 0; id < 3; id++) {
//         // fill the source
//         for (int i2 = fine_start_[2]; i2 < fine_end_[2]; i2++) {
//             for (int i1 = fine_start_[1]; i1 < fine_end_[1]; i1++) {
//                 for (int i0 = fine_start_[0]; i0 < fine_end_[0]; i0++) {
//                     real_t x       = ((real_t)i0 + 0.5) * hfine_;
//                     real_t y       = ((real_t)i1 + 0.5) * hfine_;
//                     real_t z       = ((real_t)i2 + 0.5) * hfine_;
//                     real_t pos[3]  = {x, y, z};
//                     real_t x_left  = pos[id] - hfine_ / 2.0;
//                     real_t x_right = pos[id] + hfine_ / 2.0;

//                     data_fine_[m_midx(i0, i1, i2, 0, block_fine_)] = (order5_int(x_right) - order5_int(x_left)) / hfine_;
//                 }
//             }
//         }

//         // do the coarsening
//         Wavelet<5>* interp   = new Wavelet<5>();
//         lid_t       shift[3] = {0};
//         interp->Interpolate(1, shift, block_fine_, data_fine_, block_coarse_, data_coarse_);

//         //check the result
//         for (int i2 = coarse_start_[2]; i2 < coarse_end_[2]; i2++) {
//             for (int i1 = coarse_start_[1]; i1 < coarse_end_[1]; i1++) {
//                 for (int i0 = coarse_start_[0]; i0 < coarse_end_[0]; i0++) {
//                     real_t x       = ((real_t)i0 + 0.5) * hcoarse_;
//                     real_t y       = ((real_t)i1 + 0.5) * hcoarse_;
//                     real_t z       = ((real_t)i2 + 0.5) * hcoarse_;
//                     real_t pos[3]  = {x, y, z};
//                     real_t x_left  = pos[id] - hcoarse_ / 2.0;
//                     real_t x_right = pos[id] + hcoarse_ / 2.0;

//                     ASSERT_NEAR(data_coarse_[m_midx(i0, i1, i2, 0, block_coarse_)], (order5_int(x_right) - order5_int(x_left)) / hcoarse_, DOUBLE_TOL) << "while testing direction " << id;
//                 }
//             }
//         }
//     }
// }

// //==============================================================================================================================
// TEST_F(valid_Wavelet, detail_order_3) {
//     // fill the source
//     for (int i2 = fine_start_[2] - block_fine_->gs(); i2 < fine_end_[2] + block_fine_->gs(); i2++) {
//         for (int i1 = fine_start_[1] - block_fine_->gs(); i1 < fine_end_[1] + block_fine_->gs(); i1++) {
//             for (int i0 = fine_start_[0] - block_fine_->gs(); i0 < fine_end_[0] + block_fine_->gs(); i0++) {
//                 real_t x = ((real_t)i0 + 0.5) * hfine_;
//                 real_t y = ((real_t)i1 + 0.5) * hfine_;
//                 real_t z = ((real_t)i2 + 0.5) * hfine_;

//                 data_fine_[m_midx(i0, i1, i2, 0, block_fine_)] = x + y + z;
//             }
//         }
//     }

//     // do the interpolation (-1 is refinement)
//     Wavelet<3>* interp     = new Wavelet<3>();
//     real_t      detail_max[8] = {0};
//     interp->Details(block_fine_, data_fine_, detail_max);
    
//     ASSERT_NEAR(detail_max[0], hfine_, DOUBLE_TOL); // d_x = 1
//     ASSERT_NEAR(detail_max[1], hfine_, DOUBLE_TOL); // d_y = 1
//     ASSERT_NEAR(detail_max[2], hfine_, DOUBLE_TOL); // d_z = 1
// }
