
#include "core/macros.hpp"
#include "core/types.hpp"
#include "error.hpp"
#include "gtest/gtest.h"
#include "interpolating_wavelet.hpp"
#include "subblock.hpp"

/**
 * @brief test the wavelet's filters
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

    void SetUp() override{};
    void TearDown() override{};
};

TEST_F(valid_Wavelet_Kernel, filter_length) {
    InterpolatingWavelet interp;
    m_log("testing wavelet %d.%d",interp.N(),interp.Nt());

    if (interp.N() == 2 && interp.Nt() == 0) {
        // refine
        ASSERT_EQ(interp.nghost_front_refine(), 0);
        ASSERT_EQ(interp.nghost_back_refine(), 1);
        // coarsen
        ASSERT_EQ(interp.nghost_front_coarsen(), 0);
        ASSERT_EQ(interp.nghost_back_coarsen(), 0);
        // criterion
        ASSERT_EQ(interp.nghost_front_criterion_smooth(), 0);
        ASSERT_EQ(interp.nghost_back_criterion_smooth(), 1);
        // details
        ASSERT_EQ(interp.ndetail_citerion_extend_front(), 0);
        ASSERT_EQ(interp.ndetail_citerion_extend_back(), 1);
        // details
        ASSERT_EQ(interp.ndetail_smooth_extend_front(), 0);
        ASSERT_EQ(interp.ndetail_smooth_extend_back(), 1);
    } else if (interp.N() == 2 && interp.Nt() == 2) {
        // refine
        ASSERT_EQ(interp.nghost_front_refine(), 0);
        ASSERT_EQ(interp.nghost_back_refine(), 1);
        // coarsen
        ASSERT_EQ(interp.nghost_front_coarsen(), 2);
        ASSERT_EQ(interp.nghost_back_coarsen(), 1);
        // criterion
        ASSERT_EQ(interp.nghost_front_criterion_smooth(), 2);
        ASSERT_EQ(interp.nghost_back_criterion_smooth(), 3);
        // details
        ASSERT_EQ(interp.ndetail_citerion_extend_front(), 1);
        ASSERT_EQ(interp.ndetail_citerion_extend_back(), 2);
        // details
        ASSERT_EQ(interp.ndetail_smooth_extend_front(), 1);
        ASSERT_EQ(interp.ndetail_smooth_extend_back(), 2);
    } else if (interp.N() == 4 && interp.Nt() == 0) {
        // refine
        ASSERT_EQ(interp.nghost_front_refine(), 1);
        ASSERT_EQ(interp.nghost_back_refine(), 2);
        // coarsen
        ASSERT_EQ(interp.nghost_front_coarsen(), 0);
        ASSERT_EQ(interp.nghost_back_coarsen(), 0);
        // criterion
        ASSERT_EQ(interp.nghost_front_criterion_smooth(), 4);
        ASSERT_EQ(interp.nghost_back_criterion_smooth(), 5);
        // details
        ASSERT_EQ(interp.ndetail_citerion_extend_front(), 2);
        ASSERT_EQ(interp.ndetail_citerion_extend_back(), 3);
        // details
        ASSERT_EQ(interp.ndetail_smooth_extend_front(), 2);
        ASSERT_EQ(interp.ndetail_smooth_extend_back(), 3);
    } else if (interp.N() == 4 && interp.Nt() == 2) {
        // refine
        ASSERT_EQ(interp.nghost_front_refine(), 1);
        ASSERT_EQ(interp.nghost_back_refine(), 2);
        // coarsen
        ASSERT_EQ(interp.nghost_front_coarsen(), 4);
        ASSERT_EQ(interp.nghost_back_coarsen(), 3);
        // criterion
        ASSERT_EQ(interp.nghost_front_criterion_smooth(), 6);
        ASSERT_EQ(interp.nghost_back_criterion_smooth(), 7);
        // details
        ASSERT_EQ(interp.ndetail_citerion_extend_front(), 3);
        ASSERT_EQ(interp.ndetail_citerion_extend_back(), 4);
        // details
        ASSERT_EQ(interp.ndetail_smooth_extend_front(), 3);
        ASSERT_EQ(interp.ndetail_smooth_extend_back(), 4);
    } else if (interp.N() == 6 && interp.Nt() == 0) {
        // refine
        ASSERT_EQ(interp.nghost_front_refine(), 2);
        ASSERT_EQ(interp.nghost_back_refine(), 3);
        // coarsen
        ASSERT_EQ(interp.nghost_front_coarsen(), 0);
        ASSERT_EQ(interp.nghost_back_coarsen(), 0);
        // criterion
        ASSERT_EQ(interp.nghost_front_criterion_smooth(), 8);
        ASSERT_EQ(interp.nghost_back_criterion_smooth(), 9);
        // details
        ASSERT_EQ(interp.ndetail_citerion_extend_front(), 3);
        ASSERT_EQ(interp.ndetail_citerion_extend_back(), 4);
        // details
        ASSERT_EQ(interp.ndetail_smooth_extend_front(), 3);
        ASSERT_EQ(interp.ndetail_smooth_extend_back(), 4);
    } else if (interp.N() == 6 && interp.Nt() == 2) {
        // refine
        ASSERT_EQ(interp.nghost_front_refine(), 2);
        ASSERT_EQ(interp.nghost_back_refine(), 3);
        // coarsen
        ASSERT_EQ(interp.nghost_front_coarsen(), 6);
        ASSERT_EQ(interp.nghost_back_coarsen(), 5);
        // criterion
        ASSERT_EQ(interp.nghost_front_criterion_smooth(), 10);
        ASSERT_EQ(interp.nghost_back_criterion_smooth(), 11);
        // details
        ASSERT_EQ(interp.ndetail_citerion_extend_front(), 5);
        ASSERT_EQ(interp.ndetail_citerion_extend_back(), 6);
        // details
        ASSERT_EQ(interp.ndetail_smooth_extend_front(), 5);
        ASSERT_EQ(interp.ndetail_smooth_extend_back(), 6);
    }
}