
#include "core/macros.hpp"
#include "core/types.hpp"
#include "error.hpp"
#include "grid/subblock.hpp"
#include "gtest/gtest.h"
#include "wavelet/interpolating_wavelet.hpp"

/**
 * @brief test the wavelet's filters
 * 
 */
class ValidWaveletKernel : public ::testing::Test {
//    protected:
//     SubBlock* block_coarse_;
//     SubBlock* block_fine_;
//     real_p    data_coarse_;
//     real_p    data_fine_;
//     real_t    hcoarse_;
//     real_t    hfine_;

//     lid_t coarse_start_[3];
//     lid_t coarse_end_[3];
//     lid_t fine_start_[3];
//     lid_t fine_end_[3];

    void SetUp() override{};
    void TearDown() override{};
};

TEST_F(ValidWaveletKernel, filter_length) {
    InterpolatingWavelet interp;
    m_log("testing wavelet %d.%d", interp.N(), interp.Nt());

// if (interp.N() == 2 && interp.Nt() == 0) {
#if (M_WAVELET_NT == 2 && M_WAVELET_N == 0)
    // refine
    EXPECT_EQ(interp.nghost_front_refine(), 0);
    EXPECT_EQ(interp.nghost_back_refine(), 1);
    // coarsen
    EXPECT_EQ(interp.nghost_front_coarsen(), 0);
    EXPECT_EQ(interp.nghost_back_coarsen(), 0);
    // criterion
    EXPECT_EQ(interp.nghost_front_criterion_smooth(), 0);
    EXPECT_EQ(interp.nghost_back_criterion_smooth(), 1);
    // details
    EXPECT_EQ(interp.ndetail_citerion_extend_front(), 0);
    EXPECT_EQ(interp.ndetail_citerion_extend_back(), 1);
    // details
    EXPECT_EQ(interp.ndetail_smooth_extend_front(), 0);
    EXPECT_EQ(interp.ndetail_smooth_extend_back(), 1);
#endif
// } else if (interp.N() == 2 && interp.Nt() == 2) {
#if (M_WAVELET_NT == 2 && M_WAVELET_N == 2)
    // refine
    EXPECT_EQ(interp.nghost_front_refine(), 0);
    EXPECT_EQ(interp.nghost_back_refine(), 1);
    // coarsen
    EXPECT_EQ(interp.nghost_front_coarsen(), 2);
    EXPECT_EQ(interp.nghost_back_coarsen(), 1);
    // criterion
    EXPECT_EQ(interp.nghost_front_criterion_smooth(), 2);
    EXPECT_EQ(interp.nghost_back_criterion_smooth(), 3);
    // details
    EXPECT_EQ(interp.ndetail_citerion_extend_front(), 1);
    EXPECT_EQ(interp.ndetail_citerion_extend_back(), 2);
    // details
    EXPECT_EQ(interp.ndetail_smooth_extend_front(), 1);
    EXPECT_EQ(interp.ndetail_smooth_extend_back(), 2);
#endif
    // } else if (interp.N() == 4 && interp.Nt() == 0) {
#if (M_WAVELET_NT == 4 && M_WAVELET_N == 0)
    // refine
    EXPECT_EQ(interp.nghost_front_refine(), 1);
    EXPECT_EQ(interp.nghost_back_refine(), 2);
    // coarsen
    EXPECT_EQ(interp.nghost_front_coarsen(), 0);
    EXPECT_EQ(interp.nghost_back_coarsen(), 0);
    // criterion
    EXPECT_EQ(interp.nghost_front_criterion_smooth(), 4);
    EXPECT_EQ(interp.nghost_back_criterion_smooth(), 5);
    // details
    EXPECT_EQ(interp.ndetail_citerion_extend_front(), 2);
    EXPECT_EQ(interp.ndetail_citerion_extend_back(), 3);
    // details
    EXPECT_EQ(interp.ndetail_smooth_extend_front(), 2);
    EXPECT_EQ(interp.ndetail_smooth_extend_back(), 3);
#endif
    // } else if (interp.N() == 4 && interp.Nt() == 2) {
#if (M_WAVELET_NT == 4 && M_WAVELET_N == 2)
    // refine
    EXPECT_EQ(interp.nghost_front_refine(), 1);
    EXPECT_EQ(interp.nghost_back_refine(), 2);
    // coarsen
    EXPECT_EQ(interp.nghost_front_coarsen(), 4);
    EXPECT_EQ(interp.nghost_back_coarsen(), 3);
    // criterion
    EXPECT_EQ(interp.nghost_front_criterion_smooth(), 6);
    EXPECT_EQ(interp.nghost_back_criterion_smooth(), 7);
    // details
    EXPECT_EQ(interp.ndetail_citerion_extend_front(), 3);
    EXPECT_EQ(interp.ndetail_citerion_extend_back(), 4);
    // details
    EXPECT_EQ(interp.ndetail_smooth_extend_front(), 3);
    EXPECT_EQ(interp.ndetail_smooth_extend_back(), 4);
#endif
    // } else if (interp.N() == 6 && interp.Nt() == 0) {
#if (M_WAVELET_NT == 6 && M_WAVELET_N == 0)
    // refine
    EXPECT_EQ(interp.nghost_front_refine(), 2);
    EXPECT_EQ(interp.nghost_back_refine(), 3);
    // coarsen
    EXPECT_EQ(interp.nghost_front_coarsen(), 0);
    EXPECT_EQ(interp.nghost_back_coarsen(), 0);
    // criterion
    EXPECT_EQ(interp.nghost_front_criterion_smooth(), 8);
    EXPECT_EQ(interp.nghost_back_criterion_smooth(), 9);
    // details
    EXPECT_EQ(interp.ndetail_citerion_extend_front(), 4);
    EXPECT_EQ(interp.ndetail_citerion_extend_back(), 5);
    // details
    EXPECT_EQ(interp.ndetail_smooth_extend_front(), 4);
    EXPECT_EQ(interp.ndetail_smooth_extend_back(), 5);
#endif
    // } else if (interp.N() == 6 && interp.Nt() == 2) {
#if (M_WAVELET_NT == 6 && M_WAVELET_N == 2)
    // refine
    EXPECT_EQ(interp.nghost_front_refine(), 2);
    EXPECT_EQ(interp.nghost_back_refine(), 3);
    // coarsen
    EXPECT_EQ(interp.nghost_front_coarsen(), 6);
    EXPECT_EQ(interp.nghost_back_coarsen(), 5);
    // criterion
    EXPECT_EQ(interp.nghost_front_criterion_smooth(), 10);
    EXPECT_EQ(interp.nghost_back_criterion_smooth(), 11);
    // details
    EXPECT_EQ(interp.ndetail_citerion_extend_front(), 5);
    EXPECT_EQ(interp.ndetail_citerion_extend_back(), 6);
    // details
    EXPECT_EQ(interp.ndetail_smooth_extend_front(), 5);
    EXPECT_EQ(interp.ndetail_smooth_extend_back(), 6);
#endif
    // }
}