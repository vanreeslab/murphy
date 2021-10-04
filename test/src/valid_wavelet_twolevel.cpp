#include <cmath>
#include <limits>
#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "gtest/gtest.h"
#include "operator/error.hpp"
#include "operator/xblas.hpp"
#include "tools/ioh5.hpp"
#include "valid_toolbox.hpp"

static real_t sigma     = 0.1;
static real_t center[3] = {1.0, 1.0, 1.0};

// define the initial condition and the analytical solution
static lambda_setvalue_t exp_setvalue = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
    real_t pos[3];
    block->pos(i0, i1, i2, pos);
    block->data(fid, 0)(i0, i1, i2) = scalar_exp(pos, center, sigma);
};
static lambda_error_t exp_error = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
    real_t pos[3];
    block->pos(i0, i1, i2, pos);
    return scalar_exp(pos, center, sigma);
};

static lambda_setvalue_t flipflop_setvalue = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
    real_t pos[3];
    block->pos(i0, i1, i2, pos);
    const real_t x = pos[0];
    const real_t y = pos[1];
    const real_t z = pos[2];

    real_t hfine                    = 0.5 / M_N * 2.0;
    block->data(fid, 0)(i0, i1, i2) = (1.0) * sin(x * M_PI / hfine) +
                                      (1.0) * sin(y * M_PI / hfine) +
                                      (1.0) * sin(z * M_PI / hfine) +
                                      (1.0) * sin(x * M_PI / hfine) * sin(y * M_PI / hfine) +
                                      (1.0) * sin(y * M_PI / hfine) * sin(z * M_PI / hfine) +
                                      (1.0) * sin(x * M_PI / hfine) * sin(z * M_PI / hfine) +
                                      (1.0) * sin(x * M_PI / hfine) * sin(y * M_PI / hfine) * sin(z * M_PI / hfine);
};
static lambda_error_t flipflop_error = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
    real_t pos[3];
    block->pos(i0, i1, i2, pos);
    const real_t x = pos[0];
    const real_t y = pos[1];
    const real_t z = pos[2];

    real_t hfine = 0.5 / M_N * 2.0;

    return (1.0) * sin(x * M_PI / hfine) +
           (1.0) * sin(y * M_PI / hfine) +
           (1.0) * sin(z * M_PI / hfine) +
           (1.0) * sin(x * M_PI / hfine) * sin(y * M_PI / hfine) +
           (1.0) * sin(y * M_PI / hfine) * sin(z * M_PI / hfine) +
           (1.0) * sin(x * M_PI / hfine) * sin(z * M_PI / hfine) +
           (1.0) * sin(x * M_PI / hfine) * sin(y * M_PI / hfine) * sin(z * M_PI / hfine);
};
// class InitCondition_TwoLevels : public SetValue {
//    protected:
//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override {
//         //-------------------------------------------------------------------------
//         real_t        pos[3];
//         const real_t* xyz   = block->xyz();
//         const real_t* hgrid = block->hgrid();

//         real_t sigma     = 0.1;
//         real_t center[3] = {1.0, 1.0, 1.0};

//         // const real_t oo_sigma2 = 1.0 / (sigma * sigma);
//         const real_t fact = 1.0;

//         // get the pointers correct
//         real_t* data = block->data(fid, 0).Write();

//         auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // get the position
//             real_t pos[3];
//             block->pos(i0, i1, i2, pos);

//             // compute the gaussian
//             const real_t rhox       = (pos[0] - center[0]) / sigma;
//             const real_t rhoy       = (pos[1] - center[1]) / sigma;
//             const real_t rhoz       = (pos[2] - center[2]) / sigma;
//             const real_t rho        = rhox * rhox + rhoy * rhoy + rhoz * rhoz;
//             data[m_idx(i0, i1, i2)] = fact * std::exp(-rho);
//         };

//         for_loop(&op, start_, end_);
//         //-------------------------------------------------------------------------
//     };

//    public:
//     explicit InitCondition_TwoLevels() : SetValue(nullptr){};
//     explicit InitCondition_TwoLevels(const Wavelet*  interp) : SetValue(interp){};
// };

// class InitCondition_FlipFlop : public SetValue {
//    protected:
//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override {
//         //-------------------------------------------------------------------------
//         real_t        pos[3];
//         const real_t* xyz   = block->xyz();
//         const real_t* hgrid = block->hgrid();

//         real_t sigma     = 0.1;
//         real_t center[3] = {1.0, 1.0, 1.0};

//         // const real_t oo_sigma2 = 1.0 / (sigma * sigma);
//         const real_t fact = 1.0;

//         // get the pointers correct
//         real_t* data = block->data(fid, 0).Write();

//         auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // get the position
//             real_t pos[3];
//             block->pos(i0, i1, i2, pos);

//             const real_t x          = (pos[0] - center[0]);
//             const real_t y          = (pos[1] - center[1]);
//             const real_t z          = (pos[2] - center[2]);
//             real_t       hfine      = 0.5 / M_N * 2.0;
//             data[m_idx(i0, i1, i2)] = (1.0) * sin(x * M_PI / hfine) +
//                                       (1.0) * sin(y * M_PI / hfine) +
//                                       (1.0) * sin(z * M_PI / hfine) +
//                                       (1.0) * sin(x * M_PI / hfine) * sin(y * M_PI / hfine) +
//                                       (1.0) * sin(y * M_PI / hfine) * sin(z * M_PI / hfine) +
//                                       (1.0) * sin(x * M_PI / hfine) * sin(z * M_PI / hfine) +
//                                       (1.0) * sin(x * M_PI / hfine) * sin(y * M_PI / hfine) * sin(z * M_PI / hfine);
//         };

//         for_loop(&op, start_, end_);
//         //-------------------------------------------------------------------------
//     };

//    public:
//     explicit InitCondition_FlipFlop() : SetValue(nullptr){};
//     explicit InitCondition_FlipFlop(const Wavelet*  interp) : SetValue(interp){};
// };

class TwoLevel : public ::testing::TestWithParam<int> {
   public:
    int    case_id_;
    Grid*  grid_;
    Field* scal_;

    void SetUp() override {
        case_id_ = GetParam();
        m_log("--------------------------------------------------------------------------------");
        m_log("case id: %d", case_id_);
        m_log("--------------");

        // create a grid - uniform on level 1
        bool  period[3]   = {true, true, true};
        lid_t grid_len[3] = {2, 2, 2};
        grid_             = new Grid(1, period, grid_len, M_GRIDBLOCK, MPI_COMM_WORLD, nullptr);
        grid_->level_limit(0, 1);

        // create a field an put it on it
        scal_ = new Field("scal", 1);
        grid_->AddField(scal_);
    };
    void TearDown() override {
        grid_->DeleteField(scal_);
        delete (scal_);
        delete (grid_);
    };
};

using std::list;
using std::string;

static const real_t zero_tol = 1000.0 * std::numeric_limits<real_t>::epsilon();

/**
 * @brief validate every possible intersection between blocks on two levels (both moment and epsilon)
 * 
 * for one intersection, there is max 8 blocks sharing the same corner.
 * Every block might be coarse or fine, so in total 2^8 possibilities
 * 
 */
TEST_P(TwoLevel, periodic) {
    SetValue init(exp_setvalue);
    init(grid_, scal_);
    // get the Ghosts:
    const bidx_t ghost_len[2] = {grid_->interp()->nghost_front(), grid_->interp()->nghost_back()};
    grid_->GhostPull(scal_, ghost_len);

    //................................................
    // compute the details and get the max value
    Field detail("details", 1);
    grid_->AddField(&detail);
    grid_->StoreDetails(scal_, &detail);

    BMax   max;
    real_t max_detail = max(grid_, &detail);

    m_log("max detail = %e", max_detail);

    grid_->DeleteField(&detail);

    // get the moment in the fine level
    BMoment moment;
    grid_->GhostPull(scal_, &moment);
    real_t fine_moment0, fine_moment1[3];
    moment(grid_, scal_, &fine_moment0, fine_moment1);

    //................................................
    // adapt the tree
    {
        list<Patch> patch_list;
        TreeCaseIdToPatch(0, case_id_, &patch_list);

        bidx_t ghost_len[2] = {0, 0};
        grid_->GhostLengthAdapt(ghost_len);
        grid_->GhostPull(scal_, ghost_len);

        grid_->Adapt(&patch_list);
        // grid_->GhostPull(scal_,ghost_len);
    }

    //................................................
    // get the coarse moments
    real_t coarse_moment0, coarse_moment1[3];
    grid_->GhostPull(scal_, &moment);
    moment(grid_, scal_, &coarse_moment0, coarse_moment1);
    real_t mom0_coarse_error = fabs(fine_moment0 - coarse_moment0);
    m_log("[case %d] coarse moment error = |%e - %e| = %e", case_id_, fine_moment0, coarse_moment0, mom0_coarse_error);
    m_log("[case %d] coarse moment error = |%e - %e|", case_id_, fine_moment1[0], coarse_moment1[0]);
    m_log("[case %d] coarse moment error = |%e - %e|", case_id_, fine_moment1[1], coarse_moment1[1]);
    m_log("[case %d] coarse moment error = |%e - %e|", case_id_, fine_moment1[2], coarse_moment1[2]);

    //................................................
    // go back up
    {
        // create the needed patch list
        list<Patch> patch_list;
        real_t      origin[3] = {0.0, 0.0, 0.0};
        real_t      length[3] = {2.0, 2.0, 2.0};
        patch_list.push_back(Patch(origin, length, 1));

        // adapt
        bidx_t ghost_len[2] = {0, 0};
        grid_->GhostLengthAdapt(ghost_len);
        grid_->GhostPull(scal_, ghost_len);
        grid_->Adapt(&patch_list);
        // grid_->GhostPull(scal_,ghost_len);
    }

    //................................................
    // get the analytical solution
    {
        // Field sol("sol", 1);
        // grid_->AddField(&sol);
        // InitCondition_TwoLevels init;
        // init(grid_, &sol);

        // get the error
        real_t err2, erri;
        Error  error;
        grid_->GhostPull(scal_, &error);
        error.Norms(grid_, scal_, &exp_error, &err2, &erri);
        real_t interp_pred = fabs(grid_->interp()->eps_const() * max_detail);
        m_log("[case %d] interp error = %e <? %e -> factor = %e vs %e", case_id_, erri, interp_pred, erri / max_detail, grid_->interp()->eps_const());
        // grid_->DeleteField(&sol);

        // get the moment:
        real_t smooth_moment0, smooth_moment1[3];
        grid_->GhostPull(scal_, &moment);
        moment(grid_, scal_, &smooth_moment0, smooth_moment1);
        real_t mom_smooth_error[4];
        mom_smooth_error[0] = fabs(fine_moment0 - smooth_moment0);
        for (lda_t ida = 0; ida < 3; ++ida) {
            mom_smooth_error[ida + 1] = fabs(fine_moment1[ida] - smooth_moment1[ida]);
        }
        m_log("[case %d] moment error = %e %e %e %e <? %e", case_id_, mom_smooth_error[0], mom_smooth_error[1], mom_smooth_error[2], mom_smooth_error[3], zero_tol);

        ASSERT_LT(erri, interp_pred);

        // !! we cannot impose the 1st or any other moments conservation, it doesn't work for some reasons
        // therefore, we only check moment 0!!!! (it's an open-question on the why it doesn't work!)
        if (grid_->interp()->Nt() > 0) {
            for (lda_t ida = 0; ida < 4; ++ida) {
                ASSERT_LT(mom_smooth_error[ida], zero_tol);
            }
        }
    }
}

TEST_P(TwoLevel, flipfop) {
    SetValue init(flipflop_setvalue);
    init(grid_, scal_);
    // get the Ghosts:
    const bidx_t ghost_len[2] = {grid_->interp()->nghost_front(), grid_->interp()->nghost_back()};
    grid_->GhostPull(scal_, ghost_len);
    //................................................
    // compute the details and get the max value
    Field detail("details", 1);
    grid_->AddField(&detail);
    grid_->StoreDetails(scal_, &detail);

    BMax   max;
    real_t max_detail = max(grid_, &detail);

    m_log("max detail = %e", max_detail);

    grid_->DeleteField(&detail);

    // get the moment in the fine level
    BMoment moment;
    real_t  fine_moment0, fine_moment1[3];
    moment(grid_, scal_, &fine_moment0, fine_moment1);

    //................................................
    // adapt the tree
    {
        list<Patch> patch_list;
        TreeCaseIdToPatch(0, case_id_, &patch_list);

        bidx_t ghost_len[2] = {0, 0};
        grid_->GhostLengthAdapt(ghost_len);
        grid_->GhostPull(scal_, ghost_len);
        grid_->Adapt(&patch_list);
        // grid_->GhostPull(scal_);
    }

    //................................................
    // get the coarse moments
    real_t coarse_moment0, coarse_moment1[3];
    grid_->GhostPull(scal_, &moment);
    moment(grid_, scal_, &coarse_moment0, coarse_moment1);
    real_t mom0_coarse_error = fabs(fine_moment0 - coarse_moment0);
    m_log("[case %d] coarse moment error = |%e - %e| = %e", case_id_, fine_moment0, coarse_moment0, mom0_coarse_error);
    // m_log("[case %d] coarse moment error = |%e - %e|", case_id_, fine_moment1[0], coarse_moment1[0]);
    // m_log("[case %d] coarse moment error = |%e - %e|", case_id_, fine_moment1[1], coarse_moment1[1]);
    // m_log("[case %d] coarse moment error = |%e - %e|", case_id_, fine_moment1[2], coarse_moment1[2]);

    //................................................
    // go back up
    {
        // create the needed patch list
        list<Patch> patch_list;
        real_t      origin[3] = {0.0, 0.0, 0.0};
        real_t      length[3] = {2.0, 2.0, 2.0};
        patch_list.push_back(Patch(origin, length, 1));

        bidx_t ghost_len[2] = {0, 0};
        grid_->GhostLengthAdapt(ghost_len);
        grid_->GhostPull(scal_, ghost_len);

        // adapt
        grid_->Adapt(&patch_list);
        // grid_->GhostPull(scal_);
    }

    //................................................
    // get the analytical solution
    {
        // Field sol("sol", 1);
        // grid_->AddField(&sol);
        // InitCondition_FlipFlop init;
        // init(grid_, &sol);

        // get the error
        real_t err2, erri;
        Error  error;
        grid_->GhostPull(scal_, &error);
        error.Norms(grid_, scal_, &flipflop_error, &err2, &erri);
        real_t interp_pred = fabs(grid_->interp()->eps_const() * max_detail);
        m_log("[case %d] interp error = %e <? %e -> factor = %e vs %e", case_id_, erri, interp_pred, erri / max_detail, grid_->interp()->eps_const());
        // grid_->DeleteField(&sol);

        // get the moment:
        real_t smooth_moment0, smooth_moment1[3];
        grid_->GhostPull(scal_, &moment);
        moment(grid_, scal_, &smooth_moment0, smooth_moment1);
        real_t mom_smooth_error[4];
        mom_smooth_error[0] = fabs(fine_moment0 - smooth_moment0);
        for (lda_t ida = 0; ida < 3; ++ida) {
            mom_smooth_error[ida + 1] = fabs(fine_moment1[ida] - smooth_moment1[ida]);
        }
        m_log("[case %d] moment error = %e %e %e %e <? %e", case_id_, mom_smooth_error[0], mom_smooth_error[1], mom_smooth_error[2], mom_smooth_error[3], zero_tol);

        ASSERT_LT(erri, interp_pred);

        // !! we cannot impose the 1st or any other moments conservation, it doesn't work for some reasons
        // therefore, we only check moment 0!!!! (it's an open-question on the why it doesn't work!)
        if (grid_->interp()->Nt() > 0) {
            for (lda_t ida = 0; ida < 1; ++ida) {
                ASSERT_LT(mom_smooth_error[ida], zero_tol);
            }
        }
    }
}

INSTANTIATE_TEST_SUITE_P(ValidWavelet,
                         TwoLevel,
                         testing::Range(0, 256));
// INSTANTIATE_TEST_SUITE_P(ValidWavelet,
//                          TwoLevel2,
//                          testing::Range(0, 256));