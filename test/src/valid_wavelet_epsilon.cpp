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

// class InitCond_Epsilon : public SetValue {
//    protected:
//     void FillGridBlock(const qid_t*  qid, GridBlock*  block, Field*  fid) override {
//         //-------------------------------------------------------------------------
//         real_t        pos[3];
//         const real_t* xyz   = block->xyz();
//         const real_t* hgrid = block->hgrid();

//         real_t sigma     = 0.05;
//         real_t center[3] = {0.5, 0.5, 0.5};

//         // const real_t oo_sigma2 = 1.0 / (sigma * sigma);
//         const real_t fact = 1.0;

//         // get the pointers correct
//         real_t* data = block->data(fid, 0).Write();

//         auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // get the position
//             real_t pos[3];
//             block->pos(i0, i1, i2, pos);

//             // compute the gaussian
//             const real_t rhox = (pos[0] - center[0]) / sigma;
//             const real_t rhoy = (pos[1] - center[1]) / sigma;
//             const real_t rhoz = (pos[2] - center[2]) / sigma;
//             const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

//             data[m_idx(i0, i1, i2)] = fact * std::exp(-rho);
//         };

//         for_loop(&op, start_, end_);
//         //-------------------------------------------------------------------------
//     };

//    public:
//     explicit InitCond_Epsilon() : SetValue(nullptr){};
// };

static real_t sigma     = 0.05;
static real_t center[3] = {0.5, 0.5, 0.5};

using std::string;
using std::to_string;

// define the initial condition and the analytical solution
static lambda_error_t lambda_error = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
    // get the position
    real_t pos[3];
    block->pos(i0, i1, i2, pos);

    // call the function
    return scalar_exp(pos, center, sigma);
};

class Epsilon : public ::testing::TestWithParam<double> {
    void SetUp() override{};
    void TearDown() override{};
};

static const real_t zero_tol = 1000.0 * std::numeric_limits<real_t>::epsilon();

/**
 * @brief Check the epsilon + the moment condition on a "real-case" grid (periodic)
 * 
 */
TEST_P(Epsilon, periodic) {
    //-------------------------------------------------------------------------
    // get the parameters
    real_t  epsilon     = GetParam();
    level_t level_start = 4;

    // let's go
    bool  period[3]   = {true, true, true};
    lid_t grid_len[3] = {1, 1, 1};
    Grid  grid(level_start, period, grid_len, MPI_COMM_WORLD, nullptr);
    grid.level_limit(0, level_start);

    Field scal("scalar", 1);
    grid.AddField(&scal);
    scal.bctype(M_BC_EXTRAP);

    // InitCond_Epsilon init;
    // SetValue init(lambda_error,)
    // init(&grid, &scal);
    // initial condition
    lambda_setvalue_t lambda_initcond = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        real_t* data = block->data(fid).Write(i0, i1, i2);
        // set value
        data[0] = scalar_exp(pos, center, sigma);
    };
    SetValue init(lambda_initcond);
    init(&grid, &scal);

    grid.SetTol(epsilon * 1e+20, epsilon);

    // get the ghosting
    grid.GhostPull(&scal);

    // register the moments
    BMoment moment;
    real_t  fine_moment0, fine_moment1[3];
    moment(&grid, &scal, &fine_moment0, fine_moment1);

    // coarsen
    short_t count     = 1;
    level_t min_level = level_start;
    do {
        min_level = grid.MinLevel();
        grid.Coarsen(&scal);

        level_t tmp_min_lvl = grid.MinLevel();
        level_t tmp_max_lvl = grid.MaxLevel();
        m_log("Coarsening: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

        grid.GhostPull(&scal);
        real_t coarse_moment0, coarse_moment1[3];
        moment(&grid, &scal, &coarse_moment0, coarse_moment1);
        m_log("moments after coarsening: %e vs %e -> error = %e", fine_moment0, coarse_moment0, abs(fine_moment0 - coarse_moment0));

    } while (grid.MinLevel() < min_level && grid.MinLevel() > 0);

    // measure the moments
    grid.GhostPull(&scal);
    real_t coarse_moment0, coarse_moment1[3];
    moment(&grid, &scal, &coarse_moment0, coarse_moment1);
    m_log("moments after coarsening: %e vs %e -> error = %e", fine_moment0, coarse_moment0, abs(fine_moment0 - coarse_moment0));

    // and go back up
    std::list<Patch> patch;
    real_t           origin[3] = {0.0, 0.0, 0.0};
    real_t           length[3] = {(real_t)grid_len[0], (real_t)grid_len[1], (real_t)grid_len[2]};
    patch.push_back(Patch(origin, length, level_start));
    do {
        min_level = grid.MinLevel();
        // force the field refinement using a patch
        grid.GhostPull(&scal);
        grid.AdaptMagic(nullptr, &patch, nullptr, &cback_StatusCheck, nullptr, &cback_UpdateDependency, nullptr);

        level_t tmp_min_lvl = grid.MinLevel();
        level_t tmp_max_lvl = grid.MaxLevel();
        m_log("Refinement: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

    } while (grid.MinLevel() < level_start);

    // set the solution field
    // Field sol("solution", 1);
    // grid.AddField(&sol);
    // init(&grid, &sol);
    // sol.bctype(M_BC_EXTRAP);

    // grid.GhostPull(&sol);
    grid.GhostPull(&scal);

    // measure the error
    real_t normi;
    Error  error;
    error.Normi(&grid, &scal, &lambda_error, &normi);
    m_log("epsilon = %e, error = %.12e", epsilon, normi);

    // measure the moments
    real_t moment0, moment1[3];
    moment(&grid, &scal, &moment0, moment1);
    m_log("analytical moments after refinement: %e vs %e -> error = %.12e", fine_moment0, moment0, abs(fine_moment0 - moment0));
    real_t mom_smooth_error[4];
    mom_smooth_error[0] = fabs(fine_moment0 - moment0);
    for (lda_t ida = 0; ida < 3; ++ida) {
        mom_smooth_error[ida + 1] = fabs(fine_moment1[ida] - moment1[ida]);
    }

    ASSERT_LT(normi, 3.0 * epsilon);
    if (grid.interp()->Nt() > 0) {
        // for (lda_t ida = 0; ida < 4; ++ida) {
        //     ASSERT_LT(mom_smooth_error[ida], zero_tol);
        // }
        ASSERT_LT(mom_smooth_error[0], zero_tol);
        ASSERT_LT(mom_smooth_error[1], zero_tol);
        ASSERT_LT(mom_smooth_error[2], zero_tol);
        ASSERT_LT(mom_smooth_error[3], zero_tol);
    }

    // cleanup the fields
    // m_log("free fields");
    // grid.DeleteField(&sol);
    m_log("free fields");
    grid.DeleteField(&scal);
};

/**
 * @brief tests the epsilon condition, including at the boundary conditions (no moment check here!)
 * 
 */
TEST_P(Epsilon, extrap) {
    //-------------------------------------------------------------------------
    // get the parameters
    real_t epsilon = GetParam();
    level_t level_start = 4;

    // let's go
    bool  period[3]   = {false, false, false};
    lid_t grid_len[3] = {1, 1, 1};
    Grid  grid(level_start, period, grid_len, MPI_COMM_WORLD, nullptr);
    grid.level_limit(0, level_start);


    Field scal("scalar", 1);
    grid.AddField(&scal);
    scal.bctype(M_BC_EXTRAP);

    // InitCond_Epsilon init;
    // init(&grid, &scal);
    lambda_setvalue_t lambda_initcond = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        real_t* data = block->data(fid).Write(i0, i1, i2);
        // set value
        data[0] = scalar_exp(pos, center, sigma);
    };
    SetValue init(lambda_initcond);
    init(&grid, &scal);

    grid.SetTol(epsilon * 1e+20, epsilon);

    // get the ghosting
    grid.GhostPull(&scal);

    // // register the moments
    // BMoment moment;
    // real_t  fine_moment0, fine_moment1[3];
    // moment(&grid, &scal, &fine_moment0, fine_moment1);

    // coarsen
    short_t count     = 1;
    level_t min_level = level_start;
    do {
        min_level = grid.MinLevel();
        grid.Coarsen(&scal);

        level_t tmp_min_lvl = grid.MinLevel();
        level_t tmp_max_lvl = grid.MaxLevel();
        m_log("Coarsening: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

        grid.GhostPull(&scal);
        // real_t coarse_moment0, coarse_moment1[3];
        // moment(&grid, &scal, &coarse_moment0, coarse_moment1);
        // m_log("moments after coarsening: %e vs %e -> error = %e", fine_moment0, coarse_moment0, abs(fine_moment0 - coarse_moment0));

    } while (grid.MinLevel() < min_level && grid.MinLevel() > 0);

    // measure the moments
    grid.GhostPull(&scal);
    // real_t coarse_moment0, coarse_moment1[3];
    // moment(&grid, &scal, &coarse_moment0, coarse_moment1);
    // m_log("moments after coarsening: %e vs %e -> error = %e", fine_moment0, coarse_moment0, abs(fine_moment0 - coarse_moment0));

    // and go back up
    std::list<Patch> patch;
    real_t           origin[3] = {0.0, 0.0, 0.0};
    real_t           length[3] = {(real_t)grid_len[0], (real_t)grid_len[1], (real_t)grid_len[2]};
    patch.push_back(Patch(origin, length, level_start));
    do {
        min_level = grid.MinLevel();
        // force the field refinement using a patch
        grid.GhostPull(&scal);
        grid.AdaptMagic(nullptr, &patch, nullptr, &cback_StatusCheck, nullptr, &cback_UpdateDependency, nullptr);

        level_t tmp_min_lvl = grid.MinLevel();
        level_t tmp_max_lvl = grid.MaxLevel();
        m_log("Refinement: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

    } while (grid.MinLevel() < level_start);

    // set the solution field
    // Field sol("solution", 1);
    // grid.AddField(&sol);
    // init(&grid, &sol);
    // sol.bctype(M_BC_EXTRAP);

    // grid.GhostPull(&sol);
    grid.GhostPull(&scal);

    // measure the error
    real_t normi;
    Error  error;
    error.Normi(&grid, &scal, &lambda_error, &normi);
    m_log("epsilon = %e, error = %.12e", epsilon, normi);

    // measure the moments
    // real_t moment0, moment1[3];
    // moment(&grid, &scal, &moment0, moment1);
    // m_log("analytical moments after refinement: %e vs %e -> error = %.12e", fine_moment0, moment0, abs(fine_moment0 - moment0));
    // real_t mom_smooth_error[4];
    // mom_smooth_error[0] = fabs(fine_moment0 - moment0);
    // for (lda_t ida = 0; ida < 3; ++ida) {
    //     mom_smooth_error[ida + 1] = fabs(fine_moment1[ida] - moment1[ida]);
    // }

    ASSERT_LT(normi, 3.0 * epsilon);
    // if (grid.interp()->Nt() > 0) {
    //     for (lda_t ida = 0; ida < 4; ++ida) {
    //         ASSERT_LT(mom_smooth_error[ida], zero_tol);
    //     }
    // }

    // cleanup the fields
    // m_log("free fields");
    // grid.DeleteField(&sol);
    // m_log("free fields");
    grid.DeleteField(&scal);
}

INSTANTIATE_TEST_SUITE_P(ValidWavelet,
                         Epsilon,
                         testing::Values(1e-2, 1e-4, 1e-6));

// INSTANTIATE_TEST_SUITE_P(ValidWaveletExtrap,
//                          EpsilonExtrap,
//                          testing::Values(1e-2, 1e-4, 1e-6));