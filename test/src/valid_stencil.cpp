#include <cmath>
#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "gtest/gtest.h"
#include "operator/advection.hpp"
#include "operator/error.hpp"
#include "operator/xblas.hpp"
#include "tools/ioh5.hpp"
#include "valid_toolbox.hpp"

#define DOUBLE_TOL 1e-12

class Adapt : public ::testing::TestWithParam<int> {
    void SetUp() override{};
    void TearDown() override{};
};

using std::list;
using std::string;



TEST_P(Adapt, weno_periodic_cosinus_expr) {
    int case_id = GetParam();
    m_log("--------------------------------------------------------------------------------");
    m_log("case id: %d", case_id);

    // init the errors
    real_t erri_adv_weno_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_weno_5[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_cons_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_cons_5[3] = {0.0, 0.0, 0.0};

    // setup the mesh
    bool  period[3] = {true, true, true};
    lid_t L[3]      = {2, 2, 2};

    // see if we run the tests
    // bool do_weno_3 = false;
    // bool do_weno_5 = false;

    real_t rand_vel[3];
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        rand_vel[0] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[1] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[2] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
    }
    MPI_Bcast(rand_vel, 3, M_MPI_REAL, 0, MPI_COMM_WORLD);
    m_log("velocity is %e %e %e", rand_vel[0], rand_vel[1], rand_vel[2]);

    for (level_t il = 0; il < 2; ++il) {
        level_t clevel = il;
        Grid    grid(clevel, period, L, MPI_COMM_WORLD, nullptr);
        // adapt
        {
            list<Patch> patch_list;
            TreeCaseIdToPatch(clevel, case_id, &patch_list);
            grid.Adapt(&patch_list);
        }

        // do_weno_3 = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
        // do_weno_5 = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;

        // create the needed fields
        string fieldName = "field" + std::to_string(il);
        // string solName   = "sol" + std::to_string(il);
        string velName  = "vel" + std::to_string(il);
        string diffName = "deriv" + std::to_string(il);
        Field  test(fieldName, 1);
        // Field  sol(solName, 1);
        Field vel(velName, 3);
        Field dtest(diffName, 1);
        grid.AddField(&test);
        // grid.AddField(&sol);
        // grid.AddField(&vel);
        grid.AddField(&dtest);

        // set a constant velocity
        {
            lambda_i3_t<real_t,lda_t> set_vel = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const lda_t ida) -> real_t {
                return rand_vel[ida];
            };

            vel.is_expr(true);
            grid.SetExpr(&vel,set_vel);
            // const bidx_t ghost_len[2] = {3, 3};
            // SetValue     field_init(pol_op, ghost_len);
            // field_init(&grid, &vel);
        }

        // set the test field
        {
            // // cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x)
            // const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            // const real_t freq[3]    = {2.0, 2.0, 2.0};
            // const real_t alpha[3]   = {1.0, 1.0, 1.0};
            // SetCosinus   field_init(sin_len, freq, alpha);
            // field_init(&grid, &test);
            lambda_setvalue_t sin_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid, 0)(i0, i1, i2) = sin(2.0 * M_PI * 2.0 / ((real_t)L[0]) * pos[0]) +
                                                  sin(2.0 * M_PI * 2.0 / ((real_t)L[1]) * pos[1]) +
                                                  sin(2.0 * M_PI * 2.0 / ((real_t)L[2]) * pos[2]);
            };
            SetValue field_init(sin_op);
            field_init(&grid, &test);
        }

        // set the solution
        // {
        // // -> the solution: u* df/dx + v * df/dy + w*df/dz
        // const real_t freq[3]        = {2.0, 2.0, 2.0};
        // const real_t sin_len[3]     = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
        // const real_t alpha_sol_0[3] = {rand_vel[0] * 2.0 * M_PI * freq[0] / L[0],
        //                                rand_vel[1] * 2.0 * M_PI * freq[1] / L[1],
        //                                rand_vel[2] * 2.0 * M_PI * freq[2] / L[2]};
        // SetSinus     sol_init(sin_len, freq, alpha_sol_0, grid.interp());
        // sol_init(&grid, &sol);
        lambda_error_t cos_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            // call the function
            return (-rand_vel[0]) * (4.0 * M_PI / ((real_t)L[0])) * cos(2.0 * M_PI * 2.0 / ((real_t)L[0]) * pos[0]) +
                   (-rand_vel[1]) * (4.0 * M_PI / ((real_t)L[1])) * cos(2.0 * M_PI * 2.0 / ((real_t)L[1]) * pos[1]) +
                   (-rand_vel[2]) * (4.0 * M_PI / ((real_t)L[2])) * cos(2.0 * M_PI * 2.0 / ((real_t)L[2]) * pos[2]);
        };
        //     SetValue field_init(cos_op, grid.interp());
        //     field_init(&grid, &test);
        // }

        // if (do_weno_3)
        {
            Advection<M_WENO_Z, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_3 + il);
            m_log("M_WENO_Z - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_3[il]);
        }
        // if (do_weno_5)
        {
            Advection<M_WENO_Z, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_5 + il);
            m_log("M_WENO_Z - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_5[il]);
        }
        {
            Advection<M_CONS, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_cons_3 + il);
            m_log("M_CONS - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_cons_3[il]);
        }
        // if (do_weno_5)
        {
            Advection<M_CONS, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_cons_5 + il);
            m_log("M_CONS - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_cons_5[il]);
        }
    }
    // if (do_weno_3)
    {
        real_t convi = -log(erri_adv_weno_3[1] / erri_adv_weno_3[0]) / log(2);
        m_log("M_ADV_WENO - 3: the convergence orders are: norm_i:%e -> min = 1, ideal = 3", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 1) - 0.75);
    }
    // if (do_weno_5)
    {
        real_t convi = -log(erri_adv_weno_5[1] / erri_adv_weno_5[0]) / log(2);
        m_log("M_ADV_WENO - 5: the convergence orders are: norm_i:%e -> min = 3, ideal = 5", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.75);
    }
    {
        real_t convi = -log(erri_adv_cons_3[1] / erri_adv_cons_3[0]) / log(2);
        m_log("M_ADV_CONS - 3: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.75);
    }
    // if (do_weno_5)
    {
        real_t convi = -log(erri_adv_cons_5[1] / erri_adv_cons_5[0]) / log(2);
        m_log("M_ADV_CONS - 5: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 5) - 0.75);
    }
}


TEST_P(Adapt, weno_periodic_cosinus_field) {
    int case_id = GetParam();
    m_log("--------------------------------------------------------------------------------");
    m_log("case id: %d", case_id);

    // init the errors
    real_t erri_adv_weno_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_weno_5[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_cons_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_cons_5[3] = {0.0, 0.0, 0.0};

    // setup the mesh
    bool  period[3] = {true, true, true};
    lid_t L[3]      = {2, 2, 2};

    // see if we run the tests
    // bool do_weno_3 = false;
    // bool do_weno_5 = false;

    real_t rand_vel[3];
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        rand_vel[0] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[1] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[2] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
    }
    MPI_Bcast(rand_vel, 3, M_MPI_REAL, 0, MPI_COMM_WORLD);
    m_log("velocity is %e %e %e", rand_vel[0], rand_vel[1], rand_vel[2]);

    for (level_t il = 0; il < 2; ++il) {
        level_t clevel = il;
        Grid    grid(clevel, period, L, MPI_COMM_WORLD, nullptr);
        // adapt
        {
            list<Patch> patch_list;
            TreeCaseIdToPatch(clevel, case_id, &patch_list);
            grid.Adapt(&patch_list);
        }

        // do_weno_3 = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
        // do_weno_5 = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;

        // create the needed fields
        string fieldName = "field" + std::to_string(il);
        // string solName   = "sol" + std::to_string(il);
        string velName  = "vel" + std::to_string(il);
        string diffName = "deriv" + std::to_string(il);
        Field  test(fieldName, 1);
        // Field  sol(solName, 1);
        Field vel(velName, 3);
        Field dtest(diffName, 1);
        grid.AddField(&test);
        // grid.AddField(&sol);
        grid.AddField(&vel);
        grid.AddField(&dtest);

        // set a constant velocity
        {
            // const lid_t  deg[3]   = {0, 0, 0};
            // const real_t shift[3] = {0.0, 0.0, 0.0};
            // for (lda_t ida = 0; ida < 3; ++ida) {
            //     const real_t dir[3] = {rand_vel[ida], 0.0, 0.0};
            //     SetPolynom   vel_init(deg, dir, shift, grid.interp());
            //     vel_init(&grid, &vel, ida);  // put 1.0 in the indicated direction only
            // }
            // grid.GhostPull(&vel);
            lambda_setvalue_t pol_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid, 0)(i0, i1, i2) = rand_vel[0];
                block->data(fid, 1)(i0, i1, i2) = rand_vel[1];
                block->data(fid, 2)(i0, i1, i2) = rand_vel[2];
            };
            const bidx_t ghost_len[2] = {3, 3};
            SetValue     field_init(pol_op, ghost_len);
            field_init(&grid, &vel);
        }

        // set the test field
        {
            // // cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x)
            // const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            // const real_t freq[3]    = {2.0, 2.0, 2.0};
            // const real_t alpha[3]   = {1.0, 1.0, 1.0};
            // SetCosinus   field_init(sin_len, freq, alpha);
            // field_init(&grid, &test);
            lambda_setvalue_t sin_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid, 0)(i0, i1, i2) = sin(2.0 * M_PI * 2.0 / ((real_t)L[0]) * pos[0]) +
                                                  sin(2.0 * M_PI * 2.0 / ((real_t)L[1]) * pos[1]) +
                                                  sin(2.0 * M_PI * 2.0 / ((real_t)L[2]) * pos[2]);
            };
            SetValue field_init(sin_op);
            field_init(&grid, &test);
        }

        // set the solution
        // {
        // // -> the solution: u* df/dx + v * df/dy + w*df/dz
        // const real_t freq[3]        = {2.0, 2.0, 2.0};
        // const real_t sin_len[3]     = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
        // const real_t alpha_sol_0[3] = {rand_vel[0] * 2.0 * M_PI * freq[0] / L[0],
        //                                rand_vel[1] * 2.0 * M_PI * freq[1] / L[1],
        //                                rand_vel[2] * 2.0 * M_PI * freq[2] / L[2]};
        // SetSinus     sol_init(sin_len, freq, alpha_sol_0, grid.interp());
        // sol_init(&grid, &sol);
        lambda_error_t cos_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            // call the function
            return (-rand_vel[0]) * (4.0 * M_PI / ((real_t)L[0])) * cos(2.0 * M_PI * 2.0 / ((real_t)L[0]) * pos[0]) +
                   (-rand_vel[1]) * (4.0 * M_PI / ((real_t)L[1])) * cos(2.0 * M_PI * 2.0 / ((real_t)L[1]) * pos[1]) +
                   (-rand_vel[2]) * (4.0 * M_PI / ((real_t)L[2])) * cos(2.0 * M_PI * 2.0 / ((real_t)L[2]) * pos[2]);
        };
        //     SetValue field_init(cos_op, grid.interp());
        //     field_init(&grid, &test);
        // }

        // if (do_weno_3)
        {
            Advection<M_WENO_Z, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_3 + il);
            m_log("M_WENO_Z - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_3[il]);
        }
        // if (do_weno_5)
        {
            Advection<M_WENO_Z, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_5 + il);
            m_log("M_WENO_Z - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_5[il]);
        }
        {
            Advection<M_CONS, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_cons_3 + il);
            m_log("M_CONS - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_cons_3[il]);
        }
        // if (do_weno_5)
        {
            Advection<M_CONS, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_cons_5 + il);
            m_log("M_CONS - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_cons_5[il]);
        }
    }
    // if (do_weno_3)
    {
        real_t convi = -log(erri_adv_weno_3[1] / erri_adv_weno_3[0]) / log(2);
        m_log("M_ADV_WENO - 3: the convergence orders are: norm_i:%e -> min = 1, ideal = 3", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 1) - 0.75);
    }
    // if (do_weno_5)
    {
        real_t convi = -log(erri_adv_weno_5[1] / erri_adv_weno_5[0]) / log(2);
        m_log("M_ADV_WENO - 5: the convergence orders are: norm_i:%e -> min = 3, ideal = 5", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.75);
    }
    {
        real_t convi = -log(erri_adv_cons_3[1] / erri_adv_cons_3[0]) / log(2);
        m_log("M_ADV_CONS - 3: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.75);
    }
    // if (do_weno_5)
    {
        real_t convi = -log(erri_adv_cons_5[1] / erri_adv_cons_5[0]) / log(2);
        m_log("M_ADV_CONS - 5: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 5) - 0.75);
    }
}

INSTANTIATE_TEST_SUITE_P(ValidStencil,
                         Adapt,
                         testing::Range(0, 256));

class ValidStencilUniform : public ::testing::Test {
    void SetUp() override{};
    void TearDown() override{};
};

TEST_F(ValidStencilUniform, weno_extrap_cosinus) {
    // init the errors
    real_t erri_adv_weno_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_weno_5[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_cons_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_cons_5[3] = {0.0, 0.0, 0.0};

    // setup the mesh
    // bool  period[3] = {false, false, false};
    bool  period[3] = {true, true, true};
    lid_t L[3]      = {2, 2, 2};

    // see if we run the tests
    // bool do_weno_3 = false;
    // bool do_weno_5 = false;

    real_t rand_vel[3];
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        rand_vel[0] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[1] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[2] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
    }
    MPI_Bcast(rand_vel, 3, M_MPI_REAL, 0, MPI_COMM_WORLD);
    m_log("velocity is %e %e %e", rand_vel[0], rand_vel[1], rand_vel[2]);

    for (level_t il = 0; il < 2; ++il) {
        Grid grid(il, period, L, MPI_COMM_WORLD, nullptr);

        // do_weno_3 = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
        // do_weno_5 = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;

        // create the test file
        string fieldName = "field" + std::to_string(il);
        // string solName   = "sol" + std::to_string(il);
        string velName  = "vel" + std::to_string(il);
        string diffName = "deriv" + std::to_string(il);
        Field  test(fieldName, 1);
        // Field  sol(solName, 1);
        Field vel(velName, 3);
        Field dtest(diffName, 1);
        grid.AddField(&test);
        // grid.AddField(&sol);
        grid.AddField(&vel);
        grid.AddField(&dtest);

        test.bctype(M_BC_EXTRAP);

        // set a constant velocity
        {
            // const lid_t  deg[3]   = {0, 0, 0};
            // const real_t shift[3] = {0.0, 0.0, 0.0};
            // for (lda_t ida = 0; ida < 3; ++ida) {
            //     const real_t dir[3] = {rand_vel[ida], 0.0, 0.0};
            //     SetPolynom   vel_init(deg, dir, shift, grid.interp());
            //     vel_init(&grid, &vel, ida);  // put 1.0 in the indicated direction only
            // }
            // grid.GhostPull(&vel);
            lambda_setvalue_t pol_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid, 0)(i0, i1, i2) = rand_vel[0];
                block->data(fid, 1)(i0, i1, i2) = rand_vel[1];
                block->data(fid, 2)(i0, i1, i2) = rand_vel[2];
            };
            const bidx_t ghost_len[2] = {3, 3};
            SetValue     field_init(pol_op, ghost_len);
            field_init(&grid, &vel);
        }

        // set the test field
        {
            // // cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x)
            // const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            // const real_t freq[3]    = {2.0, 2.0, 2.0};
            // const real_t alpha[3]   = {1.0, 1.0, 1.0};
            // SetCosinus   field_init(sin_len, freq, alpha);
            // field_init(&grid, &test);
            lambda_setvalue_t sin_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid,0)(i0, i1, i2) = sin(4.0 * M_PI / ((real_t)L[0]) * pos[0]) +
                                               sin(4.0 * M_PI / ((real_t)L[1]) * pos[1]) +
                                               sin(4.0 * M_PI / ((real_t)L[2]) * pos[2]);
            };
            SetValue field_init(sin_op);
            field_init(&grid, &test);
        }

        // set the solution
        // {
        // // -> the solution: u* df/dx + v * df/dy + w*df/dz
        // const real_t freq[3]        = {2.0, 2.0, 2.0};
        // const real_t sin_len[3]     = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
        // const real_t alpha_sol_0[3] = {rand_vel[0] * 2.0 * M_PI * freq[0] / L[0],
        //                                rand_vel[1] * 2.0 * M_PI * freq[1] / L[1],
        //                                rand_vel[2] * 2.0 * M_PI * freq[2] / L[2]};
        // SetSinus     sol_init(sin_len, freq, alpha_sol_0, grid.interp());
        // sol_init(&grid, &sol);
        lambda_error_t cos_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            // call the function
            return rand_vel[0] * (-4.0 * M_PI / ((real_t)L[0])) * cos(4.0 * M_PI / ((real_t)L[0]) * pos[0]) +
                   rand_vel[1] * (-4.0 * M_PI / ((real_t)L[1])) * cos(4.0 * M_PI / ((real_t)L[1]) * pos[1]) +
                   rand_vel[2] * (-4.0 * M_PI / ((real_t)L[2])) * cos(4.0 * M_PI / ((real_t)L[2]) * pos[2]);
        };
        //     SetValue field_init(cos_op, grid.interp());
        //     field_init(&grid, &test);
        // }

        // if (do_weno_3)
        {
            Advection<M_WENO_Z, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_3 + il);
            m_log("M_WENO_Z - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_3[il]);

            // check the moment
            real_t sum;
            BAvg   sum_grid;
            sum_grid(&grid, &dtest, &sum);
            m_log("sum of dtest = %e", sum);

            ASSERT_LT(sum / (L[0] * L[1] * L[2]), DOUBLE_TOL);
        }
        // if (do_weno_5)
        {
            Advection<M_WENO_Z, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_5 + il);
            m_log("M_WENO_Z - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_5[il]);

            // check the moment
            real_t sum;
            BAvg   sum_grid;
            sum_grid(&grid, &dtest, &sum);
            m_log("sum of dtest = %e", sum);

            ASSERT_LT(sum / (L[0] * L[1] * L[2]), DOUBLE_TOL);
        }
        {
            Advection<M_CONS, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_cons_3 + il);
            m_log("M_CONS - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_cons_3[il]);
        }
        // if (do_weno_5)
        {
            Advection<M_CONS, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            Error error;
            error.Normi(&grid, &dtest, &cos_sol, erri_adv_cons_5 + il);
            m_log("M_CONS - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_cons_5[il]);
        }
    }
    // if (do_weno_3)
    {
        real_t convi = -log(erri_adv_weno_3[1] / erri_adv_weno_3[0]) / log(2);
        m_log("M_ADV_WENO - 3: the convergence orders are: norm_i:%e -> target: min = 1, ideal = 3", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 1) - 0.5);
    }
    // if (do_weno_5)
    {
        real_t convi = -log(erri_adv_weno_5[1] / erri_adv_weno_5[0]) / log(2);
        m_log("M_ADV_WENO - 5: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.5);
    }
    {
        real_t convi = -log(erri_adv_cons_3[1] / erri_adv_cons_3[0]) / log(2);
        m_log("M_ADV_WENO - 3: the convergence orders are: norm_i:%e -> target: min = 1, ideal = 3", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 1) - 0.5);
    }
    // if (do_weno_5)
    {
        real_t convi = -log(erri_adv_cons_5[1] / erri_adv_cons_5[0]) / log(2);
        m_log("M_ADV_WENO - 5: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.5);
    }
}

// TEST_F(ValidStencilUniform, weno_extrap_poly) {
//     // init the errors
//     real_t erri_adv_weno_3[3] = {0.0, 0.0, 0.0};
//     real_t erri_adv_weno_5[3] = {0.0, 0.0, 0.0};

//     // setup the mesh
//     bool  period[3] = {false, false, false};
//     // bool  period[3] = {true, true, true};
//     lid_t L[3]      = {2, 2, 2};

//     // see if we run the tests
//     bool do_weno_3 = false;
//     bool do_weno_5 = false;

//     const real_t rand_vel[3] = {1.0, 1.0, 1.0};  //{-1.0 + ((real_t)std::rand() / (real_t)RAND_MAX) * 2.0,
//                                                  //-1.0 + ((real_t)std::rand() / (real_t)RAND_MAX) * 2.0,
//                                                  //-1.0 + ((real_t)std::rand() / (real_t)RAND_MAX) * 2.0};
//     m_log("velocity is %e %e %e", rand_vel[0], rand_vel[1], rand_vel[2]);

//     for (level_t il = 0; il < 2; ++il) {
//         Grid grid(il, period, L, MPI_COMM_WORLD, nullptr);

//         do_weno_3 = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
//         do_weno_5 = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;

//         // create the test file
//         string fieldName = "field" + std::to_string(il);
//         // string solName   = "sol" + std::to_string(il);
//         string velName  = "vel" + std::to_string(il);
//         string diffName = "deriv" + std::to_string(il);
//         Field  test(fieldName, 1);
//         // Field  sol(solName, 1);
//         Field vel(velName, 3);
//         Field dtest(diffName, 1);
//         grid.AddField(&test);
//         // grid.AddField(&sol);
//         grid.AddField(&vel);
//         grid.AddField(&dtest);

//         test.bctype(M_BC_EXTRAP);

//         // set a constant velocity
//         {
//             // const lid_t  deg[3]   = {0, 0, 0};
//             // const real_t shift[3] = {0.0, 0.0, 0.0};
//             // for (lda_t ida = 0; ida < 3; ++ida) {
//             //     const real_t dir[3] = {rand_vel[ida], 0.0, 0.0};
//             //     SetPolynom   vel_init(deg, dir, shift, grid.interp());
//             //     vel_init(&grid, &vel, ida);  // put 1.0 in the indicated direction only
//             // }
//             // grid.GhostPull(&vel);
//             lambda_setvalue_t pol_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
//                 // get the position
//                 real_t pos[3];
//                 block->pos(i0, i1, i2, pos);

//                 // call the function
//                 block->data(fid, 0).Write(i0, i1, i2)[0] = rand_vel[0];
//                 block->data(fid, 1).Write(i0, i1, i2)[0] = rand_vel[1];
//                 block->data(fid, 2).Write(i0, i1, i2)[0] = rand_vel[2];
//             };
//             SetValue field_init(pol_op,grid.interp());
//             field_init(&grid, &vel);
//         }

//         // set the test field
//         {
//             // // cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x)
//             // const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
//             // const real_t freq[3]    = {2.0, 2.0, 2.0};
//             // const real_t alpha[3]   = {1.0, 1.0, 1.0};
//             // SetCosinus   field_init(sin_len, freq, alpha);
//             // field_init(&grid, &test);
//             lambda_setvalue_t sin_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
//                 // get the position
//                 real_t pos[3];
//                 block->pos(i0, i1, i2, pos);

//                 // call the function
//                 block->data(fid).Write(i0, i1, i2)[0] = pow(pos[0], 1)+ // m_max(M_WAVELET_N - 1, 0)) +
//                                                         pow(pos[1], 1)+ // m_max(M_WAVELET_N - 1, 0)) +
//                                                         pow(pos[2], 1); // m_max(M_WAVELET_N - 1, 0));
//             };
//             SetValue field_init(sin_op);
//             field_init(&grid, &test);
//         }

//         // set the solution
//         // {
//         // // -> the solution: u* df/dx + v * df/dy + w*df/dz
//         // const real_t freq[3]        = {2.0, 2.0, 2.0};
//         // const real_t sin_len[3]     = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
//         // const real_t alpha_sol_0[3] = {rand_vel[0] * 2.0 * M_PI * freq[0] / L[0],
//         //                                rand_vel[1] * 2.0 * M_PI * freq[1] / L[1],
//         //                                rand_vel[2] * 2.0 * M_PI * freq[2] / L[2]};
//         // SetSinus     sol_init(sin_len, freq, alpha_sol_0, grid.interp());
//         // sol_init(&grid, &sol);
//         lambda_error_t cos_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
//             // get the position
//             real_t pos[3];
//             block->pos(i0, i1, i2, pos);

//             // call the function
//             // return (-rand_vel[0]) * m_max(M_WAVELET_N - 1, 0) * pow(pos[0], m_max(M_WAVELET_N - 2, 0)) +
//             //    (-rand_vel[1]) * m_max(M_WAVELET_N - 1, 0) * pow(pos[1], m_max(M_WAVELET_N - 2, 0)) +
//             //    (-rand_vel[2]) * m_max(M_WAVELET_N - 1, 0) * pow(pos[2], m_max(M_WAVELET_N - 2, 0));

//             // const real_t sol =
//             return (-rand_vel[0]) +
//                    (-rand_vel[1]) +
//                    (-rand_vel[2]);
//         };
//         //     SetValue field_init(cos_op, grid.interp());
//         //     field_init(&grid, &test);
//         // }

//         if (do_weno_3) {
//             Advection<M_WENO_Z, 3> adv(&vel);
//             adv(&grid, &test, &dtest);
//             // now, we need to check
//             Error error;
//             error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_3 + il);
//             m_log("M_WENO_Z - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_3[il]);

//             // check the moment
//             real_t sum;
//             BAvg  sum_grid;
//             sum_grid(&grid, &dtest, &sum);
//             m_log("sum of dtest = %e", sum);

//             ASSERT_LT(sum / (L[0] * L[1] * L[2]), DOUBLE_TOL);
//         }
//         if (do_weno_5) {
//             Advection<M_WENO_Z, 5> adv(&vel);
//             adv(&grid, &test, &dtest);
//             // now, we need to check
//             Error error;
//             error.Normi(&grid, &dtest, &cos_sol, erri_adv_weno_5 + il);
//             m_log("M_WENO_Z - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_5[il]);

//             // check the moment
//             real_t sum;
//             BAvg  sum_grid;
//             sum_grid(&grid, &dtest, &sum);
//             m_log("sum of dtest = %e", sum);

//             ASSERT_LT(sum / (L[0] * L[1] * L[2]), DOUBLE_TOL);
//         }
//     }
//     if (do_weno_3) {
//         real_t convi = -log(erri_adv_weno_3[1] / erri_adv_weno_3[0]) / log(2);
//         m_log("M_ADV_WENO - 3: the convergence orders are: norm_i:%e -> target: min = 1, ideal = 3", convi);
//         ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 1) - 0.5);
//     }
//     if (do_weno_5) {
//         real_t convi = -log(erri_adv_weno_5[1] / erri_adv_weno_5[0]) / log(2);
//         m_log("M_ADV_WENO - 5: the convergence orders are: norm_i:%e", convi);
//         ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.5);
//     }
// }
