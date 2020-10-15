
#include <mpi.h>

#include "doop.hpp"
#include "error.hpp"
#include "grid.hpp"
#include "gtest/gtest.h"
#include "ioh5.hpp"
#include "murphy.hpp"
#include "setvalues.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define ORDER_TOL 0.3
#define BLVL 2

using std::list;

class valid_Wavelet_Convergence : public ::testing::Test {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

//==============================================================================================================================
TEST_F(valid_Wavelet_Convergence, ghost_reconstruction_periodic_sin) {
    // init the errors
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        // setup the mesh
        bool  period[3] = {true, true, true};
        lid_t L[3]      = {1, 1, 1};
        L[id]           = 3;

        for (level_t il = 0; il < 1; ++il) {
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il + BLVL + 1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            // create the test file
            string fieldName = "sinus" + std::to_string(id) + "__" + std::to_string(il);
            Field  test(fieldName, 1);
            grid.AddField(&test);

            // put a sinus
            const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t freq[3]    = {3.0, 1.0, 2.0};
            const real_t alpha[3]   = {1.0, 1.0, 1.0};
            SetSinus     field_init(sin_len, freq, alpha);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            Field sol("sol", 1);
            grid.AddField(&sol);
            SetSinus field_sol(sin_len, freq, alpha, grid.interp());
            field_sol(&grid, &sol);

            // now, we need to check
            ErrorCalculator error(grid.interp());
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDER_TOL);
    }
}
//==============================================================================================================================
TEST_F(valid_Wavelet_Convergence, ghost_reconstruction_periodic_cos) {
    // init the errors
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        // setup the mesh
        bool  period[3] = {true, true, true};
        lid_t L[3]      = {1, 1, 1};
        L[id]           = 3;

        for (level_t il = 0; il < 2; ++il) {
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il + BLVL + 1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            // create the test file
            string fieldName = "cosinus" + std::to_string(id) + "__" + std::to_string(il);
            Field  test(fieldName, 3);
            grid.AddField(&test);

            // put a sinus
            const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t freq[3]    = {2.0, 3.0, 1.0};
            const real_t alpha[3]   = {1.0, 1.0, 1.0};
            SetCosinus   field_init(sin_len, freq, alpha);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            Field sol("sol", 3);
            grid.AddField(&sol);
            SetCosinus field_sol(sin_len, freq, alpha, grid.interp());
            field_sol(&grid, &sol);

            // now, we need to check
            ErrorCalculator error(grid.interp());
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDER_TOL);
    }
}

//==============================================================================================================================
TEST_F(valid_Wavelet_Convergence, ghost_reconstruction_extrap_cos) {
    // init the errors
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        for (level_t il = 0; il < 2; ++il) {
            bool  period[3] = {false, false, false};
            lid_t L[3]      = {1, 1, 1};
            L[id]           = 3;
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1.0;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il + BLVL + 1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            Field test("test", 1);
            grid.AddField(&test);
            test.bctype(M_BC_EXTRAP);

            const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t freq[3]    = {2.0, 2.0, 2.0};
            const real_t alpha[3]   = {1.0, 1.0, 1.0};
            SetCosinus     field_init(sin_len, freq, alpha);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            Field sol("sol", 1);
            grid.AddField(&sol);
            SetCosinus field_sol(sin_len, freq, alpha, grid.interp());
            field_sol(&grid, &sol);

            // mask both the sol and the result
            // MaskPhysBC mask(L);
            // mask(&grid, &test);
            // mask(&grid, &sol);

            // now, we need to check
            real_t          norm2, normi;
            ErrorCalculator error(grid.interp());
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDER_TOL);
    }
}

//==============================================================================================================================
TEST_F(valid_Wavelet_Convergence, ghost_reconstruction_perper_dirichlet0_polynom) {
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        for (level_t il = 0; il < 2; il++) {
            bool period[3] = {true, true, true};
            period[id]     = false;
            lid_t L[3]     = {1, 1, 1};
            L[id]          = 3;
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1.0;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il + BLVL + 1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            Field test("test", 1);
            grid.AddField(&test);
            // test.bctype(M_BC_EXTRAP);
            // put the DIRICHLET BC in the direction of interest
            test.bctype(M_BC_DIR, 0, 2 * id);
            test.bctype(M_BC_DIR, 0, 2 * id + 1);

            // create the initial field
            lid_t  deg[3]     = {0.0, 0.0, 0.0};
            real_t dir[3]     = {0.0, 0.0, 0.0};
            real_t shift[3]   = {0.0, 0.0, 0.0};
            deg[id]           = M_WAVELET_N + 2;
            dir[id]           = -1.0;
            dir[(id + 1) % 3] = pow(L[id] / 2.0, deg[id]);
            shift[id]         = L[id] / 2.0;

            SetPolynom field_init(deg, dir, shift);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            Field sol("sol", 1);
            grid.AddField(&sol);
            SetPolynom field_sol(deg, dir, shift, grid.interp());
            field_sol(&grid, &sol);

            // mask both the sol and the result
            // MaskPhysBC mask(L);
            // mask(&grid, &test);
            // mask(&grid, &sol);

            // now, we need to check
            real_t          norm2, normi;
            ErrorCalculator error(grid.interp());
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDER_TOL);
    }
}

//==============================================================================================================================
TEST_F(valid_Wavelet_Convergence, ghost_reconstruction_perper_neuman0_cos) {
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        for (level_t il = 0; il < 2; il++) {
            bool period[3] = {true, true, true};
            period[id]     = false;
            lid_t L[3]     = {1, 1, 1};
            L[id]          = 3;
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1.0;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il + BLVL + 1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            Field test("test", 1);
            grid.AddField(&test);
            // test.bctype(M_BC_EXTRAP);
            // put the DIRICHLET BC in the direction of interest
            test.bctype(M_BC_NEU, 0, 2 * id);
            test.bctype(M_BC_NEU, 0, 2 * id + 1);

            // create the initial field
            real_t freq[3]   = {0.0, 0.0, 0.0};
            real_t coslen[3] = {L[0], L[1], L[2]};
            real_t alpha[3]  = {0.0, 0.0, 0.0};
            freq[id]         = 2.0;
            alpha[id]        = 1.0;

            SetCosinus field_init(coslen, freq, alpha);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            Field sol("sol", 1);
            grid.AddField(&sol);
            SetCosinus field_sol(coslen, freq, alpha, grid.interp());
            field_sol(&grid, &sol);

            // mask both the sol and the result
            // MaskPhysBC mask(L);
            // mask(&grid, &test);
            // mask(&grid, &sol);

            // now, we need to check
            real_t          norm2, normi;
            ErrorCalculator error(grid.interp());
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDER_TOL);
    }
}
