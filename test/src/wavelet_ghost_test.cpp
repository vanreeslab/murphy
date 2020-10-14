
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

#define DOUBLE_TOL 1e-08 // we have to pay attention that the error get's worse with the order -> we extrapolate further and further

using std::list;

class valid_Wavelet_Ghost : public ::testing::Test {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

TEST_F(valid_Wavelet_Ghost, ghost_exact_extrap) {
    for (lda_t id = 0; id < 3; id++) {
        bool  period[3] = {false, false, false};
        lid_t L[3]      = {1, 1, 1};
        L[id]           = 3;
        Grid grid(0, period, L, MPI_COMM_WORLD, nullptr);

        // create the patch refinement to refine the middle tree
        real_t origin[3]      = {0.0, 0.0, 0.0};
        origin[id]            = 1.0;
        real_t      length[3] = {1.0, 1.0, 1.0};
        Patch       p1(origin, length, 1);
        list<Patch> patch{p1};
        grid.Adapt(&patch);

        Field test("test", 1);
        grid.AddField(&test);
        test.bctype(M_BC_EXTRAP);

        // create the initial field
        const lid_t  deg[3] = {M_WAVELET_N-1, M_WAVELET_N-1, M_WAVELET_N-1};
        const real_t dir[3] = {M_SQRT2, -M_PI, M_E};
        SetPolynom   field_init(deg, dir);
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
        SetPolynom field_sol(deg, dir, grid.interp());
        field_sol(&grid, &sol);

        // mask both the sol and the result
        // MaskPhysBC mask(L);
        // mask(&grid, &test);
        // mask(&grid, &sol);

        // now, we need to check
        real_t          norm2, normi;
        ErrorCalculator error(grid.interp());
        error.Norms(&grid, &test, &sol, &norm2, &normi);

        m_log("checking in dim %d: the two norms: %e %e", id, normi, norm2);
        ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
        ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    }
}

TEST_F(valid_Wavelet_Ghost, ghost_convergence_extrap) {
    real_t ei[2];
    real_t e2[2];
    for (lda_t id = 0; id < 3; id++) {
        for (level_t il = 1; il < 3; il++) {
            bool  period[3] = {false, false, false};
            lid_t L[3]      = {1, 1, 1};
            L[id]           = 3;
            Grid grid(il, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1.0;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il + 1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            Field test("test", 1);
            grid.AddField(&test);
            test.bctype(M_BC_EXTRAP);

            // create the initial field
            const lid_t  deg[3] = {M_WAVELET_N + 2, M_WAVELET_N + 1, M_WAVELET_N};
            const real_t dir[3] = {M_SQRT2, -M_PI, M_E};
            SetPolynom   field_init(deg, dir);
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
            SetPolynom field_sol(deg, dir, grid.interp());
            field_sol(&grid, &sol);

            // mask both the sol and the result
            // MaskPhysBC mask(L);
            // mask(&grid, &test);
            // mask(&grid, &sol);

            // now, we need to check
            real_t          norm2, normi;
            ErrorCalculator error(grid.interp());
            error.Norms(&grid, &test, &sol, e2 + il, ei + il);

            m_log("checking in dim %d, level = %d: the two norms: %e %e", id, il, ei[il], e2[il]);
        }

        real_t conv2 = -log(e2[2] / e2[1]) / log(2);
        real_t convi = -log(ei[2] / ei[1]) / log(2);
        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);

        ASSERT_GE(convi, M_WAVELET_N - 0.1);
        ASSERT_GE(conv2, M_WAVELET_N - 0.1);
    }
}
