
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

#define DOUBLE_TOL 1e-11

using std::list;

class test_Wavelet_Convergence_Periodic : public ::testing::Test {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

//==============================================================================================================================
TEST_F(test_Wavelet_Convergence_Periodic, convergence_sinus) {
    for (lda_t id = 0; id < 3; id++) {
        // init the errors
        real_t erri[3] = {0.0, 0.0, 0.0};
        real_t err2[3] = {0.0, 0.0, 0.0};

        // setup the mesh
        bool  period[3] = {true, true, true};
        lid_t L[3]      = {1, 1, 1};
        L[id]           = 3;

        for (level_t il = 1; il < 3; ++il) {
            Grid grid(il, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin[3]      = {0.0, 0.0, 0.0};
            origin[id]            = 1;
            real_t      length[3] = {1.0, 1.0, 1.0};
            Patch       p1(origin, length, il+1);
            list<Patch> patch{p1};
            grid.Adapt(&patch);

            // create the test file
            string fieldName = "sinus" + std::to_string(id) + "__" + std::to_string(il);
            Field  test(fieldName, 1);
            grid.AddField(&test);

            // put a sinus
            const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t freq[3]    = {3.0, 3.0, 3.0};
            SetSinus     field_init(sin_len, freq);
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
            SetSinus field_sol(sin_len, freq, &grid);
            field_sol(&grid, &sol);

            // now, we need to check
            ErrorCalculator error(&grid);
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[2] / err2[1]) / log(2);
        real_t convi = -log(erri[2] / erri[1]) / log(2);

        m_log("the convergence orders are: norm_2:%e norm_i:%e", convi, conv2);
        ASSERT_GE(conv2, M_WAVELET_N-0.1);
        ASSERT_GE(convi, M_WAVELET_N-0.1);
    }
}