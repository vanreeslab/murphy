#include <mpi.h>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "doop.hpp"
#include "error.hpp"
#include "grid.hpp"
#include "gtest/gtest.h"
#include "ioh5.hpp"
#include "setvalues.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-13
#define ORDER2_TOL 1.0
#define ORDERI_TOL 1.2
#define BLVL 1

using std::list;
using std::string;

class ValidWaveletInterpolation : public ::testing::Test {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

//==============================================================================================================================
TEST_F(ValidWaveletInterpolation, ghost_reconstruction_periodic_sin) {
    // init the errors
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        // setup the mesh
        bool  period[3] = {true, true, true};
        lid_t L[3]      = {3, 3, 3};

        for (level_t il = 0; il < 2; ++il) {
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin11[3] = {0.0, 1.0, 1.0};
            real_t length11[3] = {3.0, 1.0, 1.0};
            Patch  p11(origin11, length11, il + BLVL + 1);
            real_t origin12[3] = {1.0, 0.0, 1.0};
            real_t length12[3] = {1.0, 3.0, 1.0};
            Patch  p12(origin12, length12, il + BLVL + 1);
            real_t origin2[3] = {0.0, 0.0, 0.0};
            real_t length2[3] = {1.0, 1.0, 1.0};
            Patch  p2(origin2, length2, il + BLVL);

            // list<Patch> patch{p11,p2};
            list<Patch> patch{p11, p12, p2};
            grid.Adapt(&patch);

            // create the test file
            string fieldName = "sinus" + std::to_string(id) + "__" + std::to_string(il);
            Field  test(fieldName, 1);
            grid.AddField(&test);

            // put a sinus
            const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t freq[3]    = {3.0, 3.0, 3.0};
            const real_t alpha[3]   = {1.0, 1.0, 1.0};
            SetSinus     field_init(sin_len, freq, alpha, grid.interp());
            // SetSinus field_init(sin_len, freq, alpha);
            field_init(&grid, &test);

            // we have perfect ghost -> what is not set will be perfect
            // we force the ghost recomputation
            test.ghost_status(false);

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
            Error error(grid.interp());

            // sanity check
            error.Norms(&grid, il + 1 + BLVL, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);
            m_log("checking in dim %d on HIGH: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL + 1), erri[il], err2[il]);
            error.Norms(&grid, il + BLVL, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);
            m_log("checking in dim %d on LOW: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
// if NT==0, the lowest level has NO error
#if (M_WAVELET_NT == 0)
            ASSERT_NEAR(erri[il], 0.0, DOUBLE_TOL);
            ASSERT_NEAR(err2[il], 0.0, DOUBLE_TOL);
#endif

            // and now check for the every level
            error.Norms(&grid, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);
            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);

            grid.DeleteField(&test);
            grid.DeleteField(&sol);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
    }
}
//==============================================================================================================================
TEST_F(ValidWaveletInterpolation, ghost_reconstruction_periodic_cos) {
    // init the errors
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        // setup the mesh
        bool  period[3] = {true, true, true};
        lid_t L[3]      = {3, 3, 3};

        for (level_t il = 0; il < 2; ++il) {
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin1[3] = {1.0, 1.0, 1.0};
            real_t length1[3] = {1.0, 1.0, 1.0};
            Patch  p1(origin1, length1, il + BLVL + 1);
            real_t origin2[3] = {0.0, 0.0, 0.0};
            real_t length2[3] = {3.0, 3.0, 3.0};
            Patch  p2(origin2, length2, il + BLVL);

            list<Patch> patch{p1, p2};
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
            Error error(grid.interp());
            // sanity check
            error.Norms(&grid, il + 1 + BLVL, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);
            m_log("checking in dim %d on HIGH: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
#if (M_WAVELET_NT == 0)
            error.Norms(&grid, il + BLVL, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);
            m_log("checking in dim %d on LOW: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
            // if NT==0, the lowest level has NO error
            ASSERT_NEAR(erri[il], 0.0, DOUBLE_TOL);
            ASSERT_NEAR(err2[il], 0.0, DOUBLE_TOL);
#endif

            error.Norms(&grid, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
    }
}

// //==============================================================================================================================
// TEST_F(ValidWaveletInterpolation, ghost_reconstruction_extrap_cos) {
//     // init the errors
//     real_t erri[2];
//     real_t err2[2];
//     for (lda_t id = 0; id < 3; id++) {
//         for (level_t il = 0; il < 2; ++il) {
//             bool  period[3] = {false, false, false};
//             lid_t L[3]      = {3, 3, 3};
//             Grid  grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

//             // create the patch refinement to refine the middle tree
//             real_t origin1[3] = {1.0, 1.0, 1.0};
//             real_t length1[3] = {1.0, 1.0, 1.0};
//             Patch  p1(origin1, length1, il + BLVL + 1);
//             real_t origin2[3] = {0.0, 0.0, 0.0};
//             real_t length2[3] = {3.0, 3.0, 3.0};
//             Patch  p2(origin2, length2, il + BLVL);

//             list<Patch> patch{p1, p2};
//             grid.Adapt(&patch);

//             Field test("test", 1);
//             grid.AddField(&test);
//             test.bctype(M_BC_EXTRAP);

//             const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
//             const real_t freq[3]    = {2.0, 2.0, 2.0};
//             const real_t alpha[3]   = {1.0, 1.0, 1.0};
//             SetCosinus   field_init(sin_len, freq, alpha);
//             field_init(&grid, &test);

//             // pull the ghosts
//             grid.GhostPull(&test);

//             // IOH5 io("data_test");
//             // io(&grid, &test);
//             // io.dump_ghost(true);
//             // io(&grid, &test);

//             // create the solution field
//             Field sol("sol", 1);
//             grid.AddField(&sol);
//             SetCosinus field_sol(sin_len, freq, alpha, grid.interp());
//             field_sol(&grid, &sol);

//             // mask both the sol and the result
//             // MaskPhysBC mask(L);
//             // mask(&grid, &test);
//             // mask(&grid, &sol);

//             // now, we need to check
//             real_t          norm2, normi;
//             Error  error(grid.interp());
//             error.Norms(&grid, &test, &sol, err2 + il, erri + il);

//             m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
//         }
//         real_t conv2 = -log(err2[1] / err2[0]) / log(2);
//         real_t convi = -log(erri[1] / erri[0]) / log(2);

//         m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
//         ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
//         ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
//     }
// }

//==============================================================================================================================
TEST_F(ValidWaveletInterpolation, ghost_reconstruction_perper_dirichlet0_polynom) {
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        for (level_t il = 0; il < 2; il++) {
            bool period[3] = {true, true, true};
            period[id]     = false;
            lid_t L[3]     = {3, 3, 3};
            Grid  grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin1[3] = {1.0, 1.0, 1.0};
            real_t length1[3] = {1.0, 1.0, 1.0};
            Patch  p1(origin1, length1, il + BLVL + 1);
            real_t origin2[3] = {0.0, 0.0, 0.0};
            real_t length2[3] = {3.0, 3.0, 3.0};
            Patch  p2(origin2, length2, il + BLVL);

            list<Patch> patch{p1, p2};
            grid.Adapt(&patch);

            Field test("test", 1);
            grid.AddField(&test);
            // test.bctype(M_BC_EXTRAP);
            // put the DIRICHLET BC in the direction of interest
            test.bctype(M_BC_DIR, 0, 2 * id);
            test.bctype(M_BC_DIR, 0, 2 * id + 1);

            // create the initial field
            lid_t  deg[3]     = {0, 0, 0};
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
            real_t norm2, normi;
            Error  error(grid.interp());
            error.Norms(&grid, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
    }
}

//==============================================================================================================================
TEST_F(ValidWaveletInterpolation, ghost_reconstruction_perper_neuman0_cos) {
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        for (level_t il = 0; il < 2; il++) {
            bool period[3] = {true, true, true};
            period[id]     = false;
            lid_t L[3]     = {3, 3, 3};
            Grid  grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            // create the patch refinement to refine the middle tree
            real_t origin1[3] = {1.0, 1.0, 1.0};
            real_t length1[3] = {1.0, 1.0, 1.0};
            Patch  p1(origin1, length1, il + BLVL + 1);
            real_t origin2[3] = {0.0, 0.0, 0.0};
            real_t length2[3] = {3.0, 3.0, 3.0};
            Patch  p2(origin2, length2, il + BLVL);

            list<Patch> patch{p1, p2};
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
            freq[id]         = 1.0;
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
            real_t norm2, normi;
            Error  error(grid.interp());
            error.Norms(&grid, &test, m_ptr<const Field>(&sol), err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
    }
}
