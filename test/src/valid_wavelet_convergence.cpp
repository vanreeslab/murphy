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

class ValidWaveletInterpolation : public ::testing::TestWithParam<int> {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

//==============================================================================================================================
TEST_P(ValidWaveletInterpolation, ghost_reconstruction_periodic_sin) {
    // get the ghost length
    const bidx_t ng           = GetParam();
    const bidx_t ghost_len[2] = {ng, ng};
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
            const real_t      sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t      freq[3]    = {3.0, 3.0, 3.0};
            lambda_setvalue_t sin_op     = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid, 0)(i0, i1, i2) = sin(2.0 * M_PI * freq[0] / sin_len[0] * pos[0]) +
                                                  sin(2.0 * M_PI * freq[1] / sin_len[1] * pos[1]) +
                                                  sin(2.0 * M_PI * freq[2] / sin_len[2] * pos[2]);
            };
            SetValue field_init(sin_op);
            field_init(&grid, &test);
            // SetSinus     field_init(sin_len, freq, alpha, grid.interp());
            // SetSinus field_init(sin_len, freq, alpha);

            // we have perfect ghost -> what is not set will be perfect
            // we force the ghost recomputation
            // test.ghost_status(false);

            // pull the ghosts
            grid.GhostPull(&test, ghost_len);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // // create the solution field
            // Field sol("sol", 1);
            // grid.AddField(&sol);
            // SetSinus field_sol(sin_len, freq, alpha, grid.interp());
            // field_sol(&grid, &sol);

            // now, we need to check the ghost
            lambda_error_t lambda_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                return sin(2.0 * M_PI * freq[0] / sin_len[0] * pos[0]) +
                       sin(2.0 * M_PI * freq[1] / sin_len[1] * pos[1]) +
                       sin(2.0 * M_PI * freq[2] / sin_len[2] * pos[2]);
            };
            // the error is computed in the ghost region (AS WELL)
            Error error(ghost_len);

            // sanity check
            error.Norms(&grid, il + 1 + BLVL, &test, &lambda_sol, err2 + il, erri + il);
            m_log("checking in dim %d on HIGH: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL + 1), erri[il], err2[il]);
            error.Norms(&grid, il + BLVL, &test, &lambda_sol, err2 + il, erri + il);
            m_log("checking in dim %d on LOW: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
// if NT==0, the lowest level has NO error
#if (M_WAVELET_NT == 0)
            ASSERT_NEAR(erri[il], 0.0, DOUBLE_TOL);
            ASSERT_NEAR(err2[il], 0.0, DOUBLE_TOL);
#endif

            // and now check for the every level
            error.Norms(&grid, &test, &lambda_sol, err2 + il, erri + il);
            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);

            grid.DeleteField(&test);
            // grid.DeleteField(&sol);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
    }
}
//==============================================================================================================================
TEST_P(ValidWaveletInterpolation, ghost_reconstruction_periodic_cos) {
    // get the ghost length
    const bidx_t ng           = GetParam();
    const bidx_t ghost_len[2] = {ng, ng};
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
            Field  test(fieldName, 1);
            grid.AddField(&test);

            // put a sinus
            const real_t      cos_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t      freq[3]    = {2.0, 3.0, 1.0};
            lambda_setvalue_t cos_op     = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid, 0)(i0, i1, i2) = cos(2.0 * M_PI * freq[0] / cos_len[0] * pos[0]) +
                                                  cos(2.0 * M_PI * freq[1] / cos_len[1] * pos[1]) +
                                                  cos(2.0 * M_PI * freq[2] / cos_len[2] * pos[2]);
            };
            SetValue field_init(cos_op);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test, ghost_len);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            // Field sol("sol", 3);
            // grid.AddField(&sol);
            // SetCosinus field_sol(sin_len, freq, alpha, grid.interp());
            // field_sol(&grid, &sol);

            // now, we need to check
            // now, we need to check the ghost
            lambda_error_t sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                return cos(2.0 * M_PI * freq[0] / cos_len[0] * pos[0]) +
                       cos(2.0 * M_PI * freq[1] / cos_len[1] * pos[1]) +
                       cos(2.0 * M_PI * freq[2] / cos_len[2] * pos[2]);
            };
            Error error(ghost_len);
            // sanity check
            error.Norms(&grid, il + 1 + BLVL, &test, &sol, err2 + il, erri + il);
            m_log("checking in dim %d on HIGH: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
#if (M_WAVELET_NT == 0)
            error.Norms(&grid, il + BLVL, &test, &sol, err2 + il, erri + il);
            m_log("checking in dim %d on LOW: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
            // if NT==0, the lowest level has NO error
            ASSERT_NEAR(erri[il], 0.0, DOUBLE_TOL);
            ASSERT_NEAR(err2[il], 0.0, DOUBLE_TOL);
#endif

            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

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
TEST_P(ValidWaveletInterpolation, ghost_reconstruction_extrap_cos) {
    // get the ghost length
    const bidx_t ng           = GetParam();
    const bidx_t ghost_len[2] = {ng, ng};
    // init the errors
    real_t erri[2];
    real_t err2[2];
    for (lda_t id = 0; id < 3; id++) {
        // setup the mesh
        bool  period[3] = {false, false, false};
        lid_t L[3]      = {2, 2, 2};

        for (level_t il = 0; il < 2; ++il) {
            Grid grid(il + BLVL, period, L, MPI_COMM_WORLD, nullptr);

            // create the patch refinement to refine the middle tree
            real_t origin1[3] = {1.0, 1.0, 0.0};
            real_t length1[3] = {1.0, 1.0, 2.0};
            Patch  p1(origin1, length1, il + BLVL + 1);
            real_t origin2[3] = {0.0, 0.0, 0.0};
            real_t length2[3] = {2.0, 2.0, 2.0};
            Patch  p2(origin2, length2, il + BLVL);

            list<Patch> patch{p1, p2};
            grid.Adapt(&patch);

            // create the test file
            string fieldName = "cosinus" + std::to_string(id) + "__" + std::to_string(il);
            Field  test(fieldName, 1);
            test.bctype(M_BC_EXTRAP);
            grid.AddField(&test);

            // put a sinus
            const real_t      cos_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
            const real_t      freq[3]    = {2.0, 3.0, 1.0};
            lambda_setvalue_t cos_op     = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid,0)(i0, i1, i2) = cos(2.0 * M_PI * freq[0] / cos_len[0] * pos[0]) +
                                                        cos(2.0 * M_PI * freq[1] / cos_len[1] * pos[1]) +
                                                        cos(2.0 * M_PI * freq[2] / cos_len[2] * pos[2]);
            };
            SetValue field_init(cos_op);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test, ghost_len);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            // Field sol("sol", 3);
            // grid.AddField(&sol);
            // SetCosinus field_sol(sin_len, freq, alpha, grid.interp());
            // field_sol(&grid, &sol);

            // now, we need to check
            // now, we need to check the ghost
            lambda_error_t sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                return cos(2.0 * M_PI * freq[0] / cos_len[0] * pos[0]) +
                       cos(2.0 * M_PI * freq[1] / cos_len[1] * pos[1]) +
                       cos(2.0 * M_PI * freq[2] / cos_len[2] * pos[2]);
            };
            Error error(ghost_len);
            // sanity check
            //             error.Norms(&grid, il + 1 + BLVL, &test, &sol, err2 + il, erri + il);
            //             m_log("checking in dim %d on HIGH: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
            // #if (M_WAVELET_NT == 0)
            //             error.Norms(&grid, il + BLVL, &test, &sol, err2 + il, erri + il);
            //             m_log("checking in dim %d on LOW: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
            //             // if NT==0, the lowest level has NO error
            //             ASSERT_NEAR(erri[il], 0.0, DOUBLE_TOL);
            //             ASSERT_NEAR(err2[il], 0.0, DOUBLE_TOL);
            // #endif

            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

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
TEST_P(ValidWaveletInterpolation, ghost_reconstruction_perper_dirichlet0_polynom) {
    // get the ghost length
    const bidx_t ng           = GetParam();
    const bidx_t ghost_len[2] = {ng, ng};

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

            // create the initial field L/2^(N+2) - (L/2-x)^(N+2)
            // lid_t  deg[3]     = {0, 0, 0};
            // real_t dir[3]     = {0.0, 0.0, 0.0};
            // real_t shift[3]   = {0.0, 0.0, 0.0};
            // deg[id]           = M_WAVELET_N + 2;
            // dir[id]           = -1.0;
            // dir[(id + 1) % 3] = pow(L[id] / 2.0, deg[id]);
            // shift[id]         = L[id] / 2.0;

            // SetPolynom field_init(deg, dir, shift);
            // field_init(&grid, &test);

            lambda_setvalue_t pol_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                const short_t deg               = M_WAVELET_N + 2;
                block->data(fid, 0)(i0, i1, i2) = pow(L[id] / 2.0, deg) - pow(L[id] / 2.0 - pos[id], deg);
            };
            SetValue field_init(pol_op);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test, ghost_len);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            // Field sol("sol", 1);
            // grid.AddField(&sol);
            // SetPolynom field_sol(deg, dir, shift, grid.interp());
            // field_sol(&grid, &sol);
            lambda_error_t sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                const short_t deg = M_WAVELET_N + 2;
                return pow(L[id] / 2.0, deg) - pow(L[id] / 2.0 - pos[id], deg);
            };

            // mask both the sol and the result
            // MaskPhysBC mask(L);
            // mask(&grid, &test);
            // mask(&grid, &sol);

            // now, we need to check
            real_t norm2, normi;
            Error  error(ghost_len);
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

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
TEST_P(ValidWaveletInterpolation, ghost_reconstruction_perper_neuman0_cos) {
    // get the ghost length
    const bidx_t ng           = GetParam();
    const bidx_t ghost_len[2] = {ng, ng};

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
            // real_t freq[3]   = {0.0, 0.0, 0.0};
            // real_t coslen[3] = {(real_t)L[0], (real_t)L[1],(real_t) L[2]};
            // real_t alpha[3]  = {0.0, 0.0, 0.0};
            // freq[id]         = 1.0;
            // alpha[id]        = 1.0;

            // SetCosinus field_init(coslen, freq, alpha);
            // field_init(&grid, &test);
            lambda_setvalue_t cos_op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                block->data(fid,0)(i0, i1, i2) = cos(2.0 * M_PI / L[id] * pos[id]);
            };
            SetValue field_init(cos_op);
            field_init(&grid, &test);

            // pull the ghosts
            grid.GhostPull(&test, ghost_len);

            // IOH5 io("data_test");
            // io(&grid, &test);
            // io.dump_ghost(true);
            // io(&grid, &test);

            // create the solution field
            // Field sol("sol", 1);
            // grid.AddField(&sol);
            // SetCosinus field_sol(coslen, freq, alpha, grid.interp());
            // field_sol(&grid, &sol);
            lambda_error_t sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
                // get the position
                real_t pos[3];
                block->pos(i0, i1, i2, pos);

                // call the function
                return cos(2.0 * M_PI / L[id] * pos[id]);
            };

            // mask both the sol and the result
            // MaskPhysBC mask(L);
            // mask(&grid, &test);
            // mask(&grid, &sol);

            // now, we need to check
            real_t norm2, normi;
            Error  error(ghost_len);
            error.Norms(&grid, &test, &sol, err2 + il, erri + il);

            m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il + BLVL), erri[il], err2[il]);
        }
        real_t conv2 = -log(err2[1] / err2[0]) / log(2);
        real_t convi = -log(erri[1] / erri[0]) / log(2);

        m_log("==> the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
        ASSERT_GE(conv2, M_WAVELET_N - ORDER2_TOL);
        ASSERT_GE(convi, M_WAVELET_N - ORDERI_TOL);
    }
}

INSTANTIATE_TEST_SUITE_P(ValidGhost,
                         ValidWaveletInterpolation,
                         testing::Range(1, M_GS));