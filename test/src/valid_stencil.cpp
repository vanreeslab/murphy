#include <cmath>
#include <list>

#include "gtest/gtest.h"
#include "core/macros.hpp"

#include "grid/grid.hpp"
#include "grid/field.hpp"
#include "core/types.hpp"
#include "operator/advection.hpp"
#include "operator/error.hpp"

#include "tools/ioh5.hpp"

#define DOUBLE_TOL 1e-13

class valid_Stencil : public ::testing::Test {
    void SetUp() override {};
    void TearDown() override {};
};

using std::list;
using std::string;

TEST_F(valid_Stencil, advection_periodic_cosinus) {
    for (lda_t id = 0; id < 3; ++id) {
    // init the errors
    real_t erri_adv_center_2[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_center_4[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_center_6[3] = {0.0, 0.0, 0.0};

    real_t erri_adv_weno_3[3] = {0.0, 0.0, 0.0};
    real_t erri_adv_weno_5[3] = {0.0, 0.0, 0.0};

    // setup the mesh
    bool  period[3] = {true, true, true};
    lid_t L[3]      = {3, 3, 3};

    // see if we run the tests
    bool do_center_2 = false;  // = grid.NGhostFront() >= 1 && grid.NGhostBack() >= 1;
    bool do_center_4 = false;  // = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
    bool do_center_6 = false;  // = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;
    bool do_weno_3   = false;
    bool do_weno_5   = false;

    for (level_t il = 1; il < 3; ++il) {
        Grid grid(il, period, L, MPI_COMM_WORLD, nullptr);
        // create the patch refinement to refine the middle tree
        real_t      origin1[3] = {1.0, 1.0, 1.0};
        real_t      length1[3] = {1.0, 1.0, 1.0};
        Patch       p1(origin1, length1, il + 1);
        real_t      origin2[3] = {0.0, 0.0, 0.0};
        real_t      length2[3] = {L[0], L[1], L[2]};
        Patch       p2(origin2, length2, il);
        list<Patch> patch{p1, p2};
        grid.Adapt(&patch);

        do_center_2 = grid.NGhostFront() >= 1 && grid.NGhostBack() >= 1;
        do_center_4 = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
        do_center_6 = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;
        do_weno_3   = grid.NGhostFront() >= 2 && grid.NGhostBack() >= 2;
        do_weno_5   = grid.NGhostFront() >= 3 && grid.NGhostBack() >= 3;

        // create the test file
        string fieldName = "field" + std::to_string(il);
        string solName   = "sol" + std::to_string(il);
        string velName   = "vel" + std::to_string(il);
        string diffName  = "deriv" + std::to_string(il);
        Field  test(fieldName, 1);
        Field  sol(solName, 1);
        Field  vel(velName, 3);
        Field  dtest(diffName, 1);
        grid.AddField(&test);
        grid.AddField(&sol);
        grid.AddField(&vel);
        grid.AddField(&dtest);

        // set a constant velocity
        const lid_t  deg[3]   = {0, 0, 0};
        const real_t dir[3]   = {1.0, 0.0, 0.0};
        const real_t shift[3] = {0.0, 0.0, 0.0};
        SetPolynom   vel_init(deg, dir, shift, grid.interp());
        vel_init(&grid, &vel, id);  // put 1.0 in the indicated direction only

        // put a cos -> the field: cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x) + cos(2*pi*freq[0]/L[0] * x)
        const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
        const real_t freq[3]    = {2.0, 2.0, 2.0};
        const real_t alpha[3]   = {1.0, 1.0, 1.0};
        SetCosinus   field_init(sin_len, freq, alpha);
        field_init(&grid, &test);

        // -> the solution: u* df/dx + v * df/dy + w*df/dz
        const real_t alpha_sol_0[3] = {2.0 * M_PI * freq[0] / L[0] * (id == 0),
                                       2.0 * M_PI * freq[1] / L[1] * (id == 1),
                                       2.0 * M_PI * freq[2] / L[2] * (id == 2)};
        SetSinus     sol_init(sin_len, freq, alpha_sol_0, grid.interp());
        sol_init(&grid, &sol);

        // test advection centered order 2
        if (do_center_2) {
            Advection<M_ADV_CENTER, 2> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            ErrorCalculator error;
            error.Normi(&grid, &dtest, &sol, erri_adv_center_2 + il);
            m_log("M_ADV_CENTER - 2: checking res = %f, ei = %e", std::pow(2, il), erri_adv_center_2[il]);
        }
        if (do_center_4) {
            Advection<M_ADV_CENTER, 4> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            ErrorCalculator error;
            error.Normi(&grid, &dtest, &sol, erri_adv_center_4 + il);
            m_log("M_ADV_CENTER - 4: checking res = %f, ei = %e", std::pow(2, il), erri_adv_center_4[il]);
        }
        if (do_center_6) {
            Advection<M_ADV_CENTER, 6> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            ErrorCalculator error;
            error.Normi(&grid, &dtest, &sol, erri_adv_center_6 + il);
            m_log("M_ADV_CENTER - 6: checking res = %f, ei = %e", std::pow(2, il), erri_adv_center_6[il]);
        }
        if (do_weno_3) {
            Advection<M_ADV_WENO_VEL, 3> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            ErrorCalculator error;
            error.Normi(&grid, &dtest, &sol, erri_adv_weno_3 + il);
            m_log("M_ADV_WENO_VEL - 3: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_3[il]);
        }
        if (do_weno_5) {
            Advection<M_ADV_WENO_VEL, 5> adv(&vel);
            adv(&grid, &test, &dtest);
            // now, we need to check
            ErrorCalculator error;
            error.Normi(&grid, &dtest, &sol, erri_adv_weno_5 + il);
            m_log("M_ADV_WENO_VEL - 5: checking res = %f, ei = %e", std::pow(2, il), erri_adv_weno_5[il]);
        }
    }
    m_log("in dir = %d did tests %d %d %d %d %d", do_center_2, do_center_4, do_center_6, do_weno_3, do_weno_5);
    // real_t conv2 = -log(err2[2] / err2[1]) / log(2);
    if (do_center_2) {
        real_t convi = -log(erri_adv_center_2[2] / erri_adv_center_2[1]) / log(2);
        m_log("M_ADV_CENTER - 2: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 2) - 0.1);
    }
    if (do_center_4) {
        real_t convi = -log(erri_adv_center_4[2] / erri_adv_center_4[1]) / log(2);
        m_log("M_ADV_CENTER - 4: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 4) - 0.1);
    }
    if (do_center_6) {
        real_t convi = -log(erri_adv_center_6[2] / erri_adv_center_6[1]) / log(2);
        m_log("M_ADV_CENTER - 6: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 6) - 0.1);
    }
    if (do_weno_3) {
        real_t convi = -log(erri_adv_weno_3[2] / erri_adv_weno_3[1]) / log(2);
        m_log("M_ADV_WENO - 3: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.1);
    }
    if (do_weno_5) {
        real_t convi = -log(erri_adv_weno_5[2] / erri_adv_weno_5[1]) / log(2);
        m_log("M_ADV_WENO - 5: the convergence orders are: norm_i:%e", convi);
        ASSERT_GE(convi, m_min(M_WAVELET_N - 1, 3) - 0.1);
    }
    }
}

// TEST_F(valid_Stencil, laplacian_periodic_cosinus) {
//     for (lda_t id = 0; id < 3; ++id) {
//         // init the errors
//         real_t erri[3] = {0.0, 0.0, 0.0};
//         real_t err2[3] = {0.0, 0.0, 0.0};

//         // setup the mesh
//         bool  period[3] = {true, true, true};
//         lid_t L[3]      = {3, 3, 3};

//         for (level_t il = 1; il < 3; ++il) {
//             Grid grid(il, period, L, MPI_COMM_WORLD, nullptr);

//             // create the patch refinement to refine the middle tree
//             real_t      origin1[3] = {1.0, 1.0, 1.0};
//             real_t      length1[3] = {1.0, 1.0, 1.0};
//             Patch       p1(origin1, length1, il + 1);
//             real_t      origin2[3] = {0.0, 0.0, 0.0};
//             real_t      length2[3] = {L[0], L[1], L[2]};
//             Patch       p2(origin2, length2, il);
//             list<Patch> patch{p1, p2};
//             grid.Adapt(&patch);

//             // create the test file
//             string fieldName = "field" + std::to_string(id) + "__" + std::to_string(il);
//             string solName   = "sol" + std::to_string(id) + "__" + std::to_string(il);
//             string diffName  = "diff" + std::to_string(id) + "__" + std::to_string(il);
//             Field  test(fieldName, 3);
//             Field  sol(solName, 3);
//             Field  diff(diffName, 3);
//             grid.AddField(&test);
//             grid.AddField(&sol);
//             grid.AddField(&diff);

//             // put a sinus -> the field
//             const real_t sin_len[3] = {(real_t)L[0], (real_t)L[1], (real_t)L[2]};
//             const real_t freq[3]    = {3.0, 1.0, 1.0};
//             const real_t alpha[3]    = {1.0, 1.0, 1.0};
//             SetCosinus     field_init(sin_len, freq,alpha);
//             field_init(&grid, &test);

//             // -> the solution
//             const real_t alpha_sol[3] = {-pow((2 * M_PI * freq[0]) / L[0], 2), -pow((2 * M_PI * freq[1]) / L[1], 2), -pow((2 * M_PI * freq[2]) / L[2], 2)};
//             SetCosinus     sol_init(sin_len, freq,alpha_sol,grid.interp());
//             sol_init(&grid, &sol);

//             // compue the laplacian
//             LaplacianCross<3> diffusion;
//             diffusion(&grid,&test,&diff);

//             // grid.GhostPull(&diff);
//             // IOH5 io("data_test");
//             // io(&grid,&diff);
//             // io(&grid,&sol);

//             // now, we need to check
//             ErrorCalculator error;
//             error.Norms(&grid, &diff, &sol, err2 + il, erri + il);

//             m_log("checking in dim %d: res = %f, ei = %e e2 = %e", id, std::pow(2, il), erri[il], err2[il]);
//         }
//         real_t conv2 = -log(err2[2] / err2[1]) / log(2);
//         real_t convi = -log(erri[2] / erri[1]) / log(2);

//         m_log("the convergence orders are: norm_2:%e norm_i:%e", conv2, convi);
//         ASSERT_GE(conv2, m_min(M_WAVELET_N-2,2) - 0.1);
//         ASSERT_GE(convi, m_min(M_WAVELET_N-2,2) - 0.1);
//     }
// }

// // #endif

// // TEST_F(valid_Stencil, laplacian_o2) {
// //     vort_->bctype(M_BC_EXTRAP_3);

// //     // x^2 + y^2 + z^2
// //     real_t     dir[3]  = {1.0, 1.0, 1.0};
// //     lid_t      deg[3]  = {2, 2, 2};
// //     SetPolynom polynom = SetPolynom(deg, dir);
// //     polynom(grid_, vort_);

// //     // 2 + 2 + 2
// //     real_t     dir2[3]  = {2.0, 2.0, 2.0};
// //     lid_t      deg2[3]  = {0, 0, 0};
// //     SetPolynom polynom2 = SetPolynom(deg2, dir2);
// //     polynom2(grid_, sol_);

// //     LaplacianCross<3> lapla = LaplacianCross<3>();
// //     lapla(grid_, vort_, diff_);

// //     real_t          norm2, normi;
// //     ErrorCalculator myerr;
// //     myerr.Norms(grid_, diff_, sol_, &norm2, &normi);

// //     // ASSERT_LE(norm2, normi);
// //     ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
// //     ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
// // }

// // TEST_F(valid_Stencil, laplacian_o4) {
// //     vort_->bctype(M_BC_EXTRAP_5);

// //     // x^4 + y^4 + z^4
// //     real_t     dir[3]  = {1.0, 1.0, 1.0};
// //     lid_t      deg[3]  = {4, 4, 4};
// //     SetPolynom polynom = SetPolynom(deg, dir);
// //     polynom(grid_, vort_);

// //     // 12 x^2 + 12 y^2 + 12 z^2
// //     real_t     dir2[3]  = {12.0, 12.0, 12.0};
// //     lid_t      deg2[3]  = {2, 2, 2};
// //     SetPolynom polynom2 = SetPolynom(deg2, dir2);
// //     polynom2(grid_, sol_);

// //     LaplacianCross<5> lapla = LaplacianCross<5>();
// //     lapla(grid_, vort_, diff_);

// //     real_t          norm2, normi;
// //     ErrorCalculator myerr;
// //     myerr.Norms(grid_, diff_, sol_, &norm2, &normi);

// //     // ASSERT_LE(norm2, normi);
// //     ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
// //     ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
// // }

// // #define N_CONV 2

// // TEST_F(valid_Stencil, convergence_laplacian_o2_boundary3) {
// //     vort_->bctype(M_BC_EXTRAP_3);

// //     real_t erri[N_CONV] = {0.0};
// //     real_t origin1[3]   = {0.0, 0.0, 0.0};
// //     real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
// //     real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
// //     real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

// //     real_t normi[N_CONV];

// //     for (int il = lvl_; il < lvl_ + N_CONV; il++) {
// //         list<Patch> patches;
// //         patches.push_back(Patch(origin1, length1, il));
// //         patches.push_back(Patch(origin2, length2, il + 1));
// //         grid_->Adapt(&patches);

// //         // x^3 + y^4 + z^3
// //         real_t     dir[3]  = {1.0, 1.0, 1.0};
// //         lid_t      deg[3]  = {3, 4, 3};
// //         SetPolynom polynom = SetPolynom(deg, dir);
// //         polynom(grid_, vort_);

// //         // 2 + 2 + 2
// //         real_t     dir2[3]  = {6.0, 12.0, 6.0};
// //         lid_t      deg2[3]  = {1, 2, 1};
// //         SetPolynom polynom2 = SetPolynom(deg2, dir2);
// //         polynom2(grid_, sol_);

// //         LaplacianCross<3> lapla = LaplacianCross<3>();
// //         lapla(grid_, vort_, diff_);

// //         real_t          norm2;
// //         ErrorCalculator myerr;
// //         myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
// //         if (grid_->mpirank() == 0) {
// //             printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
// //         }
// //     }
// //     if (grid_->mpirank() == 0) {
// //         printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 1.95);
// //     }
// //     ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 1.95);
// // }

// // TEST_F(valid_Stencil, convergence_laplacian_o2_boundary4) {
// //     vort_->bctype(M_BC_EXTRAP_4);

// //     real_t erri[N_CONV] = {0.0};
// //     real_t origin1[3]   = {0.0, 0.0, 0.0};
// //     real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
// //     real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
// //     real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

// //     real_t normi[N_CONV];

// //     for (int il = lvl_; il < lvl_ + N_CONV; il++) {
// //         list<Patch> patches;
// //         patches.push_back(Patch(origin1, length1, il));
// //         patches.push_back(Patch(origin2, length2, il + 1));
// //         grid_->Adapt(&patches);

// //         // x^3 + y^4 + z^3
// //         real_t     dir[3]  = {1.0, 1.0, 1.0};
// //         lid_t      deg[3]  = {3, 4, 3};
// //         SetPolynom polynom = SetPolynom(deg, dir);
// //         polynom(grid_, vort_);

// //         // 2 + 2 + 2
// //         real_t     dir2[3]  = {6.0, 12.0, 6.0};
// //         lid_t      deg2[3]  = {1, 2, 1};
// //         SetPolynom polynom2 = SetPolynom(deg2, dir2);
// //         polynom2(grid_, sol_);

// //         LaplacianCross<3> lapla = LaplacianCross<3>();
// //         lapla(grid_, vort_, diff_);

// //         real_t          norm2;
// //         ErrorCalculator myerr;
// //         myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
// //         if (grid_->mpirank() == 0) {
// //             printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
// //         }
// //     }
// //     if (grid_->mpirank() == 0) {
// //         printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 3.99);
// //     }
// //     ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 3.99);
// // }

// // TEST_F(valid_Stencil, convergence_laplacian_o4_boundary4) {
// //     vort_->bctype(M_BC_EXTRAP_4);

// //     real_t erri[N_CONV] = {0.0};
// //     real_t origin1[3]   = {0.0, 0.0, 0.0};
// //     real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
// //     real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
// //     real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

// //     real_t normi[N_CONV];

// //     for (int il = lvl_; il < lvl_ + N_CONV; il++) {
// //         list<Patch> patches;
// //         patches.push_back(Patch(origin1, length1, il));
// //         patches.push_back(Patch(origin2, length2, il + 1));
// //         grid_->Adapt(&patches);

// //         // x^6 + y^5 + z^7
// //         real_t     dir[3]  = {1.0, 1.0, 1.0};
// //         lid_t      deg[3]  = {6, 5, 7};
// //         SetPolynom polynom = SetPolynom(deg, dir);
// //         polynom(grid_, vort_);

// //         // solution
// //         real_t     dir2[3]  = {30.0, 20.0, 42.0};
// //         lid_t      deg2[3]  = {4, 3, 5};
// //         SetPolynom polynom2 = SetPolynom(deg2, dir2);
// //         polynom2(grid_, sol_);

// //         LaplacianCross<5> lapla = LaplacianCross<5>();
// //         lapla(grid_, vort_, diff_);

// //         real_t          norm2;
// //         ErrorCalculator myerr;
// //         myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
// //         if (grid_->mpirank() == 0) {
// //             printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
// //         }
// //     }
// //     if (grid_->mpirank() == 0) {
// //         printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 3.82);
// //     }
// //     ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 3.82);
// // }

// // TEST_F(valid_Stencil, convergence_laplacian_o4_boundary5) {
// //     vort_->bctype(M_BC_EXTRAP_5);

// //     real_t erri[N_CONV] = {0.0};
// //     real_t origin1[3]   = {0.0, 0.0, 0.0};
// //     real_t length1[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};
// //     real_t origin2[3]   = {0.0, 0.0, 0.5 * l_[2]};
// //     real_t length2[3]   = {1.0 * l_[0], 1.0 * l_[1], 0.5 * l_[2]};

// //     real_t normi[N_CONV];

// //     for (int il = lvl_; il < lvl_ + N_CONV; il++) {
// //         list<Patch> patches;
// //         patches.push_back(Patch(origin1, length1, il));
// //         patches.push_back(Patch(origin2, length2, il + 1));
// //         grid_->Adapt(&patches);

// //         // x^6 + y^5 + z^7
// //         real_t     dir[3]  = {1.0, 1.0, 1.0};
// //         lid_t      deg[3]  = {6, 5, 7};
// //         SetPolynom polynom = SetPolynom(deg, dir);
// //         polynom(grid_, vort_);

// //         // solution
// //         real_t     dir2[3]  = {30.0, 20.0, 42.0};
// //         lid_t      deg2[3]  = {4, 3, 5};
// //         SetPolynom polynom2 = SetPolynom(deg2, dir2);
// //         polynom2(grid_, sol_);

// //         LaplacianCross<5> lapla = LaplacianCross<5>();
// //         lapla(grid_, vort_, diff_);

// //         real_t          norm2;
// //         ErrorCalculator myerr;
// //         myerr.Norms(grid_, diff_, sol_, &norm2, normi + (il - lvl_));
// //         if (grid_->mpirank() == 0) {
// //             printf("lvl = %d error = %e %e\n", il, norm2, normi[il - lvl_]);
// //         }
// //     }
// //     if (grid_->mpirank() == 0) {
// //         printf("tested convergence order: %e vs %e\n", normi[N_CONV - 2] / normi[N_CONV - 1], 7.77);
// //     }
// //     ASSERT_GT(normi[N_CONV - 2] / normi[N_CONV - 1], 7.77);
// // }
