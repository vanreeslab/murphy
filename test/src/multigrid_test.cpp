#include <cmath>

#include "boundary.hpp"
#include "error.hpp"
#include "field.hpp"
#include "gtest/gtest.h"
#include "ioh5.hpp"
#include "laplacian.hpp"
#include "murphy.hpp"
#include "setvalues.hpp"
#include "multigrid.hpp"
#incldue <mpi.h>

#define DOUBLE_TOL 1e-9

class valid_Multigrid : public ::testing::Test {
   protected:
    Grid* grid_;

    int  lvl_         = 2;
    bool periodic_[3] = {false, false, false};
    int  l_[3]        = {1, 1, 1};

    Field* vort_;
    Field* res_;
    Field* psi_;
    Field* analytic_;

    void SetUp() override {};
    void TearDown() override {};
};

TEST_F(valid_Multigrid, perio_uniform_order2) {
    //-----------
    // get the comm size and return if we are NOT sequential:
    int commsize;
    MPI_Comm_size(MPI_COMM_WORLD, &commsize);
    if (commsize > 1) {
        return 0;
    }

    lvl_         = 3;
    periodic_[0] = true;
    periodic_[1] = true;
    periodic_[2] = true;
    grid_        = new Grid(lvl_, periodic_, l_, MPI_COMM_WORLD, NULL);
    //add the field
    vort_     = new Field("vorticity", 1);
    res_      = new Field("residual", 1);
    psi_      = new Field("streamfunc", 1);
    analytic_ = new Field("analytic", 1);
    grid_->AddField(vort_);
    grid_->AddField(res_);
    grid_->AddField(psi_);
    grid_->AddField(analytic_);

    vort_->bctype(M_BC_EVEN);
    res_->bctype(M_BC_EVEN);
    psi_->bctype(M_BC_EVEN);
    //-------------
    real_t freq[3]   = {2.0, 1.0, 1.0};
    real_t length[3] = {1.0 * l_[0], 1.0 * l_[1], 1.0 * l_[2]};
    // set the fields
    SetCosinus cosinus = SetCosinus(length, freq);
    cosinus(grid_, analytic_);
    SetLaplaCosinus laplacos = SetLaplaCosinus(length, freq);
    laplacos(grid_, vort_);

    // initiate the solver
    Multigrid* poisson = new Multigrid(grid_, 0, vort_, psi_, res_);
    poisson->Solve();

    IOH5 dump("data");
    dump(grid_, vort_);
    dump(grid_, psi_);
    dump(grid_, analytic_);

    real_t          norm2, normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_, psi_, analytic_, &norm2, &normi);

    delete (poisson);

    printf("we have some errors: e_2 = %e, ei = %e\n", norm2, normi);

    // ASSERT_LE(norm2, normi);
    // ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    // ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);

    delete (vort_);
    delete (res_);
    delete (psi_);
    delete (analytic_);
    delete (grid_);
}

// TEST_F(valid_Multigrid, unbounded_uniform_order2) {
//     //-----------
//     //-------------
//     vort_->bctype(M_BC_ZERO);
//     res_->bctype(M_BC_ZERO);
//     psi_->bctype(M_BC_EXTRAP_4);

//     real_t shift     = M_PI / 4096.0;
//     real_t sigma[3]  = {0.05, 0.05, 0.05};
//     real_t center[3] = {0.5 + shift, 0.5 - shift, 0.5 + shift};
//     // rhs = exponential
//     SetExponential expo = SetExponential(center, sigma, -1.0);
//     expo(grid_, vort_);
//     // analytic = erf
//     SetErf erf = SetErf(center, sigma, 1.0);
//     erf(grid_, analytic_);

//     // initiate the solver
//     Multigrid* poisson = new Multigrid(grid_, 0, vort_, psi_, res_);
//     poisson->Solve();

//     IOH5 dump("data");
//     dump(grid_, vort_);
//     dump(grid_, psi_);
//     dump(grid_, analytic_);

//     real_t          norm2, normi;
//     ErrorCalculator myerr;
//     myerr.Norms(grid_, psi_, analytic_, &norm2, &normi);

//     delete (poisson);

//     printf("we have some errors: e_2 = %e, ei = %e\n", norm2, normi);

//     // ASSERT_LE(norm2, normi);
//     // ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
//     // ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
// }