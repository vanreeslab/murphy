#include <cmath>

#include "boundary.hpp"
#include "error.hpp"
#include "field.hpp"
#include "gtest/gtest.h"
#include "laplacian.hpp"
#include "murphy.hpp"
#include "setvalues.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-9

class valid_Stencil : public ::testing::Test {
   protected:
    Grid* grid_;

    int  lvl_         = 1;
    bool periodic_[3] = {false, false, false};
    int  l_[3]        = {1, 2, 3};

    Field* vort_;
    Field* diff_;
    Field* sol_;

    void SetUp() override {
        grid_ = new Grid(lvl_, periodic_, l_, MPI_COMM_WORLD, NULL);

        //add the field
        vort_ = new Field("vorticity", 1);
        diff_ = new Field("diffusion", 1);
        sol_  = new Field("solution", 1);
        grid_->AddField(vort_);
        grid_->AddField(diff_);
        grid_->AddField(sol_);
    };
    void TearDown() override {
        delete (vort_);
        delete (diff_);
        delete (sol_);
        delete (grid_);
    };
};

TEST_F(valid_Stencil, laplacian_o2) {
    vort_->bctype(M_BC_EXTRAP_3);

    // x^2 + y^2 + z^2
    real_t     dir[3]  = {1.0, 1.0, 1.0};
    lid_t      deg[3]  = {2, 2, 2};
    SetPolynom polynom = SetPolynom(deg, dir);
    polynom(grid_, vort_);

    // 2 + 2 + 2
    real_t     dir2[3]  = {2.0, 2.0, 2.0};
    lid_t      deg2[3]  = {0, 0, 0};
    SetPolynom polynom2 = SetPolynom(deg2, dir2);
    polynom2(grid_, sol_);

    LaplacianCross<3> lapla(grid_);
    lapla(vort_, diff_);

    real_t          norm2, normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_, diff_, sol_, &norm2, &normi);

    // ASSERT_LE(norm2, normi);
    ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
}


TEST_F(valid_Stencil, laplacian_o4) {
    vort_->bctype(M_BC_EXTRAP_5);

    // x^4 + y^4 + z^4
    real_t     dir[3]  = {1.0, 1.0, 1.0};
    lid_t      deg[3]  = {4, 4, 4};
    SetPolynom polynom = SetPolynom(deg, dir);
    polynom(grid_, vort_);

    // 12 x^2 + 12 y^2 + 12 z^2
    real_t     dir2[3]  = {12.0, 12.0, 12.0};
    lid_t      deg2[3]  = {2, 2, 2};
    SetPolynom polynom2 = SetPolynom(deg2, dir2);
    polynom2(grid_, sol_);

    LaplacianCross<5> lapla(grid_);
    lapla(vort_, diff_);

    real_t          norm2, normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_, diff_, sol_, &norm2, &normi);

    // ASSERT_LE(norm2, normi);
    ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
}
