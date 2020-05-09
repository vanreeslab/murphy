#include "grid.hpp"
#include "setvalues.hpp"
#include "error.hpp"
#include "gtest/gtest.h"

#define DOUBLE_TOL 1e-13

using std::pow;

class test_Grid_Uniform : public ::testing::Test {
   protected:
    Grid* grid_;

    int  lvl_         = 0;
    bool periodic_[3] = {false, true, false};
    int  l_[3]        = {1, 2, 3};

    void SetUp() override {
        grid_ = new Grid(lvl_, periodic_, l_, MPI_COMM_WORLD, NULL);
    };
    void TearDown() override {
        delete (grid_);
    };
};

TEST_F(test_Grid_Uniform, creation_non_null) {
    // create a grid
    ASSERT_NE(grid_, nullptr);
}

TEST_F(test_Grid_Uniform, creation_num_dof) {
    // create a grid
    size_t numdof = l_[0] * l_[1] * l_[2] * pow(2, 3 * lvl_) * (M_N * M_N * M_N);
    ASSERT_EQ(grid_->GlobalNumDof(), numdof);

    Field* vort = new Field("vorticity", 3);
}

TEST_F(test_Grid_Uniform, fields_num_dof) {
    ASSERT_EQ(grid_->NField(), 0);

    // add a new field
    Field* vort = new Field("vort", 3);
    grid_->AddField(vort);
    ASSERT_EQ(grid_->NField(), 1);

    // add the same field -> shouldn't change
    Field* vort2 = new Field("vort", 3);
    grid_->AddField(vort2);
    ASSERT_EQ(grid_->NField(), 1);

    // add another field
    Field* velocity = new Field("velocity", 3);
    grid_->AddField(velocity);
    ASSERT_EQ(grid_->NField(), 2);

    // remove a field
    grid_->DeleteField(vort);
    ASSERT_EQ(grid_->NField(), 1);

    delete (vort);
    delete (vort2);
    delete (velocity);
}

TEST_F(test_Grid_Uniform, refine_coarsen_uniform_1_level) {
    // create a field
    Field* vort = new Field("vorticity", 1);
    grid_->AddField(vort);
    // set a Gaussian
    real_t      center[3] = {l_[0] * 0.5, l_[1] * 0.5, l_[2] * 0.5};
    SetGaussian gaussian  = SetGaussian(0.1, center);
    gaussian(grid_, vort);
    vort->bctype(M_BC_EXTRAP_5);
    // refine the grid
    grid_->Refine(1);
    // coarsen the grid
    grid_->Coarsen(1);
    
    // compare
    Field* sol = new Field("solution", 1);
    grid_->AddField(sol);
    gaussian(grid_, sol);

    real_t norm2,normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_,vort,sol,&norm2,&normi);

    // the 2 norm is always smaller
    ASSERT_LE(norm2,normi);
    // we retrieve the field
    ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);

    // destroy the fields
    delete (vort);
    delete(sol);

}

TEST_F(test_Grid_Uniform, ref_coars_ref_coarse_uniform) {
    // create a field
    Field* vort = new Field("vorticity", 3);
    grid_->AddField(vort);
    // set a Gaussian
    real_t      center[3] = {l_[0] * 0.5, l_[1] * 0.5, l_[2] * 0.5};
    SetGaussian gaussian  = SetGaussian(0.1, center);
    gaussian(grid_, vort);
    vort->bctype(M_BC_EXTRAP_5);
    // refine the grid
    grid_->Refine(2);
    // coarsen the grid
    grid_->Coarsen(1);
    // coarsen the grid
    grid_->Refine(1);
    // coarsen the grid
    grid_->Coarsen(2);
    
    // compare
    Field* sol = new Field("solution", 3);
    grid_->AddField(sol);
    gaussian(grid_, sol);

    real_t norm2,normi;
    ErrorCalculator myerr;
    myerr.Norms(grid_,vort,sol,&norm2,&normi);

    // the 2 norm is always smaller
    ASSERT_LE(norm2,normi);
    // we retrieve the field
    ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);

    // destroy the fields
    delete (vort);
    delete(sol);

}