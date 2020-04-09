#include "grid.hpp"

#include "gtest/gtest.h"

using std::pow;

class test_Grid : public ::testing::Test {
   protected:
    Grid* grid_;

    int  lvl_         = 2;
    bool periodic_[3] = {false, true, false};
    int  l_[3]        = {1, 2, 3};

    void SetUp() override {
        grid_ = new Grid(lvl_, periodic_, l_, MPI_COMM_WORLD, NULL);
    };
    void TearDown() override {
        delete (grid_);
    };
};

TEST_F(test_Grid, creation_non_null) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if(rank == 1) grid_ = nullptr;
    // create a grid
    ASSERT_NE(grid_, nullptr);
}
TEST_F(test_Grid, creation_num_dof) {
    // create a grid
    size_t numdof = l_[0] * l_[1] * l_[2] * pow(2, 3*lvl_) * (M_N * M_N * M_N);
    ASSERT_EQ(grid_->GlobalNumDof(), numdof);
}