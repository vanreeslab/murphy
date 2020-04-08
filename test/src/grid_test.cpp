#include "gtest/gtest.h"
#include "grid.hpp"

TEST(test_grid, creation)
{
        bool periodic[3] = {false, true, false};
        int  l[3]        = {1, 2, 3};

        // create a grid
        Grid*  grid = new Grid(2, periodic, l, MPI_COMM_WORLD, NULL);
        ASSERT_NE(grid,nullptr);
        delete (grid);
    
}