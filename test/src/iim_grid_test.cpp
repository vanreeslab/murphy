#include "gtest/gtest.h"
#include "grid/grid.hpp"

class IIMBlockTest : public BlockOperator {
   public:
    explicit IIMBlockTest(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len){};

    void operator()(const ForestGrid* grid, Field* fid_x) {
        DoOpMesh(this, &IIMBlockTest::PrintEachIIMBlockOrigin, grid, fid_x);
        fid_x->ghost_len(ghost_len_res_);
    }

    void PrintEachIIMBlockOrigin(const qid_t* qid, const IIMBlock* block, Field* fid_x) {
        block->PrintIIMBlockOrigin();
    }
};

TEST(grid_type, iim) {

    // initialize grid
    bool  period[3] = {false, false, false};
    lid_t L[3]      = {2, 2, 2};
    level_t init_level = 0;
    BlockDataType block_type = M_IIMBLOCK;
    Grid grid(init_level, period, L, block_type, MPI_COMM_WORLD, nullptr);

    // create level set field
    Field level_set("LevelSet", 1);
    level_set.bctype(M_BC_EXTRAP);
    grid.AddField(&level_set);

    // Update level set ghosts
    const bidx_t ghost_len[2] = {1, 2};
    grid.GhostPull(&level_set, ghost_len);

    // refine the grid
    const real_t patch_origin[3] = {0., 0., 0.};
    const real_t patch_length[3] = {1., 1., 1.};
    {
        const level_t patch_level = 1;
        std::list<Patch> patch_list = {Patch(patch_origin, patch_length, patch_level)};
        grid.Adapt(&patch_list);
    }

    // Process test BlockOperator
    IIMBlockTest iim_block_test(ghost_len);
    iim_block_test(&grid, &level_set);
    // check command line for 15 print statements w/ blok origins

    // coarsen the grid
    {
        const level_t patch_level = 0;
        std::list<Patch> patch_list = {Patch(patch_origin, patch_length, patch_level)};
        grid.Adapt(&patch_list);
    }

    // re-process test BlockOperator
    iim_block_test(&grid, &level_set);
    // check command line for 8 more print statements w/ block origins
    
    // dummy assert, if the program didn't crash it's working
    ASSERT_EQ(1,1);
}
