#include "gtest/gtest.h"
#include "gtest/gtest-death-test.h"
#include "core/doop.hpp"
#include "grid/grid.hpp"
#include <typeinfo>

// Count blocks and check their type
class BlockChecker {
   protected:
    bool types_check_ = true;
    short_t block_count_ = 0;
    const std::type_info& correct_type_;

   public:
    BlockChecker(const std::type_info& correct_type): correct_type_(correct_type) {}

    short_t block_count() const {return block_count_;}
    bool    types_check() const {return types_check_;}

    void CheckBlocks(ForestGrid* grid) {
        DoOpMesh(this, &BlockChecker::DoMagic, grid); 
        short_t global_block_count = 0;
        bool    global_types_check = true;
        MPI_Allreduce(&(this->block_count_), &global_block_count, 1, MPI_SHORT,    MPI_SUM,  MPI_COMM_WORLD);
        MPI_Allreduce(&(this->types_check_), &global_types_check, 1, MPI_CXX_BOOL, MPI_LAND, MPI_COMM_WORLD);
        block_count_ = global_block_count;
        types_check_ = global_types_check;
    }

    void DoMagic(const qid_t* quad, CartBlock* block) {
        block_count_++;
        types_check_ &= typeid(*block) == correct_type_;
    }
};

// Does something that works on any CartBlock
class CartBlockOperator {
   protected:
    short_t class_member = 0;

   public:
    CartBlockOperator() {}

    void OperateConst(const qid_t* quad, CartBlock* block, const Field* fid) const {
        block->data(fid, 0)(0, 0, 0) += 1.0;
    }

    void OperateNonConst(const qid_t*quad, CartBlock* block, const Field* fid) {
        block->data(fid, 0)(0, 0, 0) += 1.0;
        class_member++;
    }
};

// Do something that requires a GridBlock
class GridBlockOperator {
   protected:
    short_t class_member = 0;

   public:
    GridBlockOperator() {}

    void OperateConst(const qid_t* quad, GridBlock* block, const Field* fid) const {
        block->status_level(M_ADAPT_FINER);
    }

    void OperateNonConst(const qid_t*quad, GridBlock* block, const Field* fid) {
        block->status_level(M_ADAPT_FINER);
        class_member++;
    }
};

// setup a ForestGrid with a given block type
ForestGrid* CreateForestGrid(level_t init_lvl, const bool period[3], const lid_t L[3], BlockDataType block_type) {    
    ForestGrid* grid = new ForestGrid(init_lvl, period, L, block_type, MPI_COMM_WORLD);
    p8est_iterate(grid->p4est_forest(), nullptr, nullptr, get_cback_CreateBlock(grid->block_type()), nullptr, nullptr, nullptr);
    grid->SetupP4estMeshAndGhost();
    return grid;
}

// setup a ForestGrid with a given block type
void DestroyForestGrid(ForestGrid* grid) {
    grid->DestroyP4estMeshAndGhost();
    p8est_iterate(grid->p4est_forest(), nullptr, nullptr, get_cback_DestroyBlock(grid->block_type()), nullptr, nullptr, nullptr);
    delete grid;
}

TEST(BlockTypeDeathTest, doop) {

    // GTEST_FLAG_SET(death_test_style, "threadsafe"); - this should compile, but doesn't?
    // At the moment it's difficult to run death tests with MPI, so those are commented out.
    // They can be tested safely on a single process with :
    // ./murphy_test --gtest_filter="BlockTypeDeathTest.doop" --gtest_death_test_style=threadsafe

    // create a grid of GridBlocks and a grid of CartBlocks
    const bool period[3] = {false, false, false};
    const lid_t L[3] = {2, 2, 2};
    const level_t init_lvl = 0;

    // setup two ForestGrids with separate block types
    ForestGrid* grid_gridblock = CreateForestGrid(init_lvl, period, L, M_GRIDBLOCK);
    ForestGrid* grid_cartblock = CreateForestGrid(init_lvl, period, L, M_CARTBLOCK);

    // check that each one now contains blocks of the correct type
    BlockChecker cartblock_checker(typeid(CartBlock));
    cartblock_checker.CheckBlocks(grid_cartblock);
    ASSERT_EQ(cartblock_checker.block_count(), L[0]*L[1]*L[2]);
    ASSERT_EQ(cartblock_checker.types_check(), true);

    BlockChecker gridblock_checker(typeid(GridBlock));
    gridblock_checker.CheckBlocks(grid_gridblock);
    ASSERT_EQ(gridblock_checker.block_count(), L[0]*L[1]*L[2]);
    ASSERT_EQ(gridblock_checker.types_check(), true);

    // add a field to both
    Field field("field", 1);
    //grid_gridblock->AddField(&field);
    DoOpTree(nullptr, &GridBlock::AddField, grid_gridblock, &field);
    DoOpTree(nullptr, &CartBlock::AddField, grid_cartblock, &field);

    // create objects that operate on the blocks
    CartBlockOperator cartblock_operator;
    GridBlockOperator gridblock_operator;

    // Test that permissible statements execute and non-permissible ones fail
    // const function on block and qid
    DoOpMesh(&cartblock_operator, &CartBlockOperator::OperateConst, grid_cartblock, &field);
    DoOpMesh(&cartblock_operator, &CartBlockOperator::OperateConst, grid_gridblock, &field);
    DoOpMesh(&gridblock_operator, &GridBlockOperator::OperateConst, grid_gridblock, &field);
    //ASSERT_DEATH({DoOpMesh(&gridblock_operator, &GridBlockOperator::OperateConst, grid_cartblock, &field);}, ".*");

    // non-const function on block and qid
    DoOpMesh(&cartblock_operator, &CartBlockOperator::OperateNonConst, grid_cartblock, &field);
    DoOpMesh(&cartblock_operator, &CartBlockOperator::OperateNonConst, grid_gridblock, &field);
    DoOpMesh(&gridblock_operator, &GridBlockOperator::OperateNonConst, grid_gridblock, &field);
    //ASSERT_DEATH({DoOpMesh(&gridblock_operator, &GridBlockOperator::OperateNonConst, grid_cartblock, &field);}, ".*");

    // non-const block members
    std::map<std::string, Field*> empty_field_map;
    DoOpMesh(nullptr, &CartBlock::AddFields, grid_cartblock, &empty_field_map);
    DoOpMesh(nullptr, &CartBlock::AddFields, grid_gridblock, &empty_field_map);
    DoOpMesh(nullptr, &GridBlock::StatusReset, grid_gridblock);
    //ASSERT_DEATH({DoOpMesh(nullptr, &GridBlock::StatusReset, grid_cartblock);}, ".*");

    // test complete, destroy the ForestGrids
    DestroyForestGrid(grid_gridblock);
    DestroyForestGrid(grid_cartblock);

}