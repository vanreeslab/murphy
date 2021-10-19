#include "gtest/gtest.h"

#include "src/core/types.hpp"
#include "src/grid/grid.hpp"
#include "src/operator/setvalues.hpp"

// Check that there is no nan in a field 
// Looping on the field also ensure there is no segfault. 
class GridChecker{
   public:
    GridChecker() = default;
    void operator()(ForestGrid* grid, Field* field_test) {
        if(grid->block_type() == M_GRIDBLOCK){
            DoOpMesh(this, &GridChecker::GridBlockCheck, grid, field_test); 
        }
        // Check additional information 
        // [ ... ]
    }


    void GridBlockCheck(const qid_t* qid, GridBlock* block, Field* field_test) {
        // Check the field : Look for NaN or segfault.
        const ConstMemData data_test   = block->data(field_test, 0);
        auto op = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            ASSERT_EQ(std::isnan(data_test(i0, i1, i2)), false);
        };
        bidx_t ghost_updated[2] = {field_test->get_ghost_len(0), field_test->get_ghost_len(1)};
        for_loop(&op, block->ExtendedSpan(ghost_updated));
    }
};

void SetGhostTrue(Field* field){
    const bidx_t ghost_len_actual[2] = {M_GS, M_GS};
    field->ghost_len(ghost_len_actual);
};

void AdaptAndCheck(Grid* grid, std::list<Patch>* patches, Field* field){
    grid->Adapt(patches);
    SetGhostTrue(field);

    GridChecker check;
    check(grid, field);
};


class TestPartitioner : public ::testing::TestWithParam<int> {
  public: 
    const level_t ilvl      = 1; 
    const bool    is_per[3]       = {true, true, true};
    const lid_t   domain_length[3] = {1, 1, 2};
    const real_t  patch_length[3]  = {1.0, 1.0, 2.0};
    const real_t  patch_origin[3]  = {0.0, 0.0, 0.0};
    Grid*   grid;
    Field*  test_field;

    void SetUp() override{
        
        BlockDataType blocktype = (BlockDataType) GetParam();
        grid                    = new Grid(ilvl, is_per, domain_length, blocktype, MPI_COMM_WORLD, nullptr);
        test_field              = new Field("zero_field", 1);
        grid->AddField(test_field); 
        SetGhostTrue(test_field);
    };
    void TearDown() override{
        grid->DeleteField(test_field);
        delete grid;
        delete test_field;
    };
};

// Adapt the grid from a number of blocks larger than the number of CPUs
// to a number of blocks smaller than the number of CPUs. 
TEST_P(TestPartitioner, nblocks_leq_ncpus){    
    std::list<Patch> plist;
    plist.push_back(Patch(patch_origin, patch_length, ilvl - 1));
    AdaptAndCheck(grid, &plist, test_field);

    plist.pop_back();
    plist.push_back(Patch(patch_origin, patch_length, ilvl));
    AdaptAndCheck(grid, &plist, test_field);
}

// Adapt the grid from a small number of blocks per CPUs 
// to a larger number of blocks per CPUs. 
TEST_P(TestPartitioner, nblocks_geq_ncpus){

    std::list<Patch> plist;
    plist.push_back(Patch(patch_origin, patch_length, ilvl + 1));
    AdaptAndCheck(grid, &plist, test_field);

    plist.pop_back();
    plist.push_back(Patch(patch_origin, patch_length, ilvl));
    AdaptAndCheck(grid, &plist, test_field);
}


INSTANTIATE_TEST_SUITE_P(ValidPartitioner,
                         TestPartitioner,
                         testing::Values(M_GRIDBLOCK));