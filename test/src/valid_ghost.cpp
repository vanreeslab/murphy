#include <cmath>
#include <list>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "gtest/gtest.h"
#include "operator/advection.hpp"
#include "operator/error.hpp"
#include "tools/ioh5.hpp"

#define DOUBLE_TOL 1e-13

class InitialCondition : public SetValue {
   protected:
    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override {
        //-------------------------------------------------------------------------
        // get the pointers correct
        real_t* data = block->data(fid, 0).Write();

        auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            // simply give the level
            data[m_idx(i0, i1, i2)] = pos[0] + pos[1] + pos[2];
        };

        // m_log("filling from %d to %d", start_, end_);
        for_loop(&op, start_, end_);
        //-------------------------------------------------------------------------
    };

   public:
    explicit InitialCondition() : SetValue(nullptr){};
    explicit InitialCondition(m_ptr<const Wavelet> interp) : SetValue(interp){};
};

class valid_Ghost : public ::testing::Test {
    void SetUp() override{};
    void TearDown() override{};
};

using std::list;
using std::string;

/**
 * @brief validate every possible intersection between blocks
 * 
 * for one intersection, there is max 8 blocks sharing the same corner.
 * Every block might be coarse or fine, so in total 2^8 possibilities
 * 
 */
TEST_F(valid_Ghost, valid_ghost_extrap) {
    m_assert(sizeof(int) >= 1, "the int must be of size > 8");
    for (int case_id = 0; case_id < 256; ++case_id) {
        m_log("case id: %d", case_id);
        // every bit corresponds to the status of a block: 1 is fine, 0 is coarse
        list<Patch> patch_list;
        for (int itree = 0; itree < 8; ++itree) {
            real_t origin[3] = {1.0 * (itree % 2),
                                1.0 * (itree % 4) / 2,
                                1.0 * (itree / 4)};
            real_t length[3] = {1.0, 1.0, 1.0};

            level_t level = (case_id >> itree) % 2;
            Patch   patch(origin, length, level);
            patch_list.push_back(patch);

            m_log("tree %d @ %f %f %f has level %d", itree, origin[0], origin[1], origin[2], level);
        }

        // create a grid
        bool period[3] = {false, false, false};
        // bool  period[3]   = {true, true, true};
        lid_t grid_len[3] = {2, 2, 2};
        Grid  grid(0, period, grid_len, MPI_COMM_WORLD, nullptr);
        grid.level_limit(0, 1);

        // adapt given the patch
        grid.Adapt(&patch_list);

        // create a field an put it on it
        Field scal("scal", 1);
        grid.AddField(&scal);
        InitialCondition init;
        init(&grid, &scal);

        // get the Ghosts:
        scal.bctype(M_BC_EXTRAP);
        grid.GhostPull(&scal);

        // analytical solution
        Field sol("sol", 1);
        grid.AddField(&sol);
        InitialCondition initsol(grid.interp());
        initsol(&grid, &sol);

        Field err("error", 1);
        grid.AddField(&err);

        // check the error, ghost included!
        real_t          err2, erri;
        ErrorCalculator error(grid.interp());
        error.Norms(&grid, &scal, &sol, &err, &err2, &erri);

        m_log("errors: err2 = %e erri = %e", err2, erri);

        ASSERT_DOUBLE_EQ(erri,0.0);

        // // print
        // IOH5 dump("data");
        // dump(&grid, &scal);
        // dump.dump_ghost(true);
        // dump(&grid, &scal);
        // dump.dump_ghost(true);
        // dump(&grid, &sol);
        // dump.dump_ghost(true);
        // dump(&grid, &err);
    }
}
