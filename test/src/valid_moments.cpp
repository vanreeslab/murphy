#include <cmath>
#include <list>
#include <limits>

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "gtest/gtest.h"
#include "operator/setvalues.hpp"
#include "operator/xblas.hpp"
#include "tools/ioh5.hpp"

class InitCondition : public SetValue {
   protected:
    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override {
        //-------------------------------------------------------------------------
        real_t        pos[3];
        const real_t* xyz   = block->xyz();
        const real_t* hgrid = block->hgrid();

        real_t sigma     = 0.1;
        real_t center[3] = {1.0, 1.0, 1.0};

        const real_t oo_sigma2 = 1.0 / (sigma * sigma);
        const real_t fact      = 1.0;  /// sqrt(M_PI * sigma * sigma);  //todo change that because sqrt(M_PI * sigma_ * sigma_) is the initial amplitude

        // get the pointers correct
        real_t* data = block->data(fid, 0).Write();

        auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            // compute the gaussian
            const real_t rhox = (pos[0] - center[0]) / sigma;
            const real_t rhoy = (pos[1] - center[1]) / sigma;
            const real_t rhoz = (pos[2] - center[2]) / sigma;
            const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

            data[m_idx(i0, i1, i2)] = fact * std::exp(-rho);

            // simply give the level
            // data[m_idx(i0, i1, i2)] = block->level();
        };

        for_loop(&op, start_, end_);
        //-------------------------------------------------------------------------
    };

   public:
    explicit InitCondition() : SetValue(nullptr){};
};

class valid_Moments : public ::testing::Test {
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
TEST_F(valid_Moments, valid_moment_periodic) {
    real_t eps        = 1000.0 * std::numeric_limits<real_t>::epsilon();

    m_assert(sizeof(int) >= 1, "the int must be of size > 8");
    for (int case_id = 0; case_id < 256; ++case_id) {
        m_log("--------------------------------------------------------------------------------");
        m_log("case id: %d", case_id);
        m_log("--------------");

        // create a grid - uniform on level 1
        // bool  period[3]   = {true, true, true};
        bool  period[3]   = {false, true, true};
        lid_t grid_len[3] = {2, 2, 2};
        Grid  grid(1, period, grid_len, MPI_COMM_WORLD, nullptr);
        grid.level_limit(0, 1);

        // create a field an put it on it
        Field scal("scal", 1);
        grid.AddField(&scal);
        InitCondition init;
        init(&grid, &scal);

        // get the Ghosts:
        scal.bctype(M_BC_EXTRAP);
        grid.GhostPull(&scal);

        // dump
        // IOH5 dump("data");
        // dump(&grid, &scal, 0);

        // get the moment:
        BMoment moment;
        real_t  fine_moment0, fine_moment1[3];
        moment(&grid, &scal, &fine_moment0, fine_moment1);

        // adapt the tree
        {
            // create the patch list
            list<Patch> patch_list;
            for (int itree = 0; itree < 8; ++itree) {
                real_t origin[3] = {1.0 * (itree % 2),
                                    1.0 * ((itree % 4) / 2),
                                    // 1.0 * (itree % 4) / 2,
                                    1.0 * (itree / 4)};
                real_t length[3] = {1.0, 1.0, 1.0};

                level_t level = (case_id >> itree) % 2;
                patch_list.push_back(Patch(origin, length, level));

                m_log("tree %d @ %f %f %f has level %d", itree, origin[0], origin[1], origin[2], level);
            }
            // adapt
            grid.Adapt(&patch_list);
            grid.GhostPull(&scal);
        }
        // dump(&grid, &scal, 1);

        real_t coarse_moment0, coarse_moment1[3];
        moment(&grid, &scal, &coarse_moment0, coarse_moment1);
        real_t mom0_coarse_error = fabs(fine_moment0 - coarse_moment0);
        m_log("[case %d] coarse moment error = |%e - %e| = %e", case_id, fine_moment0, coarse_moment0, mom0_coarse_error);
        ASSERT_LT(mom0_coarse_error, eps);
        // go back up
        {
            // create the patch list
            list<Patch> patch_list;
            for (int itree = 0; itree < 8; ++itree) {
                real_t origin[3] = {1.0 * (itree % 2),
                                    1.0 * ((itree % 4) / 2),
                                    // 1.0 * (itree % 4) / 2,
                                    1.0 * (itree / 4)};
                real_t length[3] = {1.0, 1.0, 1.0};

                level_t level = 1;
                patch_list.push_back(Patch(origin, length, level));

                m_log("tree %d @ %f %f %f has level %d", itree, origin[0], origin[1], origin[2], level);
            }
            // adapt
            grid.Adapt(&patch_list);
            grid.GhostPull(&scal);
        }
        // dump(&grid, &scal, 2);

        // get the moment:
        // BMoment moment;
        real_t smooth_moment0, smooth_moment1[3];
        moment(&grid, &scal, &smooth_moment0, smooth_moment1);
        real_t mom0_smooth_error = fabs(fine_moment0 - smooth_moment0);
        m_log("[case %d] smoothed moment error = |%e - %e| = %e < %e", case_id, fine_moment0, smooth_moment0, mom0_smooth_error, eps);
        ASSERT_LT(mom0_smooth_error, eps);

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
