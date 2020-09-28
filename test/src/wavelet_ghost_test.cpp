
#include <mpi.h>

#include "doop.hpp"
#include "error.hpp"
#include "grid.hpp"
#include "gtest/gtest.h"
#include "ioh5.hpp"
#include "murphy.hpp"
#include "setvalues.hpp"
#include "subblock.hpp"
#include "wavelet.hpp"

#define DOUBLE_TOL 1e-11

using std::list;

// define a mask on the edge and corner ghosts
class MaskPhysBC : public OperatorF {
   protected:
    real_t L_[3];
    void   ApplyOpF(const qid_t* qid, GridBlock* block, Field* fid) override {
        //-------------------------------------------------------------------------
        real_t pos[3];

        for (lda_t ida = 0; ida < fid->lda(); ida++) {
            data_ptr data = block->data(fid, ida);
            for (int i2 = (-M_GS); i2 < (M_N + M_GS); i2++) {
                for (int i1 = (-M_GS); i1 < (M_N + M_GS); i1++) {
                    for (int i0 = (-M_GS); i0 < (M_N + M_GS); i0++) {
                        // get the position
                        real_t pos[3];
                        m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

                        const bool is_out_0 = (pos[0] < 0.0) || (pos[0] >= L_[0]);
                        const bool is_out_1 = (pos[1] < 0.0) || (pos[1] >= L_[1]);
                        const bool is_out_2 = (pos[2] < 0.0) || (pos[2] >= L_[2]);

                        // const bool too_big_0 = (i0 < 0) || (i0 >= M_N);
                        // const bool too_big_1 = (i1 < 0) || (i1 >= M_N);
                        // const bool too_big_2 = (i2 < 0) || (i2 >= M_N);

                        const bool to_trash = ((is_out_0 + is_out_1 + is_out_2) >= 1);

                        data[m_idx(i0, i1, i2)] = data[m_idx(i0, i1, i2)] * (!to_trash);
                    }
                }
            }
        }
        //-------------------------------------------------------------------------
    };

   public:
    MaskPhysBC(const lid_t L[3]) {
        m_begin;
        //-------------------------------------------------------------------------
        for (lda_t id = 0; id < 3; ++id) {
            L_[id] = (real_t)L[id];
        }
        //-------------------------------------------------------------------------
        m_end;
    };
};

class test_Wavelet_Ghost : public ::testing::Test {
   protected:
    void SetUp() override{};
    void TearDown() override{};
};

//==============================================================================================================================
#if (M_WAVELET_N == 2)
TEST_F(test_Wavelet_Ghost, ghost_order_2_x) {
    for (lda_t id = 0; id < 3; id++) {
        bool  period[3] = {false, false, false};
        lid_t L[3]      = {1, 1, 1};
        L[id]           = 3;
        Grid grid(0, period, L, MPI_COMM_WORLD, nullptr);

        // create the patch refinement to refine the middle tree
        real_t origin[3]      = {0.0, 0.0, 0.0};
        origin[id]            = 1.0;
        real_t      length[3] = {1.0, 1.0, 1.0};
        Patch       p1(origin, length, 1);
        list<Patch> patch{p1};
        grid.Adapt(&patch);

        Field test("test", 1);
        grid.AddField(&test);
        test.bctype(M_BC_EXTRAP_3);

        // create the initial field
        const lid_t  deg[3] = {1, 1, 1};
        const real_t dir[3] = {M_SQRT2, M_PI, -M_E};
        SetPolynom   field_init(deg, dir);
        field_init(&grid, &test);

        // pull the ghosts
        grid.GhostPull(&test);

        // IOH5 io("data_test");
        // io(&grid, &test);
        // io.dump_ghost(true);
        // io(&grid, &test);

        // create the solution field
        Field sol("sol", 1);
        grid.AddField(&sol);
        SetPolynom field_sol(deg, dir, grid.NGhostFront(), grid.NGhostBack());
        field_sol(&grid, &sol);

        // mask both the sol and the result
        MaskPhysBC mask(L);
        mask(&grid, &test);
        mask(&grid, &sol);

        // now, we need to check
        real_t          norm2, normi;
        ErrorCalculator error(&grid);
        error.Norms(&grid, &test, &sol, &norm2, &normi);

        m_log("checking in dim %d: the two norms: %e %e", id, normi, norm2);
        ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
        ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    }
}
#elif (M_WAVELET_N == 4)
TEST_F(test_Wavelet_Ghost, ghost_order_4_x) {
    for (lda_t id = 0; id < 3; id++) {
        bool  period[3] = {false, false, false};
        lid_t L[3]      = {1, 1, 1};
        L[id]           = 3;
        Grid grid(0, period, L, MPI_COMM_WORLD, nullptr);

        // create the patch refinement to refine the middle tree
        real_t origin[3]      = {0.0, 0.0, 0.0};
        origin[id]            = 1.0;
        real_t      length[3] = {1.0, 1.0, 1.0};
        Patch       p1(origin, length, 1);
        list<Patch> patch{p1};
        grid.Adapt(&patch);

        Field test("test", 1);
        grid.AddField(&test);
        test.bctype(M_BC_EXTRAP_5);

        // create the initial field
        const lid_t  deg[3] = {3, 3, 3};
        const real_t dir[3] = {M_SQRT2, M_PI, -M_E};
        SetPolynom   field_init(deg, dir);
        field_init(&grid, &test);

        // pull the ghosts
        grid.GhostPull(&test);

        // IOH5 io("data_test");
        // io(&grid, &test);
        // io.dump_ghost(true);
        // io(&grid, &test);

        // create the solution field
        Field sol("sol", 1);
        grid.AddField(&sol);
        SetPolynom field_sol(deg, dir, grid.NGhostFront(), grid.NGhostBack());
        field_sol(&grid, &sol);

        // mask both the sol and the result
        MaskPhysBC mask(L);
        mask(&grid, &test);
        mask(&grid, &sol);

        // now, we need to check
        real_t          norm2, normi;
        ErrorCalculator error(&grid);
        error.Norms(&grid, &test, &sol, &norm2, &normi);

        m_log("checking in dim %d: the two norms: %e %e", id, normi, norm2);
        ASSERT_NEAR(normi, 0.0, DOUBLE_TOL);
        ASSERT_NEAR(norm2, 0.0, DOUBLE_TOL);
    }
}
#endif
