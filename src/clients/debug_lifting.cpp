#include "clients/debug_lifting.hpp"

#include <list>

#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "operator/xblas.hpp"
#include "tools/ioh5.hpp"
#include "tools/patch.hpp"

class InitialCondition : public SetValue {
   protected:
    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override {
        //-------------------------------------------------------------------------
        real_t        pos[3];
        const real_t* xyz   = block->xyz();
        const real_t* hgrid = block->hgrid();

        real_t sigma     = 0.1;
        real_t center[3] = {1.0, 0.5, 0.5};

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
            const real_t rho  = rhox * rhox;  //+ rhoy * rhoy + rhoz * rhoz;

            data[m_idx(i0, i1, i2)] = fact * std::exp(-rho);

            // simply give the level
            // data[m_idx(i0, i1, i2)] = block->level();
        };

        for_loop(&op, start_, end_);
        //-------------------------------------------------------------------------
    };

   public:
    explicit InitialCondition() : SetValue(nullptr){};
};

void DebugLifting::InitParam(ParserArguments* param) {}

void DebugLifting::Run() {
    //-------------------------------------------------------------------------
    // create a grid, put a ring on it on the fixel level
    // bool period[3] = {false, false, false};
    bool  period[3]   = {false, true, true};
    lid_t grid_len[3] = {2, 1, 1};
    Grid  grid(1, period, grid_len, MPI_COMM_WORLD, nullptr);
    grid.level_limit(0, 1);

    // add the field
    Field scal("scalar", 1);
    grid.AddField(&scal);
    scal.bctype(M_BC_EXTRAP);

    // set the initial condition
    InitialCondition ring;
    ring(&grid, &scal);

    // dump
    IOH5 dump("data");
    grid.GhostPull(&scal);
    dump(&grid, &scal, 0);

    // get the moment:
    BMoment moment;
    real_t  fine_moment0, fine_moment1[3];
    moment(&grid, &scal, &fine_moment0, fine_moment1);

    // go coarse
    {
        std::list<Patch> plist;

        real_t p1_o[3] = {0.0, 0.0, 0.0};
        real_t p1_l[3] = {1.0, 1.0, 1.0};
        plist.push_back(Patch(p1_o, p1_l, 0));
        // real_t p2_o[3] = {0.0, 0.0, 0.0};
        // real_t p2_l[3] = {1.0, 1.0, 1.0};
        // plist.push_back(Patch(p2_o, p2_l, 0));

        // adapt the grid
        grid.Adapt(&plist);

        // dump again, the ghosts should be ok
        grid.GhostPull(&scal);
        dump(&grid, &scal, 1);
        dump.dump_ghost(true);
        dump(&grid, &scal, 1);
    }
    real_t coarse_moment0, coarse_moment1[3];
    moment(&grid, &scal, &coarse_moment0, coarse_moment1);

    // move back up
    {
        std::list<Patch> plist;

        real_t p1_o[3] = {0.0, 0.0, 0.0};
        real_t p1_l[3] = {1.0, 1.0, 1.0};
        plist.push_back(Patch(p1_o, p1_l, 1));
        // real_t p2_o[3] = {0.0, 0.0, 0.0};
        // real_t p2_l[3] = {1.0, 1.0, 1.0};
        // plist.push_back(Patch(p2_o, p2_l, 0));

        // adapt the grid
        grid.Adapt(&plist);

        // dump again, the ghosts should be ok
        grid.GhostPull(&scal);
        dump(&grid, &scal, 2);
        dump.dump_ghost(true);
        dump(&grid, &scal, 2);
    }
    // get the moment:
    // BMoment moment;
    real_t adapt_moment0, adapt_moment1[3];
    moment(&grid, &scal, &adapt_moment0, adapt_moment1);

    m_log("moments: %e -> %e : error = %e", fine_moment0, adapt_moment0, fabs(fine_moment0 - adapt_moment0));

    //-------------------------------------------------------------------------
}