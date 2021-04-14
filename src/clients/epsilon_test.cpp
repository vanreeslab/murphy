#include "epsilon_test.hpp"

#include "grid/grid.hpp"
#include "grid/gridcallback.hpp"
#include "operator/error.hpp"
#include "operator/setvalues.hpp"
#include "operator/xblas.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

class InitialCondition : public SetValue {
   protected:
    void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override {
        //-------------------------------------------------------------------------
        real_t        pos[3];
        const real_t* xyz   = block->xyz();
        const real_t* hgrid = block->hgrid();

        real_t sigma     = 0.05;
        real_t center[3] = {0.5, 0.5, 0.5};

        const real_t oo_sigma2 = 1.0 / (sigma * sigma);
        const real_t fact      = 1.0/ sqrt(M_PI * sigma * sigma); //todo change that because sqrt(M_PI * sigma_ * sigma_) is the initial amplitude

        // get the pointers correct
        real_t* data = block->data(fid, 0).Write();

        auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the position
            real_t pos[3];
            m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

            // compute the gaussian
            const real_t rhox = (pos[0] - center[0]) / sigma;
            const real_t rhoy = (pos[1] - center[1]) / sigma;
            const real_t rhoz = (pos[2] - center[2]) / sigma;
            const real_t rho  = rhox * rhox + rhoy * rhoy + rhoz * rhoz;

            // data[m_idx(i0, i1, i2)] =   fact * std::exp(-rho);
            data[m_idx(i0, i1, i2)] = fact * std::exp(-pow(rhox,2) -pow(rhoy,2) );//-pow(rhoz,2));

            // const real_t rhox_1     = (pos[1] - center[1] - 0.25) / sigma;
            // const real_t rhox_2     = (pos[1] - center[1] + 0.25) / sigma;
            // data[m_idx(i0, i1, i2)] = std::exp(-pow(rhox_1, 2)) + std::exp(-pow(rhox_2, 2));
        };

        for_loop(&op, start_, end_);
        //-------------------------------------------------------------------------
    };

   public:
    InitialCondition() : SetValue(nullptr){};
};

// class CompactInitialCondition : public SetValue {
//    protected:
//     void FillGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid) override {
//         //-------------------------------------------------------------------------
//         real_t        pos[3];
//         const real_t* xyz   = block->xyz();
//         const real_t* hgrid = block->hgrid();

//         real_t R         = 0.4;
//         real_t sigma     = 0.3;
//         real_t center[3] = {0.5, 0.5, 0.5};

//         const real_t oo_sigma2 = 1.0 / (sigma * sigma);
//         const real_t oo_R2     = 1.0 / (R * R);
//         const real_t fact      = 1.0;  /// sqrt(M_PI * sigma_ * sigma_); //todo change that because sqrt(M_PI * sigma_ * sigma_) is the initial amplitude

//         // get the pointers correct
//         real_t* data = block->data(fid, 0).Write();

//         auto op = [=, &data](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // get the position
//             real_t pos[3];
//             m_pos(pos, i0, i1, i2, block->hgrid(), block->xyz());

//             // compute the gaussian
//             const real_t rx   = (pos[0] - center[0]);
//             const real_t ry   = (pos[1] - center[1]);
//             const real_t rz   = (pos[2] - center[2]);
//             const real_t r2   = (rx * rx + ry * ry + rz * rz);
//             const real_t rho2 = r2 * oo_sigma2;
//             const real_t Rho2 = r2 * oo_R2;

//             if (Rho2 < 1.0) {
//                 data[m_idx(i0, i1, i2)] = fact * exp(-rho2 / (1.0 - Rho2));
//             } else {
//                 data[m_idx(i0, i1, i2)] = 0.0;
//             }
//         };

//         for_loop(&op, start_, end_);
//         //-------------------------------------------------------------------------
//     };

//    public:
//     CompactInitialCondition() : SetValue(nullptr){};
// };

void EpsilonTest::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // get the level
    level_start_ = param->init_lvl;
    level_min_   = param->level_min;
    level_max_   = param->level_max;

    eps_start_ = param->eps_start;
    delta_eps_ = param->delta_eps;

    //-------------------------------------------------------------------------
}

void EpsilonTest::Run() {
    //-------------------------------------------------------------------------
    real_t depsilon = 1e-10;
    real_t epsilon  = 1e-2;  //epsilon_start_;
    while (epsilon >= std::pow(2.0, -26)) {
        // create a grid, put a ring on it on the fixel level
        // bool period[3] = {false, false, false};
        bool  period[3]   = {true, true, true};
        lid_t grid_len[3] = {1, 1, 1};
        Grid  grid(level_start_, period, grid_len, MPI_COMM_WORLD, nullptr);
        grid.level_limit(level_min_, level_max_);

        Field scal("scalar", 1);
        grid.AddField(&scal);
        scal.bctype(M_BC_EXTRAP);

        real_t offset    = 0.0;  //1.0 / M_PI * 1.0/(std::pow(2, level_start_) * M_N);  // this is a fraction of h
        real_t center[3] = {0.5 + offset, 0.5 + offset, 0.5 + offset};

        // ring
        lda_t         normal      = 2;
        real_t        sigma       = 0.05;
        real_t        radius      = 0.25;
        real_t        velocity[3] = {0.0, 0.0, 0.0};
        SetScalarRing ring(normal, center, sigma, radius, velocity);

        // exponential
        // real_t sigma     = 0.05;
        // real_t sigmav[3] = {sigma,sigma,sigma};
        // SetExponential ring(center, sigmav, 1.0);
        // ring.SetTime(0.0);

        // // testing polynomial
        // const lid_t  deg[3]   = {5, 2, 3};
        // const real_t dir[3]   = {M_PI, 0.0, 0.0};
        // const real_t shift[3] = {0.0, 0.23, 0.46};
        // SetPolynom   ring(deg, dir, shift);

        // custon stuffs
        // InitialCondition ring;
        // CompactInitialCondition ring;

        // apply it
        ring(&grid, &scal);
        grid.SetTol(epsilon * 1e+20, epsilon);

        grid.GhostPull(&scal);
        IOH5 dump("data");
        dump(&grid, &scal, 0);

        // compute the moment at the start
        grid.GhostPull(&scal);
        BMoment moment;
        real_t  sol_moment0, sol_moment1[3];
        moment(&grid, &scal, &sol_moment0, sol_moment1);
        // BDiscreteMoment dmoment;
        // real_t          sol_dmoment0, sol_dmoment1[3];
        // dmoment(&grid, &scal, &sol_dmoment0, sol_dmoment1);

        // coarsen
        short_t count = 1;
        level_t min_level = level_start_;
        do {
            min_level = grid.MinLevel();
            // coarsen the field if needed
            grid.Coarsen(&scal);

            level_t tmp_min_lvl = grid.MinLevel();
            level_t tmp_max_lvl = grid.MaxLevel();
            m_log("Coarsening: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

            grid.GhostPull(&scal);
            real_t coarse_moment0, coarse_moment1[3];
            moment(&grid, &scal, &coarse_moment0, coarse_moment1);
            m_log("moments after coarsening: %e vs %e -> error = %e", sol_moment0, coarse_moment0, abs(sol_moment0 - coarse_moment0));


            // {
            //     // set the solution field
            //     Field sol("solution", 1);
            //     grid.AddField(&sol);
            //     ring(&grid, &sol);
            //     sol.bctype(M_BC_EXTRAP);

            //     grid.GhostPull(&sol);
            //     grid.GhostPull(&scal);

            //     // measure the error
            //     real_t          normi;
            //     ErrorCalculator error;
            //     error.Normi(&grid, &scal, &sol, &normi);
            //     m_log("epsilon = %e, error = %.12e", epsilon, normi);

            //     grid.DeleteField(&sol);
            // }

            // dump(&grid,&scal,count);
            // count++;

            // grid.Coarsen(&scal);
            // grid.GhostPull(&scal);
            // // real_t coarse_moment0, coarse_moment1[3];
            // moment(&grid, &scal, &coarse_moment0, coarse_moment1);
            // m_log("moments after coarsening: %e vs %e -> error = %e", sol_moment0, coarse_moment0, abs(sol_moment0 - coarse_moment0));

            // dump(&grid,&scal,count);

            // set the solution field
            // Field sol("solution", 1);
            // grid.AddField(&sol);
            // ring(&grid, &sol);
            // sol.bctype(M_BC_EXTRAP);

            // grid.GhostPull(&sol);
            // grid.GhostPull(&scal);

            // // measure the error
            // real_t          normi;
            // ErrorCalculator error;
            // error.Normi(&grid, &scal, &sol, &normi);
            // m_log("epsilon = %e, error = %.12e", epsilon, normi);

        } while (grid.MinLevel() < min_level && grid.MinLevel() > m_max(0, level_min_));
        // } while (grid.MinLevel() < min_level && grid.MinLevel() > (level_start_ - 2));

        // measure the moments
        grid.GhostPull(&scal);
        real_t coarse_moment0, coarse_moment1[3];
        moment(&grid, &scal, &coarse_moment0, coarse_moment1);
        // real_t coarse_dmoment0, coarse_dmoment1[3];
        // dmoment(&grid, &scal, &coarse_dmoment0, coarse_dmoment1);

        m_log("moments after coarsening: %e vs %e -> error = %e", sol_moment0, coarse_moment0, abs(sol_moment0 - coarse_moment0));

        dump(&grid,&scal,1);

        // track the number of block, levels
        level_t grid_level_min = grid.MinLevel();
        level_t grid_level_max = grid.MaxLevel();
        long    nblock         = grid.global_num_quadrants();

        // and go back up
        std::list<Patch> patch;
        real_t           origin[3] = {0.0, 0.0, 0.0};
        real_t           length[3] = {(real_t)grid_len[0], (real_t)grid_len[1], (real_t)grid_len[2]};
        patch.push_back(Patch(origin, length, level_start_));
        do {

            // {
            //     // set the solution field
            //     Field sol("solution", 1);
            //     grid.AddField(&sol);
            //     ring(&grid, &sol);
            //     sol.bctype(M_BC_EXTRAP);

            //     grid.GhostPull(&sol);
            //     grid.GhostPull(&scal);

            //     // measure the error
            //     real_t          normi;
            //     ErrorCalculator error;
            //     error.Normi(&grid, &scal, &sol, &normi);
            //     m_log("epsilon = %e, error = %.12e", epsilon, normi);

            //     grid.DeleteField(&sol);
            // }


            min_level = grid.MinLevel();
            // force the field refinement using a patch
            grid.GhostPull(&scal);
            m_log("need to adapt the call to Adapt Magic");
            // grid.AdaptMagic(nullptr, nullptr, &cback_Patch, reinterpret_cast<void*>(&patch), cback_UpdateDependency, nullptr);

            level_t tmp_min_lvl = grid.MinLevel();
            level_t tmp_max_lvl = grid.MaxLevel();
            m_log("Refinement: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

        } while (grid.MinLevel() < level_start_);

        // set the solution field
        Field sol("solution", 1);
        grid.AddField(&sol);
        ring(&grid, &sol);
        sol.bctype(M_BC_EXTRAP);

        grid.GhostPull(&sol);
        grid.GhostPull(&scal);

        // measure the error
        real_t          normi;
        ErrorCalculator error;
        error.Normi(&grid, &scal, &sol, &normi);

        grid.DumpLevels(1, "data", string("w" + std::to_string(M_WAVELET_N) + std::to_string(M_WAVELET_NT)));

        // real_t norm2;
        // Field err("error",1);
        // grid.AddField(&err);
        // err.bctype(M_BC_EXTRAP);
        // error.Norms(&grid,&scal,&sol,&err,&norm2,&normi);
        // grid.GhostPull(&err);

        // measure the moments
        real_t moment0, moment1[3];
        moment(&grid, &scal, &moment0, moment1);
        real_t dmoment0, dmoment1[3];
        // dmoment(&grid, &scal, &dmoment0, dmoment1);
        m_log("analytical moments after refinement: %e vs %e -> error = %.12e", sol_moment0, moment0, abs(sol_moment0 - moment0));
        // m_log("analytical moments after refinement: %e vs %e -> error = %e", sol_dmoment0, dmoment0, abs(sol_dmoment0 - dmoment0));

        // m_log("moment 0: from %e to %e: error %e", sol_moment0, moment0, fabs(sol_moment0 - moment0));
        // m_log("moment 1x: from %e to %e: error %e", sol_moment1[0], moment1[0], fabs(sol_moment1[0] - moment1[0]));
        // m_log("moment 1y: from %e to %e: error %e", sol_moment1[1], moment1[1], fabs(sol_moment1[1] - moment1[1]));
        // m_log("moment 1z: from %e to %e: error %e", sol_moment1[2], moment1[2], fabs(sol_moment1[2] - moment1[2]));
        // m_log("moment 0: from %e to %e: error %e", sol_moment0, coarse_moment0, fabs(sol_moment0 - coarse_moment0));
        // m_log("moment 1x: from %e to %e: error %e", sol_moment1[0], coarse_moment1[0], fabs(sol_moment1[0] - coarse_moment1[0]));
        // m_log("moment 1y: from %e to %e: error %e", sol_moment1[1], coarse_moment1[1], fabs(sol_moment1[1] - coarse_moment1[1]));
        // m_log("moment 1z: from %e to %e: error %e", sol_moment1[2], coarse_moment1[2], fabs(sol_moment1[2] - coarse_moment1[2]));
        // m_log("discrete moment 0: from %e to %e: error %e", sol_dmoment0, coarse_dmoment0, fabs(sol_dmoment0 - coarse_dmoment0));
        // m_log("discrete moment 1x: from %e to %e: error %e", sol_dmoment1[0], coarse_dmoment1[0], fabs(sol_dmoment1[0] - coarse_dmoment1[0]));
        // m_log("discrete moment 1y: from %e to %e: error %e", sol_dmoment1[1], coarse_dmoment1[1], fabs(sol_dmoment1[1] - coarse_dmoment1[1]));
        // m_log("discrete moment 1z: from %e to %e: error %e", sol_dmoment1[2], coarse_dmoment1[2], fabs(sol_dmoment1[2] - coarse_dmoment1[2]));

        m_log("epsilon = %e, error = %.12e", epsilon, normi);
        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            FILE* file = fopen(std::string("data/moments_w" + std::to_string(M_WAVELET_N) + std::to_string(M_WAVELET_NT) + ".data").c_str(), "a+");
            fprintf(file, "%e;%ld;%d;%d;%16.16e", epsilon, nblock, grid_level_min, grid_level_max, normi);

            fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e",
                    std::fabs(coarse_moment0 - sol_moment0),
                    std::fabs(coarse_moment1[0] - sol_moment1[0]),
                    std::fabs(coarse_moment1[1] - sol_moment1[1]),
                    std::fabs(coarse_moment1[2] - sol_moment1[2]));
            fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e",
                    std::fabs(moment0 - sol_moment0),
                    std::fabs(moment1[0] - sol_moment1[0]),
                    std::fabs(moment1[1] - sol_moment1[1]),
                    std::fabs(moment1[2] - sol_moment1[2]));
            // fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e",
            //         std::fabs(coarse_dmoment0 - sol_dmoment0),
            //         std::fabs(coarse_dmoment1[0] - sol_dmoment1[0]),
            //         std::fabs(coarse_dmoment1[1] - sol_dmoment1[1]),
            //         std::fabs(coarse_dmoment1[2] - sol_dmoment1[2]));
            // fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e",
            //         std::fabs(dmoment0 - sol_dmoment0),
            //         std::fabs(dmoment1[0] - sol_dmoment1[0]),
            //         std::fabs(dmoment1[1] - sol_dmoment1[1]),
            //         std::fabs(dmoment1[2] - sol_dmoment1[2]));
            fprintf(file, "\n");
            fclose(file);
        }

        // cleanup the fields
        grid.DeleteField(&sol);
        grid.DeleteField(&scal);

        // get the new epsilon
        epsilon *= depsilon;
    }
    //-------------------------------------------------------------------------
}
