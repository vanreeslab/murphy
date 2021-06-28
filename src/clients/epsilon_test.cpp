#include "epsilon_test.hpp"

#include "grid/grid.hpp"
#include "grid/gridcallback.hpp"
#include "operator/error.hpp"
#include "operator/setvalues.hpp"
#include "operator/xblas.hpp"
#include "tools/ioh5.hpp"

constexpr real_t sigma     = 0.05;
constexpr real_t center[3] = {0.5, 0.5, 0.5};

using std::string;
using std::to_string;

// define the initial condition and the analytical solution
static lambda_error_t lambda_error = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
    // get the position
    real_t pos[3];
    block->pos(i0, i1, i2, pos);

    // call the function
    return scalar_exp(pos, center, sigma);
};

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
    real_t depsilon = delta_eps_;
    real_t epsilon  = eps_start_;
    m_log("starting with epsilon = %e", epsilon);
    while (epsilon >= std::pow(2.0, -34)) {
        // create a grid, put a ring on it on the fixel level
        // bool period[3] = {false, false, false};
        bool  period[3]   = {true, true, true};
        lid_t grid_len[3] = {1, 1, 1};
        Grid  grid(level_start_, period, grid_len, MPI_COMM_WORLD, nullptr);
        grid.level_limit(level_min_, level_max_);

        Field scal("scalar", 1);
        grid.AddField(&scal);
        scal.bctype(M_BC_EXTRAP);

        // initial condition
        lambda_setvalue_t lambda_initcond = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);
            real_t* data = block->data(fid).Write(i0, i1, i2);
            // set value
            data[0] = scalar_exp(pos, center, sigma);
        };
        SetValue init(lambda_initcond);
        init(&grid, &scal);

        // apply it
        grid.SetTol(epsilon * 1e+20, epsilon);

        // grid.GhostPull(&scal);
        // IOH5 dump("data");
        // dump(&grid, &scal, 0);

        // compute the moment at the start
        grid.GhostPull(&scal);
        BMoment moment;
        real_t  sol_moment0;
        real_t  sol_moment1[3];
        moment(&grid, &scal, &sol_moment0, sol_moment1);

        // coarsen
        short_t count     = 1;
        level_t min_level = level_start_;
        do {
            min_level = grid.MinLevel();
            // coarsen the field if needed
            grid.Coarsen(&scal);

            level_t tmp_min_lvl = grid.MinLevel();
            level_t tmp_max_lvl = grid.MaxLevel();
            m_log("Coarsening: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

            grid.GhostPull(&scal);
            real_t coarse_moment0;
            real_t coarse_moment1[3];
            moment(&grid, &scal, &coarse_moment0, coarse_moment1);
            m_log("moments after coarsening: %e vs %e -> error = %e", sol_moment0, coarse_moment0, abs(sol_moment0 - coarse_moment0));

        } while (grid.MinLevel() < min_level && grid.MinLevel() > m_max(0, level_min_));

        // measure the moments
        grid.GhostPull(&scal);
        real_t coarse_moment0;
        real_t coarse_moment1[3];
        moment(&grid, &scal, &coarse_moment0, coarse_moment1);

        m_log("moments after coarsening: %e vs %e -> error = %e", sol_moment0, coarse_moment0, abs(sol_moment0 - coarse_moment0));

        // dump(&grid,&scal,1);

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
            min_level = grid.MinLevel();
            // force the field refinement using a patch
            grid.GhostPull(&scal);
            grid.AdaptMagic(nullptr, &patch, nullptr, &cback_StatusCheck, nullptr, &cback_UpdateDependency, nullptr);

            level_t tmp_min_lvl = grid.MinLevel();
            level_t tmp_max_lvl = grid.MaxLevel();
            m_log("Refinement: level is now %d to %d", tmp_min_lvl, tmp_max_lvl);

        } while (grid.MinLevel() < level_start_);

        // set the solution field
        Field sol("solution", 1);
        grid.AddField(&sol);
        init(&grid, &sol);
        sol.bctype(M_BC_EXTRAP);

        grid.GhostPull(&sol);
        grid.GhostPull(&scal);

        // measure the error
        real_t normi;
        Error  error;
        // error.Normi(&grid, &scal, const Field* (&sol), &normi);
        error.Normi(&grid, &scal, &lambda_error, &normi);

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
        m_log("analytical moments after refinement: %e vs %e -> error = %.12e", sol_moment1[0], moment1[0], abs(sol_moment1[0] - moment1[0]));
        m_log("analytical moments after refinement: %e vs %e -> error = %.12e", sol_moment1[1], moment1[1], abs(sol_moment1[1] - moment1[1]));
        m_log("analytical moments after refinement: %e vs %e -> error = %.12e", sol_moment1[2], moment1[2], abs(sol_moment1[2] - moment1[2]));
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
        m_log("free fields");
        grid.DeleteField(&sol);
        m_log("free fields");
        grid.DeleteField(&scal);

        // get the new epsilon
        epsilon *= depsilon;
    }
    //-------------------------------------------------------------------------
}
