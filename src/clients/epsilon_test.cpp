#include "epsilon_test.hpp"

#include "grid/grid.hpp"
#include "operator/setvalues.hpp"
#include "operator/error.hpp"
#include "grid/gridcallback.hpp"
#include "xblas.hpp"

using std::string;
using std::to_string;

void EpsilonTest::InitParam(ParserArguments* param){
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // get the level
    level_start_ = param->init_lvl;

    //-------------------------------------------------------------------------
}

void EpsilonTest::Run(){
    //-------------------------------------------------------------------------
    real_t depsilon = 0.1;
    real_t epsilon = epsilon_start_;
    while(epsilon >= 1e-4){
        // create a grid, put a ring on it on the fixel level
        bool   period[3]   = {true, true, true};
        lid_t grid_len[3] = {1, 1, 1};
        Grid grid(level_start_, period, grid_len, MPI_COMM_WORLD, nullptr);

        Field scal("scalar", 1);
        grid.AddField(&scal);

        real_t center [3] = {0.5 ,0.5, 0.5};
        lda_t normal = 2;
        real_t sigma = 0.025;
        real_t radius = 0.5;
        real_t velocity [3] = {0.0,0.0,0.0};
        SetScalarRing ring(normal, center, sigma, radius, velocity);
        ring.SetTime(0.0);
        ring(&grid, &scal);

        grid.SetTol(epsilon * 1e+20, epsilon);

        // Field details("details", 1);
        // grid.DumpDetails(&scal,&details);
        grid.GhostPull(&scal);
        BMoment moment;
        real_t sol_moment0, sol_moment1[3];
        moment(&grid, &scal, &sol_moment0, sol_moment1);

        // coarsen 
        level_t min_level = level_start_;
        do {
            min_level = grid.MinLevel();
            // coarsen the field if needed
            grid.Coarsen(&scal);

            m_log("Coarsening: level is now %d to %d",grid.MinLevel(),grid.MaxLevel());

        } while (grid.MinLevel() < min_level);
        // measure the moments
        grid.GhostPull(&scal);
        real_t coarse_moment0, coarse_moment1[3];
        moment(&grid, &scal, &coarse_moment0, coarse_moment1);

        // and go back up
        std::list<Patch> patch;
        real_t      origin[3] = {0.0, 0.0, 0.0};
        real_t      length[3] = {grid_len[0], grid_len[1], grid_len[2]};
        patch.push_back(Patch(origin, length, level_start_));
        do {
            min_level = grid.MinLevel();
            // force the field refinement using a patch
            grid.GhostPull(&scal);
            grid.Adapt(nullptr, nullptr, &cback_Patch, reinterpret_cast<void*>(&patch), cback_UpdateDependency, nullptr);

            m_log("Refinement: level is now %d to %d",grid.MinLevel(),grid.MaxLevel());

        } while (grid.MinLevel() < level_start_);

        // set the solution field
        Field sol("solution", 1);
        grid.AddField(&sol);
        ring(&grid, &sol);


        grid.GhostPull(&sol);
        grid.GhostPull(&scal);

        // measure the error
        real_t          normi;
        ErrorCalculator error;
        error.Normi(&grid, &scal, &sol, &normi);

        // measure the moments
        real_t moment0, moment1[3];
        moment(&grid, &scal, &moment0, moment1);

        m_log("moment 0: coarse %e and fine %e vs %e", coarse_moment0, moment0, sol_moment0);

        m_log("epsilon = %e, error = %e", epsilon, normi);
        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0) {
            FILE* file = fopen(std::string("data/moments_w" + std::to_string(M_WAVELET_N) + std::to_string(M_WAVELET_NT) + ".data").c_str(), "a+");
            fprintf(file, "%e;%16.16e", epsilon, normi);
            // fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e", coarse_moment0,coarse_moment1[0],coarse_moment1[1],coarse_moment1[2]);
            fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e", moment0,moment1[0],moment1[1],moment1[2]);
            fprintf(file, ";%16.16e;%16.16e;%16.16e;%16.16e", sol_moment0,sol_moment1[0],sol_moment1[1],sol_moment1[2]);
            fprintf(file,"\n");
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