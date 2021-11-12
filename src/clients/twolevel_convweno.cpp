#include "twolevel_convweno.hpp"


#include "clients/convergence_weno.hpp"

#include <list>
#include <string>

#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "operator/advection.hpp"
#include "operator/error.hpp"
#include "operator/xblas.hpp"

using std::list;
using std::string;
using std::to_string;

TwoLevelConvWeno::~TwoLevelConvWeno() {
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
}

void TwoLevelConvWeno::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    ilevel_   = m_max(param->init_lvl, 0);
    adapt_    = !param->no_adapt;
    fix_weno_ = param->fix_weno;

    // get the level
    level_min_   = param->level_min;
    level_max_   = param->level_max;

    nu_ = (param->reynolds < 0.0) ? 0.0 : ((1.0) / param->reynolds);
    //-------------------------------------------------------------------------
}

// static const real_t sigma     = 0.05;
// static const real_t radius    = 0.25;
// static const real_t beta      = 3.0;
// static const auto   freq      = std::vector<short_t>{};  //std::vector<short_t>{5, 1000};
// static const auto   amp       = std::vector<real_t>{};   //std::vector<real_t>{0.2, 0.2};
// static const real_t center[3] = {1.0,1.0,};
static const bidx_t grid_len[3] = {2, 2, 2};
// lambdas

void TwoLevelConvWeno::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // get a random velocity
    real_t rand_vel[3] = {0.0, 0.0, 0.0};
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        rand_vel[0] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[1] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
        rand_vel[2] = -1.0 + (static_cast<real_t>(std::rand()) / static_cast<real_t>(RAND_MAX)) * 2.0;
    }
    MPI_Bcast(rand_vel, 3, M_MPI_REAL, 0, MPI_COMM_WORLD);

    static lambda_setvalue_t lambda_initcond = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        // call the function
        block->data(fid, 0)(i0, i1, i2) = sin(2.0 * M_PI * 2.0 / ((real_t)grid_len[0]) * pos[0]) +
                                          sin(2.0 * M_PI * 2.0 / ((real_t)grid_len[1]) * pos[1]) +
                                          sin(2.0 * M_PI * 2.0 / ((real_t)grid_len[2]) * pos[2]);
    };
    static lambda_error_t lambda_sol_rhs = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        // call the function
        return (-rand_vel[0]) * (4.0 * M_PI / ((real_t)grid_len[0])) * cos(2.0 * M_PI * 2.0 / ((real_t)grid_len[0]) * pos[0]) +
               (-rand_vel[1]) * (4.0 * M_PI / ((real_t)grid_len[1])) * cos(2.0 * M_PI * 2.0 / ((real_t)grid_len[1]) * pos[1]) +
               (-rand_vel[2]) * (4.0 * M_PI / ((real_t)grid_len[2])) * cos(2.0 * M_PI * 2.0 / ((real_t)grid_len[2]) * pos[2]) +
                (-nu_) *(
                    pow(4.0 * M_PI / ((real_t)grid_len[0]),2)* sin(4.0 * M_PI / ((real_t)grid_len[0]) * pos[0])+
                    pow(4.0 * M_PI / ((real_t)grid_len[1]),2)* sin(4.0 * M_PI / ((real_t)grid_len[1]) * pos[1])+
                    pow(4.0 * M_PI / ((real_t)grid_len[2]),2)* sin(4.0 * M_PI / ((real_t)grid_len[2]) * pos[2]));
    };

    //......................................................................
    for(level_t il=level_min_; il<= (level_max_-1); ++il){
        m_log("================================================================================");
        m_log("levels = %d and %d", il, il + 1);
        //......................................................................
        // create a uniform grid at level L+1
        bool period[3] = {true, true, true};
        Grid grid(il + 1, period, grid_len, M_GRIDBLOCK, MPI_COMM_WORLD, nullptr);
        // // adapt
        // if (adapt_) {
        //     list<Patch> patch_list;
        //     // upgrade one of the block
        //     real_t origin[3] = {0.0, 0.0, 0.0};
        //     real_t length[3] = {1.0, 1.0, 1.0};
        //     patch_list.push_back(Patch(origin, length, il + 1));
        //     grid.Adapt(&patch_list);
        // }
        //......................................................................
        // create the scalar field and adapt if needed
        Field test("scalar", 1);
        grid.AddField(&test);

        SetValue init(lambda_initcond);
        init(&grid, &test);

        // need to ghost
        const bidx_t ghost_len[2] = {m_max(3, grid.interp()->nghost_front()),
                                     m_max(3, grid.interp()->nghost_back())};
        grid.GhostPull(&test, ghost_len);

        // compute the moments
        real_t  moment0    = 0.0;
        real_t  moment1[3] = {0.0, 0.0, 0.0};
        BMoment moment;
        moment(&grid, &test, &moment0, moment1);
        m_log("moments are %e",moment0);

        // we need to adapt the grid, coarsen one block
        if (adapt_) {
            list<Patch> patch_list;
            // upgrade one of the block
            real_t origin[3] = {0.0, 0.0, 0.0};
            real_t length[3] = {1.0, 1.0, 1.0};
            patch_list.push_back(Patch(origin, length, il));
            grid.Adapt(&patch_list);
        }

        grid.GhostPull(&test, ghost_len);

        real_t  adapted_moment0    = 0.0;
        real_t  adapted_moment1[3] = {0.0, 0.0, 0.0};
        moment(&grid, &test, &adapted_moment0, adapted_moment1);
        m_log("error moment after adaptation = %e",adapted_moment0-moment0);


        // velocity
        Field vel("vel", 3);
        grid.AddField(&vel);
        lambda_setvalue_t lambda_vel = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
            block->data(fid, 0)(i0, i1, i2) = rand_vel[0];
            block->data(fid, 1)(i0, i1, i2) = rand_vel[1];
            block->data(fid, 2)(i0, i1, i2) = rand_vel[2];
        };
        // const bidx_t ghost_len[2] = {3, 3};
        SetValue     vel_init(lambda_vel, ghost_len);
        vel_init(&grid, &vel);
        
        //......................................................................
        // compute the derivative
        Field dtest("deriv", 1);
        grid.AddField(&dtest);

        // get the analytical solution
        const lda_t    normal     = 2;
        

        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const real_t  hmax            = grid.FinestH();
        const real_t  hmin            = grid.CoarsestH();
        const level_t lmin            = grid.MinLevel();
        const level_t lmax            = grid.MaxLevel();
        const long    global_num_quad = grid.global_num_quadrants();

        {
            const string id_name = "a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT)+ "_nu" + to_string(nu_);
            grid.DumpLevels(il, "data", id_name);
        }

        {  // WENO 3
            if (fix_weno_) {
                Advection<M_CONS, 3> adv(&vel, nu_);
                adv(&grid, &test, &dtest);
            } else {
                Advection<M_WENO_Z, 3> adv(&vel, nu_);
                adv(&grid, &test, &dtest);
            }
            m_log("error weno 3");
            // now, we need to check

            Error  error;
            real_t erri, err2;
            error.Norms(&grid, &dtest, &lambda_sol_rhs, &err2, &erri);

            const string scheme_name = (fix_weno_) ? "cons3" : "weno3";
            const string id_name     = scheme_name + "_a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + "_nu" + to_string(nu_);
            if (rank == 0) {
                string fname     = "data/conv_" + id_name + ".data";
                FILE*  file_diag = fopen(fname.c_str(), "a+");
                fprintf(file_diag, "%e %e %e %e %e %e %d %d %ld\n", grid.rtol(), grid.ctol(), hmin, hmax, err2, erri, lmin, lmax, global_num_quad);
                fclose(file_diag);
            }
            if (fix_weno_) {
                m_log("CONS-3: %e %e %e %e", hmin, hmax, err2, erri);
            } else {
                m_log("WENO-3: %e %e %e %e", hmin, hmax, err2, erri);
            }
        }
        {  // WENO 5
            if (fix_weno_) {
                Advection<M_CONS, 5> adv(&vel, nu_);
                adv(&grid, &test, &dtest);
            } else {
                Advection<M_WENO_Z, 5> adv(&vel, nu_);
                adv(&grid, &test, &dtest);
            }
            // now, we need to check
            Error  error;
            real_t erri, err2;
            error.Norms(&grid, &dtest, &lambda_sol_rhs, &err2, &erri);

            const string scheme_name = (fix_weno_) ? "cons5" : "weno5";
            const string id_name     = scheme_name + "_a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + "_nu" + to_string(nu_);

            if (rank == 0) {
                string fname     = "data/conv_" + id_name + ".data";
                FILE*  file_diag = fopen(fname.c_str(), "a+");
                fprintf(file_diag, "%e %e %e %e %e %e %d %d %ld\n", grid.rtol(), grid.ctol(), hmin, hmax, err2, erri, lmin, lmax, global_num_quad);
                fclose(file_diag);
            }
            if (fix_weno_) {
                m_log("CONS-5: %e %e %e %e", hmin, hmax, err2, erri);
            } else {
                m_log("WENO-5: %e %e %e %e", hmin, hmax, err2, erri);
            }
        }

        grid.DeleteField(&vel);
        grid.DeleteField(&test);
        grid.DeleteField(&dtest);
    }
    //-------------------------------------------------------------------------
    m_end;
}
