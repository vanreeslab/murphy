#include "clients/convergence_weno.hpp"

#include <list>
#include <string>

#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "operator/advection.hpp"
#include "operator/error.hpp"

using std::list;
using std::string;
using std::to_string;

#define N_PT 3

ConvergenceWeno::~ConvergenceWeno() {
    //-------------------------------------------------------------------------
    //-------------------------------------------------------------------------
}

void ConvergenceWeno::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    ilevel_   = m_max(param->init_lvl, 0);
    adapt_    = !param->no_adapt;
    fix_weno_ = param->fix_weno;

    // get the level
    level_start_ = param->init_lvl;
    level_min_   = param->level_min;
    level_max_   = param->level_max;

    eps_start_ = param->eps_start;
    delta_eps_ = param->delta_eps;
    //-------------------------------------------------------------------------
}

static const real_t sigma     = 0.05;
static const real_t radius    = 0.25;
static const real_t beta      = 3.0;
static const auto   freq      = std::vector<short_t>{};  //std::vector<short_t>{5, 1000};
static const auto   amp       = std::vector<real_t>{};   //std::vector<real_t>{0.2, 0.2};
static const real_t center[3] = {0.5, 0.5, 0.5};

// lambdas
static lambda_setvalue_t lambda_initcond = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
    // get the position
    real_t pos[3];
    block->pos(i0, i1, i2, pos);
    real_t* data = block->data(fid).Write(i0, i1, i2);
    // set value
    // data[0] = scalar_exp(pos, center, sigma);
    data[0] = scalar_compact_ring(pos, center, 2, radius, sigma, beta, freq, amp);
};
static lambda_error_t lambda_error = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
    // get the position
    real_t pos[3];
    block->pos(i0, i1, i2, pos);
    // set value
    return scalar_compact_ring(pos, center, 2, radius, sigma, beta, freq, amp);
};

void ConvergenceWeno::Run() {
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
    //......................................................................
    real_t  depsilon   = (adapt_) ? delta_eps_ : 0.0;
    real_t  epsilon    = eps_start_;
    short_t id_counter = 0;
    m_log("starting with epsilon = %e", epsilon);
    while (epsilon >= std::pow(2.0, -34)) {
        m_log("================================================================================");
        m_log("epsilon = %e", epsilon);
        //......................................................................
        // create a grid
        bool  period[3]   = {true, true, true};
        lid_t grid_len[3] = {1, 1, 1};
        Grid  grid(level_start_, period, grid_len, MPI_COMM_WORLD, nullptr);
        grid.level_limit(level_min_, level_max_);

        //......................................................................
        // create the scalar field and adapt if needed
        Field test("scalar", 1);
        {
            grid.AddField(&test);
            test.bctype(M_BC_ZERO);

            // set the values
            SetValue init(lambda_initcond);
            init(&grid, &test);

            // adapt
            if (adapt_) {
                grid.SetTol(epsilon, epsilon * 1e-2);
                grid.SetRecursiveAdapt(1);
                grid.Adapt(&test, &init);
            }
        }
        // grid is now adapted!
        //......................................................................
        // add the other fields
        Field vel("vel", 3);
        {
            grid.AddField(&vel);
            lambda_setvalue_t lambda_vel = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
                block->data(fid, 0).Write(i0, i1, i2)[0] = rand_vel[0];
                block->data(fid, 1).Write(i0, i1, i2)[0] = rand_vel[1];
                block->data(fid, 2).Write(i0, i1, i2)[0] = rand_vel[2];
            };
            const bidx_t ghost_len[2] = {3, 3};
            SetValue     vel_init(lambda_vel, ghost_len);
            vel_init(&grid, &vel);
        }

        //......................................................................
        // compute the derivative
        Field dtest("deriv", 1);
        grid.AddField(&dtest);

        // get the analytical solution
        const lda_t    normal     = 2;
        lambda_error_t lambda_sol = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
            // get the position
            real_t pos[3];
            block->pos(i0, i1, i2, pos);

            const lda_t idx = (normal + 1) % 3;
            const lda_t idy = (normal + 2) % 3;
            const lda_t idz = normal;

            const real_t gamma  = beta * sigma;
            const real_t x      = pos[idx] - center[idx];
            const real_t y      = pos[idy] - center[idy];
            const real_t z      = pos[idz] - center[idz];
            const real_t sigma2 = sigma * sigma;
            const real_t gamma2 = gamma * gamma;

            // compute the gaussian
            const real_t r_par = sqrt(pow(x, 2) + pow(y, 2)) - radius;
            const real_t r     = sqrt(pow(r_par, 2) + pow(z, 2));
            const real_t rho1  = r / sigma;
            const real_t rho2  = r / gamma;

            if (1e-15 < r && r < gamma) {
                const real_t exp_val  = exp(-pow(rho1, 2) / (1.0 - pow(rho2, 2)));
                const real_t exp_fact = (-2.0 * r / sigma2) / (1.0 - pow(rho2, 2)) + (-2.0 * pow(r, 3) / sigma2) / (gamma2 * pow(1.0 - pow(r, 2) / gamma2, 2));
                const real_t drdx     = (x * r_par) / ((r_par + radius) * r);
                const real_t drdy     = (y * r_par) / ((r_par + radius) * r);
                const real_t drdz     = (z / r);

                return -exp_val * exp_fact * (rand_vel[idx] * drdx + rand_vel[idy] * drdy + rand_vel[idz] * drdz);
            } else {
                return 0.0;
            }
        };

        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        const real_t  hmax            = grid.FinestH();
        const real_t  hmin            = grid.CoarsestH();
        const level_t lmin            = grid.MinLevel();
        const level_t lmax            = grid.MaxLevel();
        const long    global_num_quad = grid.global_num_quadrants();

        {
            const string id_name = "a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
            grid.DumpLevels(id_counter, "data", id_name);
        }

        {  // WENO 3
            if (fix_weno_) {
                Advection<M_CONS, 3> adv(&vel);
                adv(&grid, &test, &dtest);
            } else {
                Advection<M_WENO_Z, 3> adv(&vel);
                adv(&grid, &test, &dtest);
            }
            m_log("error weno 3");
            // now, we need to check

            Error  error;
            real_t erri, err2;
            error.Norms(&grid, &dtest, &lambda_sol, &err2, &erri);

            const string scheme_name = (fix_weno_) ? "cons3" : "weno3";
            const string id_name     = scheme_name + "_a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
            if (rank == 0) {
                string fname     = "data/conv_" + id_name + ".data";
                FILE*  file_diag = fopen(fname.c_str(), "a+");
                fprintf(file_diag, "%e %e %e %e %e %e %d %d %ld\n", grid.rtol(), grid.ctol(), hmin, hmax, err2, erri, lmin, lmax,global_num_quad);
                fclose(file_diag);
            }
            m_log("WENO-3: %e %e %e %e", hmin, hmax, err2, erri);
        }
        {  // WENO 5
            if (fix_weno_) {
                Advection<M_CONS, 5> adv(&vel);
                adv(&grid, &test, &dtest);
            } else {
                Advection<M_WENO_Z, 5> adv(&vel);
                adv(&grid, &test, &dtest);
            }
            // now, we need to check
            Error  error;
            real_t erri, err2;
            error.Norms(&grid, &dtest, &lambda_sol, &err2, &erri);

            const string scheme_name = (fix_weno_) ? "cons3" : "weno3";
            const string id_name     = scheme_name + "_a" + to_string(adapt_) + "_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);

            if (rank == 0) {
                string fname     = "data/conv_" + id_name + ".data";
                FILE*  file_diag = fopen(fname.c_str(), "a+");
                fprintf(file_diag, "%e %e %e %e %e %e %d %d %ld\n", grid.rtol(), grid.ctol(), hmin, hmax, err2, erri, lmin, lmax,global_num_quad);
                fclose(file_diag);
            }
            m_log("WENO-5: %e %e %e %e", hmin, hmax, err2, erri);
        }

        grid.DeleteField(&vel);
        grid.DeleteField(&test);
        grid.DeleteField(&dtest);

        // get the new epsilon
        epsilon *= depsilon;
        id_counter += 1;
    }
    //-------------------------------------------------------------------------
    m_end;
}
