#include "clients/simple_advection.hpp"

#include "operator/advection.hpp"
#include "operator/blas.hpp"
#include "operator/diagnostics.hpp"
#include "operator/error.hpp"
#include "operator/xblas.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

#define M_DIFF_ONLY 0

// ring
static const lda_t  ring_normal    = 2;
static const real_t ring_sigma     = 0.05;
static const real_t ring_radius    = 0.75;
static const real_t ring_beta      = 2;
static const real_t ring_center[3] = {1.5, 1.5, 0.5};
// exponential
static const real_t exp_sigma     = 0.2;
static const real_t exp_beta      = 6;
static const real_t exp_center[3] = {1.5, 1.5, 1.5};
// pure advection in Z
#if M_DIFF_ONLY
static const real_t velocity[3] = {0.0, 0.0, 0.0};
#else
static const real_t velocity[3] = {0.0, 0.0, 1.0};
#endif

static const lambda_expr_t lambda_velocity = [](const real_t x, const real_t y, const real_t z, const lda_t ida) -> real_t {
#if M_DIFF_ONLY
    return 0.0;
#else
    return (ida == 2);
#endif
};

SimpleAdvection::~SimpleAdvection() {
    //-------------------------------------------------------------------------
    // delete the field
    m_profStart(prof_, "cleanup");
    // grid_->DeleteField(vel_);
    grid_->DeleteField(scal_);

    delete vel_;
    delete scal_;
    delete grid_;

    m_profStop(prof_, "cleanup");

    if (!(prof_ == nullptr)) {
        prof_->Disp();
        delete prof_;
    }
    //-------------------------------------------------------------------------
}

void SimpleAdvection::InitParam(ParserArguments* param) {
    //--------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // take the no adaptation
    no_adapt_    = param->no_adapt;
    grid_on_sol_ = param->grid_on_sol;
    weno_        = param->weno;
    fix_weno_    = param->fix_weno;
    m_assert(weno_ == 3 || weno_ == 5, "the weno order must be 3 or 5");

    // cfl
    cfl_ = param->cfl_max;
    // diffusion: Re = U * L / nu if neg it means no diffusion
    nu_ = (param->reynolds < 0.0) ? 0.0 : ((1.0 * exp_sigma) / param->reynolds);
    m_log("nu = %e",nu_);

    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("SimpleAdvection") + to_string(comm_size) + string("ranks") + string("_w") + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + "_re" + to_string(param->reynolds);;
        prof_       = new Prof(name);
    }
    m_profStart(prof_, "init");

    // setup the grid
    bidx_t length[3] = {3,3,6};
    bool   period[3] = {false, false, false};
    grid_          = new Grid(param->init_lvl, period, length,M_GRIDBLOCK, MPI_COMM_WORLD, prof_);

    // set the min/max level
    grid_->level_limit(param->level_min, param->level_max);

    //-------------------------------------------------------------------------
    scal_ = new Field("scalar", 1);
    scal_->bctype(M_BC_ZERO);
    // scal_->bctype(M_BC_EXTRAP);
    grid_->AddField(scal_);

    // setup the scalar ring
    lambda_setvalue_t lambda_ring = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);

        real_t value = 0.0;
        //  value += scalar_compact_ring(pos, ring_center, ring_normal, ring_radius, ring_sigma, ring_beta);
        value += scalar_diff_exp(pos, exp_center, exp_sigma,nu_,0.0);
        block->data(fid, 0)(i0, i1, i2) = value;
    };

    const bidx_t ghost_len_interp[2] = {m_max(grid_->interp()->nghost_front(), 3),
                                        m_max(grid_->interp()->nghost_back(), 3)};
    SetValue     ring(lambda_ring, ghost_len_interp);
    ring(grid_, scal_);

    m_log("set values with center = %f %f %f, sigma = %f, beta = %f", exp_center[0], exp_center[1], exp_center[2], exp_sigma, exp_beta);

    // reinterpret the coarsen tol
    optimal_tol_       = param->optimal_tol;
    real_t coarsen_tol = (param->coarsen_tol < 0.0) ? 0.0 : ((optimal_tol_) ? (param->refine_tol / pow(2.0, M_WAVELET_N)) : (param->coarsen_tol));
    refine_only_       = (param->coarsen_tol < 0.0);

    grid_->SetTol(param->refine_tol, coarsen_tol);
    // adapt the grid
    // the refine only starts from the lossless compression
    if (!no_adapt_) {
        // if the ctol is smaller than epsilon, just put epsilon
        grid_->SetRecursiveAdapt(true);
        m_assert(!(refine_only_ && grid_->ctol() > std::numeric_limits<real_t>::epsilon()), "If I am refine only, my coarsening tolerance MUST BE lossless (i.e. 0.0)");
        grid_->Adapt(scal_, &ring);
    }

    //-------------------------------------------------------------------------
    // set the velocity field
    vel_ = new Field("velocity", 3);
    vel_->bctype(M_BC_EXTRAP);
    vel_->is_expr(true);
    grid_->SetExpr(vel_, lambda_velocity);

    tstart_ = param->time_start;
    tfinal_ = param->time_final;
    m_profStop(prof_, "init");
    //-------------------------------------------------------------------------
}

void SimpleAdvection::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // advection
    RKFunctor* advection;
    if (weno_ == 3 && !fix_weno_) {
        advection = new Advection<M_WENO_Z, 3>(vel_, nu_, prof_);
        m_log("WENO-Z order 3 (cfl = %f)", advection->cfl_rk3());
    } else if (weno_ == 5 && !fix_weno_) {
        advection = new Advection<M_WENO_Z, 5>(vel_, nu_, prof_);
        m_log("WENO-Z order 5 (cfl = %f)", advection->cfl_rk3());
    } else if (weno_ == 3 && fix_weno_) {
        advection = new Advection<M_CONS, 3>(vel_, nu_, prof_);
        m_log("CONS order 3 (cfl = %f)", advection->cfl_rk3());
    } else if (weno_ == 5 && fix_weno_) {
        advection = new Advection<M_CONS, 5>(vel_, nu_, prof_);
        m_log("CONS order 5 (cfl = %f)", advection->cfl_rk3());
    } else {
        advection = nullptr;
        m_assert(false, "weno order = %d not valid", weno_);
    }

    // time integration
    iter_t        iter = 0;
    real_t        t    = tstart_;
    const RK3_TVD rk3(grid_, scal_, advection, prof_, cfl_);

    t_deterr_ = 0.25;
    t_deterr_accum_ = 0.0;

    // let's gooo
    m_profStart(prof_, "run");
    const real_t wtime_start = MPI_Wtime();
    real_t       dt          = 0.0;
    while (t < tfinal_ && iter < iter_max()) {
        m_log("--------------------------------------------------------------------------------");
        //................................................
        // adapt the mesh
        if (iter % iter_adapt() == 0) {
            if (!no_adapt_) {
                m_log("---- adapt mesh");
                m_log_level_plus;
                m_profStart(prof_, "adapt");
                if (!grid_on_sol_) {
                    grid_->SetRecursiveAdapt(true);
                    // grid_->Adapt(scal_);
                    if (refine_only_) {
                        grid_->Refine(scal_);
                    } else {
                        grid_->Adapt(scal_);
                    }
                } else {
                    m_assert(false, "this option is not supported without a solution field");
                }
                m_profStop(prof_, "adapt");

                // reset the velocity
                m_profStart(prof_, "set velocity");
                // set the velocity field
                const bidx_t ghost_len_interp[2] = {m_max(grid_->interp()->nghost_front(), 3),
                                                    m_max(grid_->interp()->nghost_back(), 3)};
                // SetValue     set_velocity(lambda_velocity, ghost_len_interp);
                // set_velocity(grid_, vel_);
                // grid_->SetExpr(vel_,lambda_velocity);
                m_assert(vel_->ghost_status(ghost_len_interp), "the velocity ghosts must have been computed");
                
                m_profStop(prof_, "set velocity");
                m_log_level_minus;
            }
        }
        //................................................
        // diagnostics, dumps, whatever
        if (iter % iter_diag() == 0) {
            m_log("---- run diag");
            m_log_level_plus;
            m_profStart(prof_, "diagnostics");
            real_t time_now = MPI_Wtime();
            Diagnostics(t, dt, iter, time_now - wtime_start);
            m_profStop(prof_, "diagnostics");
            m_log_level_minus;
        }

        //................................................
        m_log("---- do time-step");
        m_log_level_plus;
        //................................................
        // get the time-step given the field
        m_profStart(prof_, "compute dt");
#if M_DIFF_ONLY
        dt = rk3.ComputeDt(advection, 0.0, nu_);
#else
        dt = rk3.ComputeDt(advection, 1.0, nu_);
#endif
        m_profStop(prof_, "compute dt");

        // dump some info
        real_t wtime_now = MPI_Wtime();
        //................................................
        // advance in time
        m_profStart(prof_, "do dt");
        rk3.DoDt(dt, &t);
        iter++;
        m_log("now -> time = %f/%f - step %d/%d - dt = %e - wtime = %e", t, tfinal_, iter, iter_max(), dt, wtime_now - wtime_start);
        m_profStop(prof_, "do dt");
        m_log_level_minus;
    }
    m_profStop(prof_, "run");
    // run the last diag
    if (iter % iter_diag() != 0) {
        m_profStart(prof_, "diagnostics");
        real_t time_now = MPI_Wtime();
        Diagnostics(t, dt, iter, time_now);
        m_profStop(prof_, "diagnostics");
    }

    // free the advection
    delete advection;

    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief 
 * 
 * @param time current time
 * @param dt the value of DT we have just done
 * @param iter 
 * @param wtime 
 */
void SimpleAdvection::Diagnostics(const real_t time, const real_t dt, const iter_t iter, const real_t wtime) {
    m_begin;
    m_assert(scal_->lda() == 1, "the scalar field must be scalar");
    //-------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //..........................................................................
    // update the time accum
    t_deterr_accum_ += dt;
    m_log("time is now, %e, accum = %e",time,t_deterr_accum_);

    //..........................................................................
    // if the folder does not exist, create it
    struct stat st = {0};
    if (rank == 0 && stat(folder_diag_.c_str(), &st) == -1) {
        mkdir(folder_diag_.c_str(), 0770);
    }

    //..........................................................................
    m_profStart(prof_, "cmpt moments");
    real_t  moment0;
    real_t  moment1[3];
    BMoment moments;
    grid_->GhostPull(scal_, &moments);
    moments(grid_, scal_, &moment0, moment1);
    m_profStop(prof_, "cmpt moments");
    //................................................
    // m_profStart(prof_, "cmpt mean");
    // real_t mean_val;
    // BAvg   mean;
    // mean(grid_, scal_, &mean_val);
    // m_profStop(prof_, "cmpt mean");

    //................................................
    // get the solution at the given time
    const real_t   center_shift[3] = {time * velocity[0], time * velocity[1], time * velocity[2]};
    lambda_error_t lambda_ring     = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);

        real_t new_ring_center[3] = {ring_center[0] + center_shift[0], ring_center[1] + center_shift[1], ring_center[2] + center_shift[2]};
        real_t new_exp_center[3]  = {exp_center[0] + center_shift[0], exp_center[1] + center_shift[1], exp_center[2] + center_shift[2]};

        real_t value = 0.0;
        // value += scalar_compact_ring(pos, new_ring_center, ring_normal, ring_radius, ring_sigma, ring_beta);
        // value += scalar_compact_exp(pos, new_exp_center, exp_sigma, exp_beta);
        value += scalar_diff_exp(pos, new_exp_center, exp_sigma, nu_, time);
        return value;
    };
    // compute the error
    real_t err2, erri;
    Error  error;
    m_profStart(prof_, "cmpt error");
    error.Norms(grid_, scal_, &lambda_ring, &err2, &erri);
    m_profStop(prof_, "cmpt error");

    real_t   density = 0.0;
    BDensity dense;
    dense(grid_, &density);

    // get the min and max detail coefficients
    real_t maxmin_details[2];
    grid_->MaxMinDetails(scal_, maxmin_details);

    // tag
    lid_t  adapt_freq = no_adapt_ ? 0 : iter_adapt();
    string weno_name  = fix_weno_ ? "_cons" : "_weno";
    // int    log_ratio  = (grid_->ctol() > std::numeric_limits<real_t>::epsilon()) ? (log(grid_->ctol()) - log(grid_->rtol())) : -0;
    int    log_ratio  = (grid_->ctol() > std::numeric_limits<real_t>::epsilon()) ? (grid_->rtol()/grid_->ctol()) : 0;
    string adapt_tag  = (refine_only_) ? ("_r") : ("_a");
    string tol_tag    = (optimal_tol_) ? ("_opt") : ("_tol");
    // string tag        = "w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + adapt_tag + to_string(adapt_freq) + weno_name + to_string(weno_) + tol_tag + to_string(-log_ratio);
    string tag        = "w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + adapt_tag + to_string(adapt_freq) + weno_name + to_string(weno_) + tol_tag + to_string(log_ratio);

    // open the file
    m_profStart(prof_, "dump diag");
    FILE*   file_error;
    FILE*   file_diag;
    level_t min_level       = grid_->MinLevel();
    level_t max_level       = grid_->MaxLevel();
    long    global_num_quad = grid_->global_num_quadrants();
    m_log("iter = %6.6d time = %e: levels = (%d , %d -> %e), errors = (%e , %e), wtime = %e", iter, time, min_level, max_level, density, err2, erri, wtime);
    if (rank == 0) {
        // lid_t adapt_freq = no_adapt_ ? 0 : iter_adapt();
        // string weno_name = fix_weno_ ? "_cons" : "_weno";
        string file_name = "diag_" + tag + ".data";
        file_diag        = fopen(string(folder_diag_ + "/" + file_name).c_str(), "a+");
        fprintf(file_diag, "%6.6d;%e;%e;%ld;%d;%d", iter, time, dt, global_num_quad, min_level, max_level);
        fprintf(file_diag, ";%e", wtime);
        fprintf(file_diag, ";%e;%e", grid_->rtol(), grid_->ctol());
        fprintf(file_diag, ";%e;%e", err2, erri);
        fprintf(file_diag, ";%e;%e", maxmin_details[0], maxmin_details[1]);
        fprintf(file_diag, ";%e", density);
        fprintf(file_diag, ";%e;%e;%e;%e", moment0, moment1[0], moment1[1], moment1[2]);
        fprintf(file_diag, "\n");
        fclose(file_diag);
    }
    m_profStop(prof_, "dump diag");

    m_profStart(prof_, "dump levels");
    grid_->DumpLevels(iter, folder_diag_, string("_" + tag));
    m_profStop(prof_, "dump levels");

    m_profStart(prof_, "dump det histogram");
    if ((t_deterr_accum_ - t_deterr_) > ((-1000.0) * std::numeric_limits<real_t>::epsilon()) || iter == 0 || time >= tfinal_) {
        DetailVsError distr(grid_->interp());
        distr(iter, folder_diag_, tag, grid_, scal_, &lambda_ring);
        // reset the correct increment
        const iter_t n_dump = ((time + m_min(dt, t_deterr_) / 10.0) / t_deterr_);
        t_deterr_accum_     = time - n_dump * t_deterr_;
        m_log("dump deterr at time %e - iter %d -> accum is now %e", time, iter, t_deterr_accum_);
    }
    m_profStop(prof_, "dump det histogram");

    m_profStart(prof_, "dump field");
    if ((iter % iter_dump() == 0) && iter > 0) {
        // dump the vorticity field
        IOH5 dump(folder_diag_);
        grid_->GhostPull(scal_, ghost_len_ioh5);
        dump(grid_, scal_, iter);
        // dump(grid_, sol_, iter);

        // dump the details
        if (dump_detail()) {
            Field details("detail", scal_->lda());
            details.bctype(M_BC_EXTRAP);
            grid_->AddField(&details);
            grid_->StoreDetails(scal_, &details);

            grid_->GhostPull(&details, ghost_len_ioh5);
            // IOH5 dump("data");
            dump(grid_, &details, iter);

            grid_->DeleteField(&details);
        }
    }
    m_profStop(prof_, "dump field");
    //-------------------------------------------------------------------------
    m_end;
}

// void SimpleAdvection::GridDetErr(const real_t time, const real_t dt, const iter_t iter, const real_t wtime) {
//     m_begin;
//     m_assert(scal_->lda() == 1, "the scalar field must be scalar");
//     //-------------------------------------------------------------------------
//     rank_t rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//     //................................................
//     // get the solution at the given time

//     const real_t new_center[3] = {center[0] + time * velocity[0],
//                                   center[1] + time * velocity[1],
//                                   center[2] + time * velocity[2]};

//     lambda_error_t lambda_ring = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
//         // get the position
//         real_t pos[3];
//         block->pos(i0, i1, i2, pos);
// #ifdef CASE_RING
//         return scalar_compact_ring(pos, new_center, ring_normal, radius, sigma, beta, freq, amp);
// #endif
// #ifdef CASE_TUBE
//         return scalar_compact_tube(pos, new_center, sigma, beta, ring_normal);
// #endif
//     };
    
//     // tag
//     lid_t  adapt_freq = no_adapt_ ? 0 : iter_adapt();
//     string weno_name  = fix_weno_ ? "_cons" : "_weno";
//     int    log_ratio  = (grid_->ctol() > std::numeric_limits<real_t>::epsilon()) ? (log10(grid_->ctol()) - log10(grid_->rtol())) : -0;
//     string adapt_tag  = (refine_only_) ? ("_r") : ("_a");
//     string tag        = "w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + adapt_tag + to_string(adapt_freq) + weno_name + to_string(weno_) + "_tol" + to_string(-log_ratio);

//     if ((iter % iter_adapt() == 0) ||(iter % iter_adapt() == 1)) {
//         m_profStart(prof_, "dump det histogram");
//         DetailVsError distr(grid_->interp());
//         distr(iter, folder_diag_, tag, grid_, scal_, &lambda_ring);
//         m_profStop(prof_, "dump det histogram");
//     }
//     //-------------------------------------------------------------------------
//     m_end;
// }
