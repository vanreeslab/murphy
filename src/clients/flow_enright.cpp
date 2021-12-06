#include "clients/flow_enright.hpp"

#include "operator/advection.hpp"
#include "operator/blas.hpp"
#include "operator/diagnostics.hpp"
#include "operator/error.hpp"
#include "operator/xblas.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"
#include "core/exprdata.hpp"

using std::string;
using std::to_string;

// exponential
static const real_t exp_sigma     = 0.075;
static const real_t exp_beta      = 2.0;
static const real_t exp_center[3] = {0.35, 0.35, 0.35};

static const lambda_expr_t lambda_velocity = [](const real_t x, const real_t y, const real_t z, const lda_t ida) -> real_t {
    //--------------------------------------------------------------------------
    const real_t xc = x;  //- exp_center[0];
    const real_t yc = y;  //- exp_center[1];
    const real_t zc = z;  //- exp_center[2];
    const real_t vx = (+2.0) * pow(sin(M_PI * xc), 2) * sin(2.0 * M_PI * yc) * sin(2.0 * M_PI * zc);
    const real_t vy = (-1.0) * sin(2.0 * M_PI * xc) * pow(sin(M_PI * yc), 2) * sin(2.0 * M_PI * zc);
    const real_t vz = (-1.0) * sin(2.0 * M_PI * xc) * sin(2.0 * M_PI * yc) * pow(sin(M_PI * zc), 2);

    return vx * (ida == 0) + vy * (ida == 1) + vz * (ida == 2);
    //--------------------------------------------------------------------------
};

FlowEnright::~FlowEnright() {
    m_begin;
    //--------------------------------------------------------------------------
    grid_->DeleteField(scal_);

    delete scal_;
    delete vel_;
    delete advection_;
    delete rk3_;
    //--------------------------------------------------------------------------
    m_end;
}

void FlowEnright::Setup(ParserArguments* param) {
    m_begin;
    //--------------------------------------------------------------------------
    m_profStart(prof_, "setup");

    //..........................................................................
    weno_     = param->weno;
    fix_weno_ = param->fix_weno;
    m_assert(weno_ == 3 || weno_ == 5, "the weno order must be 3 or 5");

    // cfl
    cfl_ = param->cfl_max;

    //..........................................................................
    scal_ = new Field("scalar", 1);
    scal_->bctype(M_BC_ZERO);
    grid_->AddField(scal_);

    // setup the scalar ring
    lambda_setvalue_t lambda_exp = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);

        real_t value                    = scalar_compact_exp(pos, exp_center, exp_sigma, exp_beta);
        block->data(fid, 0)(i0, i1, i2) = value;
    };

    const bidx_t ghost_len_interp[2] = {m_max(grid_->interp()->nghost_front(), 3),
                                        m_max(grid_->interp()->nghost_back(), 3)};
    SetValue     init_scal(lambda_exp, ghost_len_interp);
    init_scal(grid_, scal_);

    m_log("set values with center = %f %f %f, sigma = %f, beta = %f", exp_center[0], exp_center[1], exp_center[2], exp_sigma, exp_beta);

    // adapt the grid
    // the refine only starts from the lossless compression
    if (!no_adapt_) {
        // if the ctol is smaller than epsilon, just put epsilon
        grid_->SetRecursiveAdapt(true);
        grid_->Adapt(scal_, &init_scal);
    }

    //..........................................................................
    // set the velocity field
    vel_ = new Field("velocity", 3);
    vel_->bctype(M_BC_EXTRAP);
    vel_->is_expr(true);
    grid_->SetExpr(vel_, lambda_velocity);

    //..........................................................................
    // advection
    if (weno_ == 3 && !fix_weno_) {
        advection_ = new Advection<M_WENO_Z, 3>(vel_, 0.0, prof_);
        m_log("WENO-Z order 3 (cfl = %f)", advection_->cfl_rk3());
    } else if (weno_ == 5 && !fix_weno_) {
        advection_ = new Advection<M_WENO_Z, 5>(vel_, 0.0, prof_);
        m_log("WENO-Z order 5 (cfl = %f)", advection_->cfl_rk3());
    } else if (weno_ == 3 && fix_weno_) {
        advection_ = new Advection<M_CONS, 3>(vel_, 0.0, prof_);
        m_log("CONS order 3 (cfl = %f)", advection_->cfl_rk3());
    } else if (weno_ == 5 && fix_weno_) {
        advection_ = new Advection<M_CONS, 5>(vel_, 0.0, prof_);
        m_log("CONS order 5 (cfl = %f)", advection_->cfl_rk3());
    } else {
        advection_ = nullptr;
        m_assert(false, "weno order = %d not valid", weno_);
    }

    rk3_ = new RK3_TVD(grid_, scal_, advection_, prof_, cfl_);

    m_profStop(prof_, "setup");
    //--------------------------------------------------------------------------
    m_end;
}

void FlowEnright::DoTimeStep(real_t* time, real_t* dt) {
    m_begin;
    //--------------------------------------------------------------------------
    dt[0] = rk3_->ComputeDt(advection_, 1.0, 0.0);
    rk3_->DoDt(dt[0], time);
    //--------------------------------------------------------------------------
    m_end;
}

void FlowEnright::Adapt(const real_t time, const real_t dt) {
    m_begin;
    //--------------------------------------------------------------------------
    grid_->SetRecursiveAdapt(true);
    grid_->Adapt(scal_);

#ifndef NDEBUG
    // check the ghost status
    const bidx_t ghost_len_interp[2] = {m_max(grid_->interp()->nghost_front(), 3),
                                        m_max(grid_->interp()->nghost_back(), 3)};
    m_assert(vel_->ghost_status(ghost_len_interp), "the velocity ghosts must have been computed");
#endif

    //--------------------------------------------------------------------------
    m_end;
}

void FlowEnright::Diagnostics(const real_t time, const real_t dt, const lid_t iter, const real_t wtime) {
    m_begin;
    //--------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //..........................................................................
    // bump the diagnostic time, the time-step has been done
    time_accum_field_ += dt;

    //..........................................................................
    m_profStart(prof_, "cmpt moments");
    real_t  moment0;
    real_t  moment1[3];
    BMoment moments;
    grid_->GhostPull(scal_, &moments);
    moments(grid_, scal_, &moment0, moment1);
    m_profStop(prof_, "cmpt moments");

    //..........................................................................
    // if the folder does not exist, create it
    struct stat st = {0};
    if (rank == 0 && stat(folder_diag_.c_str(), &st) == -1) {
        mkdir(folder_diag_.c_str(), 0770);
    }

    m_profStart(prof_, "dump diag");
    FILE*   file_error;
    FILE*   file_diag;
    level_t min_level       = grid_->MinLevel();
    level_t max_level       = grid_->MaxLevel();
    long    global_num_quad = grid_->global_num_quadrants();
    m_log("iter = %6.6d time = %e: levels = (%d , %d), wtime = %e", iter, time, min_level, max_level, wtime);

    string tag = this->name() + "w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
    if (rank == 0) {
        // lid_t adapt_freq = no_adapt_ ? 0 : iter_adapt();
    // string weno_name = fix_weno_ ? "_cons" : "_weno";
        string file_name = "diag_" + tag + ".data";
        file_diag        = fopen(string(folder_diag_ + "/" + file_name).c_str(), "a+");
        fprintf(file_diag, "%6.6d;%e;%e;%ld;%d;%d", iter, time, dt, global_num_quad, min_level, max_level);
        fprintf(file_diag, ";%e", wtime);
        fprintf(file_diag, ";%e;%e", grid_->rtol(), grid_->ctol());
        fprintf(file_diag, ";%e;%e;%e;%e", moment0, moment1[0], moment1[1], moment1[2]);
        fprintf(file_diag, "\n");
        fclose(file_diag);
    }
    m_profStop(prof_, "dump diag");

    //..........................................................................
    m_profStart(prof_, "dump field");
    const bool dump_field = (iter == 0) ||
                            (time == tfinal_) ||
                            (m_fgeq(time_accum_field_, time_dump_field_));

    m_log("should I dump? %d: %e vs %e",m_fgeq(time_accum_field_, time_dump_field_),time_accum_field_,time_dump_field_);
    if (dump_field) {
        m_log("dump scalar");
        // dump the vorticity field
        IOH5 dump(folder_diag_);
        grid_->GhostPull(scal_, ghost_len_ioh5);
        dump(grid_, scal_, iter);

        // reset the correct increment
        const iter_t n_dump = ((time + m_min(dt, time_dump_field_) / 10.0) / time_dump_field_);
        time_accum_field_   = time - n_dump * time_dump_field_;
    }
    m_profStop(prof_, "dump field");

    //--------------------------------------------------------------------------
    m_end;
}
 