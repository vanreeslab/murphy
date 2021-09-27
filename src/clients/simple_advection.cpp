#include "clients/simple_advection.hpp"

#include "operator/advection.hpp"
#include "operator/blas.hpp"
#include "operator/error.hpp"
#include "operator/xblas.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

static const lda_t  ring_normal = 2;
static const real_t sigma       = 0.10;
static const real_t radius      = 1.0;
static const real_t beta        = 3;
static const auto   freq        = std::vector<short_t>{};  //std::vector<short_t>{5, 25};
static const auto   amp         = std::vector<real_t>{};   //std::vector<real_t>{0.0, 0.1};
static const real_t center[3]   = {1.5, 1.5, 0.5};
static const real_t velocity[3] = {0.0, 0.0, 1.0};

static const lambda_setvalue_t lambda_velocity = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
    m_assert(fid->lda() == 3, "the velocity field must be a vector");
    block->data(fid, 0).Write(i0, i1, i2)[0] = velocity[0];
    block->data(fid, 1).Write(i0, i1, i2)[0] = velocity[1];
    block->data(fid, 2).Write(i0, i1, i2)[0] = velocity[2];
};

SimpleAdvection::~SimpleAdvection() {
    //-------------------------------------------------------------------------
    // delete the field
    m_profStart(prof_, "cleanup");
    grid_->DeleteField(vel_);
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

    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("SimpleAdvection") + to_string(comm_size) + string("ranks") + string("_w") + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
        prof_       = new Prof(name);
    }
    m_profStart(prof_, "init");

    // setup the grid
    bidx_t length[3] = {3, 3, 2};
    grid_            = new Grid(param->init_lvl, param->period, length, MPI_COMM_WORLD, prof_);

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
        block->data(fid).Write(i0, i1, i2)[0] = scalar_compact_ring(pos, center, ring_normal, radius, sigma, beta, freq, amp);
    };
    const bidx_t ghost_len_interp[2] = {m_max(grid_->interp()->nghost_front(), 3),
                                        m_max(grid_->interp()->nghost_back(), 3)};
    SetValue     ring(lambda_ring, ghost_len_interp);
    ring(grid_, scal_);

        // adapt the grid
    if (!no_adapt_) {
        // if the ctol is smaller than epsilon, just put epsilon
        grid_->SetTol(param->refine_tol, param->coarsen_tol);
        grid_->SetRecursiveAdapt(true);
        grid_->Adapt(scal_, &ring);
    }

    //-------------------------------------------------------------------------
    // set the velocity field
    vel_ = new Field("velocity", 3);
    vel_->bctype(M_BC_EXTRAP);
    vel_->is_temp(true);
    grid_->AddField(vel_);
    SetValue set_velocity(lambda_velocity, ghost_len_interp);
    set_velocity(grid_, vel_);

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
        advection = new Advection<M_WENO_Z, 3>(vel_, prof_);
        m_log("WENO-Z order 3 (cfl = %f)", advection->cfl_rk3());
    } else if (weno_ == 5 && !fix_weno_) {
        advection = new Advection<M_WENO_Z, 5>(vel_, prof_);
        m_log("WENO-Z order 5 (cfl = %f)", advection->cfl_rk3());
    } else if (weno_ == 3 && fix_weno_) {
        advection = new Advection<M_CONS, 3>(vel_, prof_);
        m_log("CONS order 3 (cfl = %f)", advection->cfl_rk3());
    } else if (weno_ == 5 && fix_weno_) {
        advection = new Advection<M_CONS, 5>(vel_, prof_);
        m_log("CONS order 5 (cfl = %f)", advection->cfl_rk3());
    } else {
        advection = nullptr;
        m_assert(false, "weno order = %d not valid", weno_);
    }

    // time integration
    iter_t         iter = 0;
    real_t        t    = tstart_;
    const RK3_TVD rk3(grid_, scal_, advection, prof_, cfl_);

    // let's gooo
    m_profStart(prof_, "run");
    const real_t wtime_start = MPI_Wtime();
    while (t < tfinal_ && iter < iter_max()) {
        m_log("--------------------------------------------------------------------------------");
        //................................................
        // adapt the mesh
        if ((iter % iter_adapt() == 0) && (!no_adapt_)) {
            m_log("---- adapt mesh");
            m_log_level_plus;
            m_profStart(prof_, "adapt");
            if (!grid_on_sol_) {
                grid_->SetRecursiveAdapt(true);
                grid_->Adapt(scal_);

            } else {
                m_assert(false, "this option is not supported without a solution field");
            }
            m_profStop(prof_, "adapt");

            // reset the velocity
            m_profStart(prof_, "set velocity");
            // set the velocity field
            const bidx_t ghost_len_interp[2] = {m_max(grid_->interp()->nghost_front(), 3),
                                                m_max(grid_->interp()->nghost_back(), 3)};
            SetValue     set_velocity(lambda_velocity, ghost_len_interp);
            set_velocity(grid_, vel_);
            m_assert(vel_->ghost_status(ghost_len_interp), "the velocity ghosts must have been computed");
            m_profStop(prof_, "set velocity");
            m_log_level_minus;
        }
        // we run the first diagnostic if not done yet, it's usefull to get a sense of what is going on with the adaptation
        if (iter == 0) {
            m_profStart(prof_, "diagnostics");
            m_log("---- run diag");
            m_log_level_plus;
            real_t wtime_now = MPI_Wtime();
            Diagnostics(t, 0.0, iter, wtime_now - wtime_start);
            m_profStop(prof_, "diagnostics");
            m_log_level_minus;
        }

        //................................................
        m_log("---- do time-step");
        m_log_level_plus;
        //................................................
        // get the time-step given the field
        m_profStart(prof_, "compute dt");
        real_t dt = rk3.ComputeDt(advection, vel_);
        m_profStop(prof_, "compute dt");

        // dump some info
        real_t wtime_now = MPI_Wtime();
        m_log("RK3 - time = %f/%f - step %d/%d - dt = %e - wtime = %e", t, tfinal_, iter, iter_max(), dt, wtime_now - wtime_start);

        //................................................
        // advance in time
        m_profStart(prof_, "do dt");
        rk3.DoDt(dt, &t);
        iter++;
        m_profStop(prof_, "do dt");

        m_log_level_minus;

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
    }
    m_profStop(prof_, "run");
    // run the last diag
    if (iter % iter_diag() != 0) {
        m_profStart(prof_, "diagnostics");
        real_t time_now = MPI_Wtime();
        Diagnostics(t, 0.0, iter, time_now);
        m_profStop(prof_, "diagnostics");
    }

    // free the advection
    delete advection;

    //-------------------------------------------------------------------------
    m_end;
}

void SimpleAdvection::Diagnostics(const real_t time, const real_t dt, const iter_t iter, const real_t wtime) {
    m_begin;
    m_assert(scal_->lda() == 1, "the scalar field must be scalar");
    //-------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // if the folder does not exist, create it
    struct stat st = {0};
    if (rank == 0 && stat(folder_diag_.c_str(), &st) == -1) {
        mkdir(folder_diag_.c_str(), 0770);
    }

    //................................................
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

    const real_t new_center[3] = {center[0] + time * velocity[0],
                                  center[1] + time * velocity[1],
                                  center[2] + time * velocity[2]};

    lambda_error_t lambda_ring = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block) -> real_t {
        // get the position
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        return scalar_compact_ring(pos, new_center, ring_normal, radius, sigma, beta, freq, amp);
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
    int  log_ratio = log10(grid_->ctol())-log10(grid_->rtol());
    string tag        = "w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + "_a" + to_string(adapt_freq) + weno_name + to_string(weno_) + "_tol"+to_string(-log_ratio);

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
    grid_->DistributionDetails(iter, folder_diag_, tag, scal_, 128, 1.0);
    m_profStop(prof_, "dump det histogram");

    m_profStart(prof_, "dump field");
    if (iter % iter_dump() == 0 && iter != 0) {
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
