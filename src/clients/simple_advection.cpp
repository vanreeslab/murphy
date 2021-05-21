#include "clients/simple_advection.hpp"

#include "operator/advection.hpp"
#include "operator/blas.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"
#include "operator/xblas.hpp"
#include "operator/error.hpp"

using std::string;
using std::to_string;

SimpleAdvection::~SimpleAdvection() {
    //-------------------------------------------------------------------------
    if (!prof_.IsEmpty()) {
        prof_->Disp();
        prof_.Free();
    }

    ring_.Free();
    vel_field_.Free();

    vel_.Free();
    scal_.Free();
    sol_.Free();
    grid_.Free();

    m_log("Navier Stokes is dead");
    //-------------------------------------------------------------------------
}

void SimpleAdvection::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // take the no adaptation
    no_adapt_    = param->no_adapt;
    no_weno_     = param->no_weno;
    grid_on_sol_ = param->grid_on_sol;
    weno_5_      = param->weno_5;

    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("SimpleAdvection") + to_string(comm_size) + string("ranks") + string("_w") + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
        prof_.Alloc(name);
    }

    // setup the grid
    bool   period[3] = {false, false, false};
    const lid_t L[3]      = {1, 1, 1};
    grid_.Alloc(param->init_lvl, period, L, MPI_COMM_WORLD, prof_);

    m_assert(grid_.IsOwned(), "the grid must be owned");

    // set the min/max level
    grid_->level_limit(param->level_min, param->level_max);

    // get the fields
    vel_.Alloc("velocity", 3);
    vel_->bctype(M_BC_EXTRAP);
    vel_->is_temp(true);
    scal_.Alloc("scalar", 1);
    scal_->bctype(M_BC_EXTRAP);
    grid_->AddField(vel_);
    grid_->AddField(scal_);

    // add the solution as temp
    sol_.Alloc("sol", 1);
    sol_->is_temp(true);
    grid_->AddField(sol_);

    // setup the scalar ring
    real_t velocity[3] = {0.0, 0.0, 1.0};
    real_t center[3]   = {0.5, 0.5, 0.5};
    ring_.Alloc(param->vr_normal, center, param->vr_sigma, param->vr_radius, velocity, grid_->interp());
    ring_->Profile(prof_);
    ring_->SetTime(0.0);
    (*ring_)(grid_, scal_);

    // finish the grid
    if (!no_adapt_) {
        grid_->SetTol(param->refine_tol, param->coarsen_tol);
        grid_->SetRecursiveAdapt(true);
        grid_->Adapt(scal_, ring_);
    }

    // setup the velocity, 1.0 in every direction
    const lid_t  deg[3]   = {0, 0, 0};
    const real_t dir[3]   = {1.0, 0.0, 0.0};
    const real_t shift[3] = {0.0, 0.0, 0.0};
    vel_field_.Alloc(deg, dir, shift);
    (*vel_field_)(grid_, vel_, 2);
    // take the ghosts
    grid_->GhostPull(vel_);

    // IOH5 dump(folder_diag_);
    // dump(grid_(), vel_(),0);
    // dump(grid_(), scal_(),0);

    //-------------------------------------------------------------------------
}

void SimpleAdvection::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // time
    lid_t  iter    = 0;
    real_t t_start = 0.0;
    real_t t_final = 0.1;  //2.0/8.0;
    real_t t       = 0.0;

    RKFunctor* advection;
    if (no_weno_) {
        advection = new Advection<M_ADV_CONS_VEL, 3>(vel_);
        m_log("conservative advection chosen, cfl = %f", advection->cfl_rk3());
    } else if (weno_5_) {
        advection = new Advection<M_ADV_WENO_VEL, 5>(vel_);
        m_log("WENO 5 advection chosen, cfl = %f", advection->cfl_rk3());
    } else {
        advection = new Advection<M_ADV_WENO_VEL, 3>(vel_);
        m_log("WENO 3 advection chosen, cfl = %f", advection->cfl_rk3());
    }
    RK3_TVD rk3(grid_, scal_, advection, prof_);

    m_log("advection cfl is %e",advection->cfl_rk3());

    real_t time_start = MPI_Wtime();

    // let's gooo
    m_profStart(prof_(), "run");
    while (t < t_final && iter < iter_max()) {
        m_log("--------------------");
        //................................................
        // adapt the mesh
        if (iter % iter_adapt() == 0) {
            if (!no_adapt_) {
                m_log("---- adapt mesh");
                m_profStart(prof_(), "adapt");
                if (!grid_on_sol_) {
                    grid_->Adapt(scal_);
                } else {
                    // update the solution and refine
                    ring_->SetTime(t);
                    (*ring_)(grid_, sol_);

                    grid_->Adapt(sol_);
                }
                m_profStop(prof_(), "adapt");

                // reset the velocity
                (*vel_field_)(grid_, vel_, 2);
                grid_->GhostPull(vel_);
                m_assert(vel_->ghost_status(), "the velocity ghosts must have been computed");
            }
        }
        // we run the first diagnostic
        if (iter == 0) {
            m_profStart(prof_(), "diagnostics");
            m_log("---- run diag");
            real_t time_now = MPI_Wtime();
            Diagnostics(t, 0, iter, time_now - time_start);
            m_profStop(prof_(), "diagnostics");
        }

        //................................................
        // get the time-step given the field
        real_t dt = rk3.ComputeDt(advection, vel_);

        // dump some info
        m_log("RK3 - time = %f/%f - step %d/%d - dt = %e", t, t_final, iter, iter_max(), dt);

        //................................................
        // advance in time
        m_log("---- do time-step");
        m_profStart(prof_(), "do dt");
        rk3.DoDt(dt, &t);
        iter++;
        m_profStop(prof_(), "do dt");

        //................................................
        // diagnostics, dumps, whatever
        if (iter % iter_diag() == 0) {
            m_profStart(prof_(), "diagnostics");
            m_log("---- run diag");
            real_t time_now = MPI_Wtime();
            Diagnostics(t, dt, iter, time_now);
            m_profStop(prof_(), "diagnostics");
        }
    }
    m_profStop(prof_(), "run");
    // run the last diag
    if (iter % iter_diag() != 0) {
        real_t time_now = MPI_Wtime();
        Diagnostics(t, 0.0, iter, time_now);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SimpleAdvection::Diagnostics(const real_t time, const real_t dt, const lid_t iter, const real_t wtime) {
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

    // get fields moments
    real_t  moment0;
    real_t  moment1[3];
    BMoment moments;
    grid_->GhostPull(scal_);
    moments(grid_, scal_, &moment0, moment1);
    real_t          dmoment0;
    real_t          dmoment1[3];
    BDiscreteMoment dmoments;
    dmoments(grid_, scal_, &dmoment0, dmoment1);

    // compute the error
    real_t          err2, erri;
    Error  error;
    ring_->SetTime(time);
    (*ring_)(grid_, sol_);
    error.Norms(grid_, scal_, m_ptr<const Field>(sol_), &err2, &erri);

    // open the file
    FILE*   file_error;
    FILE*   file_diag;
    level_t min_level = grid_->MinLevel();
    level_t max_level = grid_->MaxLevel();
    if (rank == 0) {
        file_diag = fopen(string(folder_diag_ + "/diag_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + ".data").c_str(), "a+");
        // iter, time, dt, total quad, level min, level max
        fprintf(file_diag, "%6.6d;%e;%e;%ld;%d;%d", iter, time, dt, grid_->global_num_quadrants(), min_level, max_level);
        fprintf(file_diag, ";%e", wtime);
        fprintf(file_diag, ";%e;%e", grid_->rtol(), grid_->ctol());
        fprintf(file_diag, ";%e;%e", err2, erri);
        fprintf(file_diag, ";%e;%e;%e;%e", moment0, moment1[0], moment1[1], moment1[2]);
        fprintf(file_diag, ";%e;%e;%e;%e", dmoment0, dmoment1[0], dmoment1[1], dmoment1[2]);
        fprintf(file_diag, "\n");
        fclose(file_diag);
    }

    grid_->DumpLevels(iter, folder_diag_, string("_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT)));

    if (iter % iter_dump() == 0 && iter != 0) {
        // dump the vorticity field
        IOH5 dump(folder_diag_);
        grid_->GhostPull(scal_);
        dump(grid_(), scal_(), iter);
        dump(grid_(), sol_(), iter);

        // dump the details
        if (dump_detail()) {
            Field details("detail", scal_->lda());
            details.bctype(M_BC_EXTRAP);
            grid_->AddField(&details);
            grid_->StoreDetails(scal_, &details);

            grid_->GhostPull(&details);
            // IOH5 dump("data");
            dump(grid_(), &details, iter);

            grid_->DeleteField(&details);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
