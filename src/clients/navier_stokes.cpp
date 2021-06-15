#include "navier_stokes.hpp"

#include <argp.h>

#include <string>

#include "advection_diffusion.hpp"
#include "conservative_advection_diffusion.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"
#include "error.hpp"
#include "ioh5.hpp"
#include "time/rk3_ls.hpp"

using std::string;
using std::to_string;

NavierStokes::~NavierStokes() {
    //-------------------------------------------------------------------------
    if (prof_ != nullptr) {
        prof_->Disp();
        delete (prof_);
    }

    delete (vort_);
    delete (grid_);

    m_log("Navier Stokes is dead");
    //-------------------------------------------------------------------------
}

void NavierStokes::InitParam(ParserArguments* param) {
    m_begin;
    //-------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    //store nu and ustream
    real_t L = 1.0;
    for (lda_t id = 0; id < 3; ++id) {
        u_stream_[id] = 0.0;
    }
    u_stream_[2] = 1.0;

    real_t nu_ = (param->reynolds > 0.0) ? (L * m_max(u_stream_[0], m_max(u_stream_[1], u_stream_[2])) / param->reynolds) : 0.0;

    dump_error_    = param->dump_error;
    compute_error_ = param->compute_error;

    // setup the profiler
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("Navier-Stokes_") + to_string(comm_size) + string("ranks");
        prof_       = new Prof(name);
    }

    // setup the grid
    grid_ = new Grid(param->init_lvl, param->period, param->length, MPI_COMM_WORLD, prof_);

    // the vorticity with the vortex ring
    vort_ = new Field("vorticity", 3);
    vort_->bctype(M_BC_EXTRAP);
    grid_->AddField(vort_);
    SetVortexRing vr_init(param->vr_normal, param->vr_center, param->vr_sigma, param->vr_radius, grid_->interp());
    vr_init.Profile(prof_);
    vr_init(grid_, vort_);

    // setup the grid with the given initial condition
    grid_->SetTol(param->refine_tol, param->coarsen_tol);
    grid_->SetRecursiveAdapt(true);
    grid_->Adapt(vort_, &vr_init);

    //-------------------------------------------------------------------------
    m_end;
}

void NavierStokes::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // time
    real_t t_start = 0.0;
    real_t t_final = 1.0;
    real_t t       = 0.0;
    // iterations
    lid_t iter = 0;

    // AdvectionDiffusion<5, 3> adv_diff(nu_, u_stream_);
    // Conservative_AdvectionDiffusion<4, 3> adv_diff(nu_, u_stream_);
    // adv_diff.Profile(prof_);

    // RK3_LS rk3(1.0 / 3.0, grid_, vort_, &adv_diff, prof_);
    

    // let's gooo
    m_profStart(prof_, "Navier-Stokes run");
    while (t < t_final && iter < iter_max()) {
        //................................................
        // get the time-step given the field
        // real_t dt = rk3.ComputeDt();
        real_t dt = 1e-4;

        // dump some info
        m_log("--------------------");
        m_log("RK3 - time = %f - step %d - dt = %e", t, iter, dt);

        //................................................
        // diagnostics, dumps, whatever
        if (iter % iter_diag() == 0) {
            m_profStart(prof_, "diagnostics");
            m_log("---- run diag");
            Diagnostics(t, dt, iter);
            m_profStop(prof_, "diagnostics");
        }

        //................................................
        // advance in time
        m_log("---- do time-step");
        m_profStart(prof_, "do dt");
        // rk3.DoDt(dt, &t);
        iter++;
        m_profStop(prof_, "do dt");

        //................................................
        // adapt the mesh
        if (iter % iter_adapt() == 0) {
            m_log("---- adapt mesh");
            m_profStart(prof_, "adapt");
            grid_->Adapt(vort_);
            m_profStop(prof_, "adapt");
        }
    }
    m_profStop(prof_, "Navier-Stokes run");
    // run the last diag
    Diagnostics(t, 0.0, iter);
    //-------------------------------------------------------------------------
    m_end;
}

void NavierStokes::Diagnostics(const real_t time, const real_t dt, const lid_t iter) {
    m_begin;
    //-------------------------------------------------------------------------
    rank_t rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // if the folder does not exist, create it
    struct stat st = {0};
    if (rank == 0 && stat(folder_diag_.c_str(), &st) == -1) {
        mkdir(folder_diag_.c_str(), 0770);
    }

    // open the file
    FILE* file_error;
    FILE* file_diag;
    level_t min_level = grid_->MinLevel();
    level_t max_level = grid_->MaxLevel();
    if (rank == 0) {
        file_diag = fopen(string(folder_diag_ + "/ns-diag.data").c_str(), "a+");

        // iter, time, dt, total quad, level min, level max
        fprintf(file_diag, "%6.6d %e %e %ld %d %d\n", iter, time, dt, grid_->global_num_quadrants(), min_level, max_level);
        fclose(file_diag);

        if (compute_error_) {
            file_error = fopen(string(folder_diag_ + "/ns-error.data").c_str(), "a+");
        }
    }

    // dump the vorticity field
    IOH5 dump(folder_diag_);
    grid_->GhostPull(vort_);
    dump(grid_, vort_, iter);

    // dump the details
    if (dump_detail()) {
        Field details("detail", 3);
        details.bctype(M_BC_EXTRAP);
        grid_->AddField(&details);
        grid_->StoreDetails(vort_, &details);

        grid_->GhostPull(&details);
        // IOH5 dump("data");
        dump(grid_, &details, iter);

        grid_->DeleteField(&details);
    }

    // compute the analytical solution
    if (compute_error_) {
        Field anal("anal", 3);
        Field err("error", 3);
        grid_->AddField(&anal);
        grid_->AddField(&err);
        real_t        center[3] = {0.5 + u_stream_[0] * time, 0.5 + u_stream_[1] * time, 0.5 + u_stream_[2] * time};
        real_t        radius    = 0.25;
        real_t        sigma     = 0.025 + sqrt(4.0 * nu_ * time);
        SetVortexRing vr_init(2, center, sigma, radius, grid_->interp());
        vr_init(grid_, &anal);

        grid_->GhostPull(vort_);
        // compute the error wrt to the analytical solution
        real_t          err2 = 0.0;
        real_t          erri = 0.0;
        Error  error(grid_->interp());
        error.Norms(grid_, vort_, &anal, &err, &err2, &erri);

        // I/O the error field
        if (dump_error_) {
            err.bctype(M_BC_EXTRAP);
            grid_->GhostPull(&err);
            dump(grid_, &err, iter);
        }

        // remove the added fields
        grid_->DeleteField(&anal);
        grid_->DeleteField(&err);

        m_log("error is norm 2: %e - norm inf: %e", err2, erri);
        if (rank == 0) {
            fprintf(file_error, "%6.6d %e %e %e\n", iter, time, err2, erri);
            fclose(file_error);
        }
    }

    //-------------------------------------------------------------------------
    m_end;
}
