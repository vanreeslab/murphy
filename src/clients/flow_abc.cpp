#include "flow_abc.hpp"

#include "core/forloop.hpp"
#include "operator/advection.hpp"
#include "rk3.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

FlowABC::~FlowABC() {
    //-------------------------------------------------------------------------
    if (prof_ != nullptr) {
        prof_->Disp();
        delete (prof_);
    }

    delete (vel_);
    delete (scal_);
    delete (grid_);

    m_log("Navier Stokes is dead");
    //-------------------------------------------------------------------------
}

void FlowABC::InitParam(ParserArguments* param) {
    //-------------------------------------------------------------------------
    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("Navier-Stokes_") + to_string(comm_size) + string("ranks");
        prof_       = new Prof(name);
    }
    iter_max_   = param->iter_max;
    iter_diag_  = param->iter_diag;
    iter_adapt_ = param->iter_adapt;
    iter_dump_  = param->iter_dump;

    // setup the grid
    bool period[3]    = {true, true, true};
    grid_             = new Grid(param->init_lvl, period, param->length, MPI_COMM_WORLD, prof_);
    const real_t L[3] = {1.0, 1.0, 1.0};

    // get the fields
    vel_  = new Field("velocity", 3);
    scal_ = new Field("scalar", 1);
    grid_->AddField(vel_);
    grid_->AddField(scal_);

    // setup the flow ring
    SetScalarTube flow_ring(param->vr_normal, param->vr_center, param->vr_sigma, param->vr_radius, grid_->interp());
    flow_ring.Profile(prof_);
    flow_ring(grid_, scal_);

    grid_->SetTol(param->refine_tol, param->coarsen_tol);
    grid_->SetRecursiveAdapt(true);
    grid_->Adapt(scal_, &flow_ring);
    //-------------------------------------------------------------------------
}

void FlowABC::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // time
    lid_t  iter    = 0;
    real_t t_start = 0.0;
    real_t t_final = 1.0;
    real_t t       = 0.0;

    SetABSVelocity flow_vel(1.0, 0.5, 0.25, grid_->interp());

    ConsAdvection<4> adv(vel_);
    RungeKutta3      rk3(grid_, scal_, &adv, prof_);
    adv.Profile(prof_);

    // let's gooo
    m_profStart(prof_, "ABS Flow run");
    while (t < t_final && iter < iter_max_) {
        m_log("--------------------");
        //................................................
        // adapt the mesh
        if (iter % iter_adapt_ == 0) {
            m_log("---- adapt mesh");
            m_profStart(prof_, "adapt");
            grid_->Adapt(scal_);
            m_profStop(prof_, "adapt");

            // reset the velocity
            flow_vel(grid_, vel_);
            m_assert(vel_->ghost_status(), "the velocity ghosts must have been computed");
        }

        //................................................
        // get the time-step given the field
        real_t dt = rk3.ComputeDt(2.0);

        // dump some info
        m_log("RK3 - time = %f - step %d - dt = %e", t, iter, dt);

        //................................................
        // advance in time
        m_log("---- do time-step");
        m_profStart(prof_, "do dt");
        rk3.DoDt(dt, &t);
        iter++;
        m_profStop(prof_, "do dt");

        //................................................
        // diagnostics, dumps, whatever
        if (iter % iter_diag_ == 0) {
            m_profStart(prof_, "diagnostics");
            m_log("---- run diag");
            Diagnostics(t, dt, iter);
            m_profStop(prof_, "diagnostics");
        }
    }
    m_profStop(prof_, "ABS Flow run");
    // run the last diag
    Diagnostics(t, 0.0, iter);
    //-------------------------------------------------------------------------
    m_end;
}

void FlowABC::Diagnostics(const real_t time, const real_t dt, const lid_t iter) {
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
    FILE*   file_error;
    FILE*   file_diag;
    level_t min_level = grid_->MinLevel();
    level_t max_level = grid_->MaxLevel();
    if (rank == 0) {
        file_diag = fopen(string(folder_diag_ + "/ns-diag.data").c_str(), "a+");

        // iter, time, dt, total quad, level min, level max
        fprintf(file_diag, "%6.6d %e %e %ld %d %d\n", iter, time, dt, grid_->global_num_quadrants(), min_level, max_level);
        fclose(file_diag);
    }

    if (iter % iter_dump_ == 0) {
        // dump the vorticity field
        IOH5 dump(folder_diag_);
        grid_->GhostPull(scal_);
        grid_->GhostPull(vel_);
        dump(grid_, scal_, iter);
        dump(grid_, vel_, iter);
    }
    //-------------------------------------------------------------------------
    m_end;
}
