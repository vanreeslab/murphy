#include "clients/weak_scalability.hpp"

#include "mpi.h"
#include "operator/advection.hpp"
#include "time/rk3_tvd.hpp"
#include "tools/ioh5.hpp"

using std::string;
using std::to_string;

static const lda_t  normal_tube = 0;
static const real_t sigma_tube  = 0.05;
static const real_t center[3]   = {0.5, 0.5, 0.5};
static const real_t velocity[3] = {0.0, 0.0, 1.0};

static const lambda_setvalue_t lambda_velocity = [](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
    m_assert(fid->lda() == 3, "the velocity field must be a vector");
    block->data(fid, 0)(i0, i1, i2) = velocity[0];
    block->data(fid, 1)(i0, i1, i2) = velocity[1];
    block->data(fid, 2)(i0, i1, i2) = velocity[2];
};

WeakScalability::~WeakScalability() {
    //-------------------------------------------------------------------------
    if (!(prof_ == nullptr)) {
        prof_->Disp();
        delete prof_;
    }
    //-------------------------------------------------------------------------
}

void WeakScalability::InitParam(ParserArguments* param) {
    m_begin;
    //--------------------------------------------------------------------------
    // call the general testcase parameters
    this->TestCase::InitParam(param);

    // the the standard stuffs
    if (param->profile) {
        // get the comm size
        rank_t comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("WeakScalability_") + to_string(comm_size) + string("ranks") + string("_w") + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT);
        prof_       = new Prof(name);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void WeakScalability::Run() {
    m_begin;
    //--------------------------------------------------------------------------
    m_profStart(prof_, "init");
    rank_t comm_size;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    // setup the grid
    bool   period[3] = {true, false, false};
    bidx_t length[3] = {m_max(1,comm_size/32), 1, 1};
    m_log("length is %d %d %d",length[0],length[1],length[2]);
    Grid   grid(1, period, length, MPI_COMM_WORLD, prof_);
    grid.level_limit(0, P8EST_QMAXLEVEL);

    //..........................................................................
    Field scal("scalar", 1);
    scal.bctype(M_BC_ZERO);
    grid.AddField(&scal);

    // setup the scalar tube
    lambda_setvalue_t lambda_tube = [=](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* const block, const Field* const fid) -> void {
        real_t pos[3];
        block->pos(i0, i1, i2, pos);
        block->data(fid, 0)(i0, i1, i2) = scalar_tube(pos, center, sigma_tube, normal_tube);
    };
    const bidx_t ghost_len_interp[2] = {m_max(grid.interp()->nghost_front(), 2),
                                        m_max(grid.interp()->nghost_back(), 2)};
    SetValue     tube(lambda_tube, ghost_len_interp);
    tube(&grid, &scal);

    // grid.GhostPull(&scal, ghost_len_ioh5);
    // IOH5 dump(folder_diag_);
    // dump(&grid, &scal, 0);

    // adapt the grid
    grid.SetTol(1e-5, 1e-7);
    grid.SetRecursiveAdapt(true);
    grid.Adapt(&scal, &tube);

    //..........................................................................
    // set the velocity field
    Field vel("velocity", 3);
    vel.bctype(M_BC_EXTRAP);
    vel.is_temp(true);
    grid.AddField(&vel);
    SetValue set_velocity(lambda_velocity, ghost_len_interp);
    set_velocity(&grid, &vel);

    m_profStop(prof_, "init");

    //..........................................................................
    // advection
    Advection<M_CONS, 3> adv_stencil(&vel, prof_);

    // time integration
    iter_t        iter = 0;
    real_t        t    = 0.0;
    const RK3_TVD rk3(&grid, &scal, &adv_stencil, prof_, 0.5);

    // let's gooo
    const iter_t iterfinal   = 50;
    const real_t tfinal      = 10000.0;
    const real_t wtime_start = MPI_Wtime();
    m_profStart(prof_, "run");
    while (t < tfinal && iter < iterfinal) {
        m_log("--------------------------------------------------------------------------------");
        //......................................................................
        // adapt the mesh
        if ((iter % iter_adapt() == 0)) {
            m_log("---- adapt mesh");
            m_log_level_plus;

            m_profStart(prof_, "adapt");
            grid.SetRecursiveAdapt(true);
            grid.Adapt(&scal);
            m_profStop(prof_, "adapt");

            // reset the velocity
            m_profStart(prof_, "set velocity");
            const bidx_t ghost_len_interp[2] = {m_max(grid.interp()->nghost_front(), 2),
                                                m_max(grid.interp()->nghost_back(), 2)};
            SetValue     set_velocity(lambda_velocity, ghost_len_interp);
            set_velocity(&grid,&vel);
            m_assert(vel.ghost_status(ghost_len_interp), "the velocity ghosts must have been computed");
            m_profStop(prof_, "set velocity");
            m_log_level_minus;
        }
        //......................................................................
        m_log("---- do time-step");
        m_log_level_plus;
        //......................................................................
        // get the time-step given the field
        m_profStart(prof_, "compute dt");
        real_t dt = rk3.ComputeDt(&adv_stencil, &vel);
        m_profStop(prof_, "compute dt");

        //......................................................................
        m_profStart(prof_, "dump IO");
        // get the rank
        rank_t rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        // get the info
        real_t  wtime_now       = MPI_Wtime();
        level_t min_level       = grid.MinLevel();
        level_t max_level       = grid.MaxLevel();
        long    global_num_quad = grid.global_num_quadrants();
        // dump
        if (rank == 0) {
            FILE*  file;
            string tag       = "w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + "_" + to_string(comm_size) + "ranks";
            string file_name = "diag_" + tag + ".data";
            file             = fopen(string(folder_diag_ + "/" + file_name).c_str(), "a+");
            fprintf(file, "%6.6d;%e;%ld;%ld;%d;%d;%e\n", iter, t, global_num_quad, global_num_quad / comm_size, min_level, max_level, wtime_now - wtime_start);
            fclose(file);
        }
        m_profStop(prof_, "dump IO");
        m_log("RK3 - time = %f/%f - step %d/%d - dt = %e - wtime = %e", t, tfinal, iter, iterfinal, dt, wtime_now - wtime_start);

        //......................................................................
        // advance in time
        m_profStart(prof_, "do dt");
        rk3.DoDt(dt, &t);
        iter++;
        m_profStop(prof_, "do dt");

        // grid.GhostPull(&scal, ghost_len_ioh5);
        // IOH5 dump(folder_diag_);
        // dump(&grid, &scal, iter);

        m_log_level_minus;
    }
    m_profStop(prof_, "run");
    const real_t wtime_final = MPI_Wtime();
    m_log("total wtime = %f", wtime_final - wtime_start);

    // delete the field
    // m_profStart(prof_, "cleanup");
    grid.DeleteField(&vel);
    grid.DeleteField(&scal);
    // m_profStop(prof_, "cleanup");

    //-------------------------------------------------------------------------
    m_end;
}
