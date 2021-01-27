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

    // the the standard stuffs
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("SimpleAdvection") + to_string(comm_size) + string("ranks");
        prof_.Alloc(name);
    }

    // setup the grid
    bool period[3] = {true, true, true};
    grid_.Alloc(param->init_lvl, period, param->length, MPI_COMM_WORLD, prof_);
    const real_t L[3] = {1.0, 1.0, 1.0};
    m_assert(grid_.IsOwned(), "the grid must be owned");

    // get the fields
    vel_.Alloc("velocity", 3);
    vel_->is_temp(true);
    scal_.Alloc("scalar", 1);
    grid_->AddField(vel_);
    grid_->AddField(scal_);

    // add the solution as temp
    sol_.Alloc("sol", 1);
    sol_->is_temp(true);
    grid_->AddField(sol_);

    // setup the scalar ring
    real_t velocity[3] = {0.0, 0.0, 1.0};
    ring_.Alloc(param->vr_normal, param->vr_center, param->vr_sigma, param->vr_radius, velocity, grid_->interp());
    ring_->Profile(prof_);
    ring_->SetTime(0.0);
    (*ring_)(grid_, scal_);

    // finish the grid
    grid_->SetTol(param->refine_tol, param->coarsen_tol);
    grid_->SetRecursiveAdapt(true);
    grid_->Adapt(scal_, ring_);

    // setup the velocity, 1.0 in every direction
    const lid_t  deg[3]   = {0, 0, 0};
    const real_t dir[3]   = {1.0, 0.0, 0.0};
    const real_t shift[3] = {0.0, 0.0, 0.0};
    vel_field_.Alloc(deg, dir, shift);
    (*vel_field_)(grid_, vel_, 2);
    // take the ghosts
    grid_->GhostPull(vel_);

    IOH5 dump(folder_diag_);
    dump(grid_(), vel_());
    //-------------------------------------------------------------------------
}

void SimpleAdvection::Run() {
    m_begin;
    //-------------------------------------------------------------------------
    // time
    lid_t  iter    = 0;
    real_t t_start = 0.0;
    real_t t_final = 1.0;
    real_t t       = 0.0;

    // test the moments with 1.0 in the z direction
    (*vel_field_)(grid_, vel_, 2);
    grid_->GhostPull(vel_);
    real_t  moment0[3];
    real_t  moment1[9];
    BMoment moments;
    moments(grid_, vel_, moment0, moment1);
    m_log("moments u = %e %e %e %e", moment0[0], moment1[0], moment1[1], moment1[2]);
    m_log("moments u = %e %e %e %e", moment0[1], moment1[3], moment1[4], moment1[5]);
    m_log("moments u = %e %e %e %e", moment0[2], moment1[6], moment1[7], moment1[8]);

    Advection<M_ADV_WENO_VEL, 3> adv(vel_);
    RK3_TVD                      rk3(grid_, scal_, &adv, prof_);
    adv.Profile(prof_);

    // let's gooo
    m_profStart(prof_(), "run");
    while (t < t_final && iter < iter_max()) {
        m_log("--------------------");
        //................................................
        // adapt the mesh
        if (iter % iter_adapt() == 0) {
            m_log("---- adapt mesh");
            m_profStart(prof_(), "adapt");
            grid_->Adapt(scal_);
            m_profStop(prof_(), "adapt");

            // reset the velocity
            (*vel_field_)(grid_, vel_, 2);
            grid_->GhostPull(vel_);
            m_assert(vel_->ghost_status(), "the velocity ghosts must have been computed");
        }
        // we run the first diagnostic
        if (iter == 0) {
            m_profStart(prof_(), "diagnostics");
            m_log("---- run diag");
            Diagnostics(t, 0, iter);
            m_profStop(prof_(), "diagnostics");
        }

        //................................................
        // get the time-step given the field
        real_t dt = rk3.ComputeDt(&adv, vel_);

        // dump some info
        m_log("RK3 - time = %f - step %d/%d - dt = %e", t, iter, iter_max(), dt);

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
            Diagnostics(t, dt, iter);
            m_profStop(prof_(), "diagnostics");
        }
    }
    m_profStop(prof_(), "run");
    // run the last diag
    if (iter % iter_diag() != 0) {
        Diagnostics(t, 0.0, iter);
    }
    //-------------------------------------------------------------------------
    m_end;
}

void SimpleAdvection::Diagnostics(const real_t time, const real_t dt, const lid_t iter) {
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

    // compute the error
    real_t          err2, erri;
    ErrorCalculator error;
    ring_->SetTime(time);
    (*ring_)(grid_, sol_);
    error.Norms(grid_, scal_, sol_, &err2, &erri);

    // open the file
    FILE*   file_error;
    FILE*   file_diag;
    level_t min_level = grid_->MinLevel();
    level_t max_level = grid_->MaxLevel();
    if (rank == 0) {
        file_diag = fopen(string(folder_diag_ + "/diag_w" + to_string(M_WAVELET_N) + to_string(M_WAVELET_NT) + ".data").c_str(), "a+");
        // iter, time, dt, total quad, level min, level max
        fprintf(file_diag, "%6.6d %e %e %ld %d %d", iter, time, dt, grid_->global_num_quadrants(), min_level, max_level);
        fprintf(file_diag, " %e %e", err2, erri);
        fprintf(file_diag, " %e %e %e %e", moment0, moment1[0], moment1[1], moment1[2]);
        fprintf(file_diag, "\n");
        fclose(file_diag);
    }

    if (iter % iter_dump() == 0) {
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
            grid_->DumpDetails(scal_, &details);

            grid_->GhostPull(&details);
            // IOH5 dump("data");
            dump(grid_(), &details, iter);

            grid_->DeleteField(&details);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
