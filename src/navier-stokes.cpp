#include <argp.h>

#include <string>

#include "advection_diffusion.hpp"
#include "defs.hpp"
#include "ioh5.hpp"
#include "navier_stokes.hpp"
#include "rk3.hpp"
#include "error.hpp"

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
    reynolds_ = param->reynolds;
    for (lda_t id = 0; id < 3; ++id) {
        u_stream_[id] = 0.0;
    }
    u_stream_[2] = 1.0;

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

    // dump the field
    grid_->GhostPull(vort_);
    IOH5 dump("data");
    dump(grid_, vort_);
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
    lid_t iter     = 0;
    lid_t iter_max = 100;

    // Advection/diffusion - rk3
    real_t nu = (reynolds_ > 0.0) ? (1.0 / reynolds_) : 0.0;

    AdvectionDiffusion<5, 3> adv_diff(nu, u_stream_);
    RungeKutta3              rk3(grid_, vort_, &adv_diff, prof_);

    // the IO for the initial condition
    IOH5 dump("data");
    {
        char name[128];
        sprintf(name, "vort_%5.5d", iter);
        grid_->GhostPull(vort_);
        dump(grid_, vort_, name);
    }

    // let's gooo
    while (t < t_final && iter < iter_max) {
        //................................................
        // get the time-step given the field
        real_t dt = rk3.ComputeDt();

        // dump some info
        m_log("--------------------");
        m_log("RK3 - time = %f - step %d - dt = %e", t, iter, dt);

        //................................................
        // advance in time
        m_log("---- do time-step");
        rk3.DoDt(dt, &t);
        iter++;

        // //................................................
        // // diagnostics, dumps, whatever
        // if (iter % 10 == 0) {
        //     char name[128];
        //     sprintf(name, "vort_%5.5d", iter);
        //     grid_->GhostPull(vort_);
        //     dump(grid_, vort_, name);
        // }
        // adapt the mesh
        // if (iter % (M_N / 2) == 0) {
        if (iter % 1 == 0) {
            m_log("---- adapt mesh");

            char name[128];
            // sprintf(name, "vort_ba_%5.5d", iter);
            // grid_->GhostPull(vort_);
            // dump(grid_, vort_, name);

            {
                Field anal("anal",3);
                Field err("error",3);
                grid_->AddField(&anal);
                grid_->AddField(&err);
                real_t        center[3] = {0.5, 0.5, 0.5 + t};
                SetVortexRing vr_init(2, center, 0.05, 0.25, grid_->interp());
                vr_init(grid_, &anal);

                // compute the error wrt to the analytical solution
                real_t          err2 = 0.0;
                real_t          erri = 0.0;
                ErrorCalculator error;
                error.Norms(grid_, vort_, &anal,&err, &err2, &erri);

                err.bctype(M_BC_EXTRAP);
                grid_->GhostPull(&err);
                sprintf(name, "err_%5.5d", iter);
                dump(grid_, &err, name);

                grid_->DeleteField(&anal);
                grid_->DeleteField(&err);

                m_log("error is norm 2: %e - norm inf: %e", err2, erri);
            }

            grid_->Adapt(vort_);

            // sprintf(name, "vort_aa_%5.5d", iter);
            // grid_->GhostPull(vort_);
            // dump(grid_, vort_, name);
        }
    }
    //-------------------------------------------------------------------------
    m_end;
}
