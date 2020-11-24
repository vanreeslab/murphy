#include <argp.h>

#include <string>

#include "navier_stokes.hpp"
#include "parser.hpp"

using std::string;
using std::to_string;

NavierStokes::~NavierStokes() {
    //-------------------------------------------------------------------------
    delete (vort_);
    delete (grid_);
    delete (prof_);

    m_log("Navier Stokes is dead");
    //-------------------------------------------------------------------------
}

void NavierStokes::InitParam(ParserArguments* param) {
    m_begin;
    //-------------------------------------------------------------------------
    // setup the profiler
    if (param->profile) {
        int comm_size;
        MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
        string name = string("MURPHY_") + to_string(comm_size) + string("ranks");
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

    //-------------------------------------------------------------------------
    m_end;
}
