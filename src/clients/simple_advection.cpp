#include "clients/simple_advection.hpp"

using std::string;
using std::to_string;

SimpleAdvection::~SimpleAdvection() {
    //-------------------------------------------------------------------------
    if (!prof_.IsEmpty()) {
        prof_->Disp();
        prof_.Free();
    }

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
    scal_.Alloc("scalar", 1);
    grid_->AddField(vel_);
    grid_->AddField(scal_);

    // setup the scalar ring
    SetScalarRing flow_ring(param->vr_normal, param->vr_center, param->vr_sigma, param->vr_radius, grid_->interp());
    flow_ring.Profile(prof_);
    flow_ring(grid_, scal_);

    // setup the velocity, 1.0 in every direction
    const lid_t  deg[3]   = {0, 0, 0};
    const real_t dir[3]   = {1.0, 0.0, 0.0};
    const real_t shift[3] = {0.0, 0.0, 0.0};
    SetPolynom   vel_init(deg, dir, shift, grid_->interp());
    vel_init(grid_, vel_);

    // finish the grid
    grid_->SetTol(param->refine_tol, param->coarsen_tol);
    grid_->SetRecursiveAdapt(true);
    grid_->Adapt(scal_, &flow_ring);

    //-------------------------------------------------------------------------
}

void SimpleAdvection::Run() {
}

void SimpleAdvection::Diagnostics(const real_t time, const real_t dt, const lid_t iter) {
}
