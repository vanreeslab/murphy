#include "rk3.hpp"

#include "blas.hpp"
#include "ioh5.hpp"

/**
 * @brief Construct a new Runge Kutta 3 object
 * 
 * add a field to the grid: rk3
 * 
 * @param grid 
 */
RungeKutta3::RungeKutta3(Grid* grid, Field* state, Stencil* f, Prof* prof) {
    m_begin;
    m_assert(grid != nullptr, "the grid cannot be null");
    m_assert(grid->IsAField(state), "the field must be part of the grid");
    //-------------------------------------------------------------------------
    grid_    = grid;
    prof_    = prof;
    field_u_ = state;
    f_       = f;  // store the stencil

    // add a temp field to the grid
    field_y_ = new Field("rk3_y", field_u_->lda());
    grid_->AddField(field_y_);
    field_y_->is_temp(true);
    //-------------------------------------------------------------------------
    m_end;
}

RungeKutta3::~RungeKutta3() {
    m_begin;
    //-------------------------------------------------------------------------
    grid_->DeleteField(field_y_);
    delete (field_y_);
    //-------------------------------------------------------------------------
    m_end;
}

void RungeKutta3::DoDt(const real_t dt, real_t* time) {
    m_begin;
    m_assert(*time >= 0, "the time cannot be negative");
    m_assert(dt > 0.0, "the dt = %e cannot be negative, nor 0");
    //-------------------------------------------------------------------------

    // create the scale and the daxpy
    Scale scale;
    Daxpy daxpy;

    m_profStart(prof_, "rk3");
    //................................................
    // step 1
    // y = 0.0 * y + F(t,u)
    // u = u + 1/3 * dt * y

    m_profStart(prof_, "scale");
    scale(grid_, 0.0, field_y_);
    m_profStop(prof_, "scale");

    m_profStart(prof_, "rhs");
    (*f_)(grid_, field_u_, field_y_);
    m_profStop(prof_, "rhs");

    // IOH5 dump("data");
    // field_y_->bctype(M_BC_EXTRAP);
    // grid_->GhostPull(field_y_);
    // dump(grid_,field_y_);

    // update the u
    m_profStart(prof_, "update");
    daxpy(grid_, 1.0 / 3.0 * dt, field_y_, field_u_, field_u_);
    m_profStop(prof_, "update");

    // udpate time
    (*time) += 1.0 / 3.0 * dt;

    //................................................
    // step 2
    // y = -5/9 * y + F(t,u)
    // u = u + 15/16 * dt * y

    m_profStart(prof_, "scale");
    scale(grid_, -5.0 / 9.0, field_y_);
    m_profStop(prof_, "scale");

    m_profStart(prof_, "rhs");
    (*f_)(grid_, field_u_, field_y_);
    m_profStop(prof_, "rhs");

    // update the u
    m_profStart(prof_, "update");
    daxpy(grid_, 15.0 / 16.0 * dt, field_y_, field_u_, field_u_);
    m_profStop(prof_, "update");

    // udpate time: 3/4 - 1/3 = 5/12
    (*time) += 5.0 / 12.0 * dt;

    //................................................
    // step 3
    // y = -153/128 * y + F(t,u)
    // u = u + 8/15 * dt * y

    // update the y using the stencil operator
    m_profStart(prof_, "scale");
    scale(grid_, -153.0 / 128.0, field_y_);
    m_profStop(prof_, "scale");

    m_profStart(prof_, "rhs");
    (*f_)(grid_, field_u_, field_y_);
    m_profStop(prof_, "rhs");

    // update the u
    m_profStart(prof_, "update");
    daxpy(grid_, 8.0 / 15.0 * dt, field_y_, field_u_, field_u_);
    m_profStop(prof_, "update");

    // udpate time: 1 - 3/4 = 1/4
    (*time) += 1.0 / 4.0 * dt;

    m_profStop(prof_, "rk3");
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief compute the permited time-step
 * 
 * @return real_t 
 */
real_t RungeKutta3::ComputeDt() {
    // m_begin;
    //-------------------------------------------------------------------------
    // know the limits
    real_t cfl_limit = sqrt(3) * 0.98;  // CFL = max_vel * dt / h
    // real_t r_limit = 2.5; // r = nu * dt / h^2

    // get the finest h in the grid
    real_t h_fine = grid_->FinestH();
    m_assert(h_fine > 0.0, "the finest h = %e must be positive", h_fine);

    // get the fastest velocity
    real_t max_vel = 1.0;  //todo change
    real_t cfl_dt  = cfl_limit * h_fine / max_vel;
    m_assert(cfl_dt > 0.0, "the CFL dt = %e must be positive", cfl_dt);

    m_log("dt = %e, using h = %e and CFL limit = %e", cfl_dt, h_fine, cfl_limit);
    //-------------------------------------------------------------------------
    // m_end;
    return cfl_dt;
}