#include "rk3_tvd.hpp"

#include "blas.hpp"
#include "ioh5.hpp"

/**
 * @brief Construct a new Runge Kutta 3 object
 * 
 * add a field to the grid: rk3
 * 
 * @param grid 
 */
RK3_TVD::RK3_TVD(m_ptr<Grid> grid, m_ptr<Field> state, m_ptr<Stencil> f, m_ptr<Prof> prof) {
    m_begin;
    m_assert(!grid.IsEmpty(), "the grid cannot be null");
    m_assert(grid->IsAField(state), "the field must be part of the grid");
    //-------------------------------------------------------------------------
    grid_    = grid;
    prof_    = prof;
    field_u_ = state;
    f_       = f;  // store the stencil

    // add a temp field to the grid
    field_y1_ = new Field("rk3_y1", field_u_->lda());
    field_y2_ = new Field("rk3_y2", field_u_->lda());
    grid_->AddField(field_y1_);
    grid_->AddField(field_y2_);
    field_y1_->is_temp(true);
    field_y2_->is_temp(true);

    //-------------------------------------------------------------------------
    m_end;
}

RK3_TVD::~RK3_TVD() {
    m_begin;
    //-------------------------------------------------------------------------
    grid_->DeleteField(field_y1_);
    grid_->DeleteField(field_y2_);
    delete (field_y1_);
    delete (field_y2_);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief process one time step of the rk3 method.
 * 
 * The coeficients are obtained in the litterature.
 * The time steps are obtained by integrating du/dt = 1 using the chosen scheme.
 * 
 * @param dt 
 * @param time 
 */
void RK3_TVD::DoDt(const real_t dt, real_t* time) {
    m_begin;
    m_assert(*time >= 0, "the time cannot be negative");
    m_assert(dt > 0.0, "the dt = %e cannot be negative, nor 0", dt);
    //-------------------------------------------------------------------------

    // create the scale and the daxpy
    Dset   reset;
    Dscale scale;
    Daxpy  daxpy;

    m_profStart(prof_(), "rk3");
    //................................................
    // step 1

    // y1 = 0.0
    m_profStart(prof_(), "reset");
    reset(grid_(), 0.0, field_y1_);
    m_profStop(prof_(), "reset");

    // y1 = 0.0 + f(u)
    m_profStart(prof_(), "rhs");
    (*f_)(grid_, field_u_, field_y1_);
    m_profStop(prof_(), "rhs");

    // y1 = dt*y1 + u -> y1 = u + dt * f(u)
    m_profStart(prof_(), "update y1");
    daxpy(grid_, dt, field_y1_, field_u_, field_y1_);
    m_profStop(prof_(), "update y1");

    //................................................
    // step 2

    // y2 = 0.0
    m_profStart(prof_(), "reset");
    reset(grid_(), 0.0, field_y2_);
    m_profStop(prof_(), "reset");

    // y2 = y2 + f(y1)
    m_profStart(prof_(), "rhs");
    (*f_)(grid_, field_y1_, field_y2_);
    m_profStop(prof_(), "rhs");

    // y2 = y2 * 1/4 * dt
    m_profStart(prof_(), "scale");
    scale(grid_(), 1.0 / 4.0 * dt, field_y2_);
    m_profStop(prof_(), "scale");

    // y2 = y2 + 1/4 * y1 + 3/4 * u -> y2 = 3/4*u + 1/4*(y1 + dt * f(y1))
    m_profStart(prof_(), "update");
    daxpy(grid_, 1.0 / 4.0, field_y1_, field_y2_, field_y2_);
    daxpy(grid_, 3.0 / 4.0, field_u_, field_y2_, field_y2_);
    m_profStop(prof_(), "update");

    //................................................
    // step 3
    // y1 = 0.0
    m_profStart(prof_(), "reset");
    reset(grid_(), 0.0, field_y1_);
    m_profStop(prof_(), "reset");

    // y1 = y1 + f(y2) = f(y2)
    m_profStart(prof_(), "rhs");
    (*f_)(grid_, field_y2_, field_y1_);
    m_profStop(prof_(), "rhs");

    // u = u * 1/3
    m_profStart(prof_(), "scale");
    scale(grid_(), 1.0 / 3.0, field_u_);
    m_profStop(prof_(), "scale");

    // u = u + 2/3 * y2 + 2/3 * dt * y1 -> u = 1/3 u + 2/3 * (y2 + dt * f(y2))
    m_profStart(prof_(), "update");
    daxpy(grid_, 2.0 / 3.0, field_y2_, field_u_, field_u_);
    daxpy(grid_, 2.0 / 3.0 * dt, field_y1_, field_u_, field_u_);
    m_profStop(prof_(), "update");

    //................................................
    (*time) += dt;

    m_profStop(prof_(), "rk3");
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief compute the permited time-step
 * 
 * @return real_t 
 */
real_t RK3_TVD::ComputeDt(const real_t max_vel) {
    // m_begin;
    //-------------------------------------------------------------------------
    // know the limits
    real_t cfl_limit = 0.8;  // CFL = max_vel * dt / h
    // real_t r_limit = 2.5; // r = nu * dt / h^2

    // get the finest h in the grid
    real_t h_fine = grid_->FinestH();
    m_assert(h_fine > 0.0, "the finest h = %e must be positive", h_fine);

    // get the fastest velocity
    real_t cfl_dt = cfl_limit * h_fine / max_vel;
    m_assert(cfl_dt > 0.0, "the CFL dt = %e must be positive", cfl_dt);

    m_log("dt = %e, using h = %e and CFL limit = %e", cfl_dt, h_fine, cfl_limit);
    //-------------------------------------------------------------------------
    // m_end;
    return cfl_dt;
}