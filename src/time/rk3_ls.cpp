#include "rk3_ls.hpp"

#include "blas.hpp"
#include "ioh5.hpp"

/**
 * @brief Construct a new Runge Kutta 3 object
 * 
 * add a field to the grid: rk3
 * 
 * @param grid 
 */
RK3_LS::RK3_LS(const real_t c2, m_ptr<Grid> grid, m_ptr<Field> state, m_ptr<Stencil> f, m_ptr<Prof> prof) {
    m_begin;
    m_assert(!grid.IsEmpty(), "the grid cannot be null");
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

    // compute the magical coefficients, see Shu1997, eq. 4.17
    const real_t z1 = sqrt(36.0 * pow(c2, 4.0) + 36.0 * pow(c2, 3.0) - 135.0 * pow(c2, 2.0) + 84.0 * c2 - 12.0);
    const real_t z2 = 2.0 * pow(c2, 2.0) + c2 - 2.0;
    const real_t z3 = 12.0 * pow(c2, 4.0) - 18.0 * pow(c2, 3.0) + 18.0 * pow(c2, 2.0) - 11.0 * c2 + 2.0;
    const real_t z4 = 36.0 * pow(c2, 4.0) - 36.0 * pow(c2, 3.0) + 13.0 * pow(c2, 2.0) - 8.0 * c2 + 4.0;
    const real_t z5 = 69.0 * pow(c2, 3.0) - 62.0 * pow(c2, 2.0) + 28.0 * c2 - 8.0;
    const real_t z6 = 34.0 * pow(c2, 4.0) - 46.0 * pow(c2, 3.0) + 34.0 * pow(c2, 2.0) - 13.0 * c2 + 2.0;

    // m_log("zs: %f %f %f %f %f %f", z1, z2, z3, z4, z5, z6);

    b1_ = c2;
    b2_ = (12.0 * c2 * (c2 - 1.0) * (3.0 * z2 - z1) - pow(3.0 * z2 - z1, 2.0)) / (144.0 * c2 * (3.0 * c2 - 2.0) * pow(c2 - 1.0, 2.0));
    b3_ = (-24.0 * (3.0 * c2 - 2.0) * pow(c2 - 1.0, 2.0)) / (pow(3.0 * z2 - z1, 2.0) - 12.0 * c2 * (c2 - 1.0) * (3.0 * z2 - z1));
    a2_ = (-(6.0 * pow(c2, 2.0) - 4.0 * c2 + 1.0) * z1 + 3.0 * z3) / ((2.0 * c2 + 1.0) * z1 - 3.0 * (c2 + 2.0) * pow(2.0 * c2 - 1.0, 2.0));
    // a2_ =  ((2.0 * c2 + 1.0) * z1 - 3.0 * (c2 + 2.0) * pow(2.0 * c2 - 1.0, 2));
    a3_ = (-z4 * z1 + 108.0 * (2.0 * c2 - 1.0) * pow(c2, 5.0) - 3.0 * (2.0 * c2 - 1.0) * z5) / (24.0 * z1 * c2 * pow(c2 - 1.0, 4.0) + 72 * c2 * z6 + 72 * pow(c2, 6.0) * (2 * c2 - 13.0));

    // time coefficient obtained by doing the scheme with dT/dt = 1
    t1_ = b1_;
    t2_ = b2_ * (1.0 + a2_);
    t3_ = b3_ * (1.0 + a3_ * (1.0 + a2_));
    //-------------------------------------------------------------------------
    m_log("RK3 LS with b1 = %e, b2 = %e, b3 = %e -- a2 = %e, a3 = %e -- t1 = %f, t2 = %f, t3 = %f", b1_, b2_, b3_, a2_, a3_,t1_,t2_,t3_);
    m_end;
}

RK3_LS::~RK3_LS() {
    m_begin;
    //-------------------------------------------------------------------------
    grid_->DeleteField(field_y_);
    delete (field_y_);
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
void RK3_LS::DoDt(const real_t dt, real_t* time) {
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
    m_profStart(prof_(), "reset");
    reset(grid_(), 0.0, field_y_);
    m_profStop(prof_(), "reset");

    m_profStart(prof_(), "rhs");
    (*f_)(grid_, field_u_, field_y_);
    m_profStop(prof_(), "rhs");

    // update the u
    m_profStart(prof_(), "update");
    daxpy(grid_, b1_ * dt, field_y_, field_u_, field_u_);
    m_profStop(prof_(), "update");

    // udpate time
    (*time) += t1_ * dt;

    //................................................
    // step 2
    m_profStart(prof_(), "scale");
    scale(grid_(), a2_, field_y_);
    m_profStop(prof_(), "scale");

    m_profStart(prof_(), "rhs");
    (*f_)(grid_, field_u_, field_y_);
    m_profStop(prof_(), "rhs");

    // update the u
    m_profStart(prof_(), "update");
    daxpy(grid_, b2_ * dt, field_y_, field_u_, field_u_);
    m_profStop(prof_(), "update");

    // time update
    (*time) += t2_ * dt;

    //................................................
    // step 3
    // update the y using the stencil operator
    m_profStart(prof_(), "scale");
    scale(grid_(), a3_, field_y_);
    m_profStop(prof_(), "scale");

    m_profStart(prof_(), "rhs");
    (*f_)(grid_, field_u_, field_y_);
    m_profStop(prof_(), "rhs");

    // update the u
    m_profStart(prof_(), "update");
    daxpy(grid_, b3_ * dt, field_y_, field_u_, field_u_);
    m_profStop(prof_(), "update");

    // udpate time
    (*time) += t3_ * dt;

    m_profStop(prof_(), "rk3");
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief compute the permited time-step
 * 
 * @return real_t 
 */
real_t RK3_LS::ComputeDt(const real_t max_vel) {
    // m_begin;
    //-------------------------------------------------------------------------
    // know the limits
    real_t cfl_limit = 1;  // CFL = max_vel * dt / h
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