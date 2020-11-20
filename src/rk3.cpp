#include "rk3.hpp"

#include "blas.hpp"

/**
 * @brief Construct a new Runge Kutta 3 object
 * 
 * add a field to the grid: rk3
 * 
 * @param grid 
 */
RungeKutta3::RungeKutta3(Grid* grid, Field* state, Stencil* f) {
    m_begin;
    m_assert(grid != nullptr, "the grid cannot be null");
    m_assert(grid->IsAField(state), "the field must be part of the grid");
    //-------------------------------------------------------------------------
    grid_    = grid;
    field_u_ = state;
    f_       = f;  // store the stencil

    // add a temp field to the grid
    field_y_ = new Field("rk3_y", field_u_->lda());
    grid_->AddField(field_y_);
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
    m_assert(*time >= 0, "the time cannot be negative");
    m_assert(dt > 0, "the dt cannot be negative, nor 0");
    //-------------------------------------------------------------------------
    //................................................
    // step 1
    // y = 0.0 * y + F(t,u)
    // u = u + 1/3 * dt * y
    {
        // update the y using the stencil operator
        Scale scale(0.0);
        DoOpMesh(&scale, &Scale::ComputeScaleGridBlock, grid_, field_y_);
        (*f_)(grid_, field_u_, field_y_);

        // update the u
        Daxpy update_u(1.0 / 3.0 * dt);
        DoOpMesh(&update_u, &Daxpy::ComputeDaxpyGridBlock, grid_, field_y_, field_u_, field_u_);

        // udpate time
        (*time) += 1.0 / 3.0 * dt;
    }

    //................................................
    // step 2
    // y = -5/9 * y + F(t,u)
    // u = u + 15/16 * dt * y
    {
        // update the y using the stencil operator
        Scale scale(-5.0 / 9.0);
        DoOpMesh(&scale, &Scale::ComputeScaleGridBlock, grid_, field_y_);
        (*f_)(grid_, field_u_, field_y_);

        // update the u
        Daxpy update_u(15.0 / 16.0 * dt);
        DoOpMesh(&update_u, &Daxpy::ComputeDaxpyGridBlock, grid_, field_y_, field_u_, field_u_);

        // udpate time: 3/4 - 1/3 = 5/12
        (*time) += 5.0 / 12.0 * dt;
    }

    //................................................
    // step 3
    // y = -153/128 * y + F(t,u)
    // u = u + 8/15 * dt * y
    {
        // update the y using the stencil operator
        Scale scale(-153.0 / 128.0);
        DoOpMesh(&scale, &Scale::ComputeScaleGridBlock, grid_, field_y_);
        (*f_)(grid_, field_u_, field_y_);

        // update the u
        Daxpy update_u(8.0 / 15.0 * dt);
        DoOpMesh(&update_u, &Daxpy::ComputeDaxpyGridBlock, grid_, field_y_, field_u_, field_u_);

        // udpate time: 1 - 3/4 = 1/4
        (*time) += 1.0 / 4.0 * dt;
    }

    //-------------------------------------------------------------------------
}