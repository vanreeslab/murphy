#ifndef SRC_RK3_HPP_
#define SRC_RK3_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "prof.hpp"
#include "stencil.hpp"

/**
 * @brief provides an implementation of a RK3 low storage (cfr Williamson). The rhs is evaluated as a stencil
 * 
 * This low-storage scheme allows to only keep one state vector while proceeding to the time integration
 * 
 * u = u0
 * t = t0
 * 
 * step 1:
 *      y = F(t,u)
 *      u = u + 1/3 * dt * y
 *      t = t0 + 1/3 * dt
 * 
 * step 2:
 *      y = -5/9 * y + F(t,u)
 *      u = u + 15/16 * dt * y
 *      t = t0 + 3/4 * dt
 * 
 * step 3
 *      y = -153/128 * y + F(t,u)
 *      u = u + 8/15 * dt * y 
 *      t = t0 + dt
 * 
 */
class RungeKutta3 {
    Grid*    grid_    = nullptr;
    Field*   field_u_ = nullptr;
    Field*   field_y_ = nullptr;
    Stencil* f_       = nullptr;
    Prof*    prof_    = nullptr;

   public:
    explicit RungeKutta3(Grid* grid, Field* state, Stencil* f, Prof* prof);
    virtual ~RungeKutta3();

    void DoDt(const real_t dt, real_t* time);

    real_t ComputeDt();
};

#endif  // SRC_RK3_HPP_