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
 * This low-storage scheme allows to only keep one state vector while proceeding to the time integration.
 * This class is based on Shu 1997
 * 
 * u = u0
 * t = t0
 * 
 * step 1:
 *      y = F(t,u)
 *      u = u + B1 * dt * y
 *      t = t0 + B1 * dt
 * 
 * step 2:
 *      y = A2 * y + F(t,u)
 *      u = u + B2 * dt * y
 *      t = t1 + B2 * (1 + A2) * dt
 * 
 * step 3
 *      y = A3 * y + F(t,u)
 *      u = u + B3 * dt * y 
 *      t = t0 + B3 * (1 + A3 * (1 + A2)) * dt
 * 
 */
class RungeKutta3 {
    real_t b1_, a2_, b2_, a3_, b3_;
    real_t t1_,t2_,t3_;

    m_ptr<Grid>    grid_    = nullptr;
    m_ptr<Field>   field_u_ = nullptr;
    m_ptr<Stencil> f_       = nullptr;
    m_ptr<Prof>    prof_    = nullptr;

    Field* field_y_ = nullptr;

   public:
    explicit RungeKutta3(const real_t c2, m_ptr<Grid> grid, m_ptr<Field> state, m_ptr<Stencil> f, m_ptr<Prof> prof);
    virtual ~RungeKutta3();

    void DoDt(const real_t dt, real_t* time);

    real_t ComputeDt(const real_t max_vel = 1.0);
};

#endif  // SRC_RK3_HPP_