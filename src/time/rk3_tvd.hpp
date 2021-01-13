#ifndef SRC_RK3_HPP
#define SRC_RK3_HPP

#include "core/macros.hpp"
#include "core/types.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "prof.hpp"
#include "stencil.hpp"

/**
 * @brief provides an implementation of a RK3 -TVD using two registers
 * 
 * For more detail, cfr
 *      Total Variation Diminishing Runge-Kutta Schemes, Gottlieb and Shu, 1998
 *      Strong Stability-Preserving High-Order Time Discretization Methods, Gottlieb et al., 2001
 * 
 */
class RK3_TVD {
    m_ptr<Grid>    grid_    = nullptr;
    m_ptr<Field>   field_u_ = nullptr;
    m_ptr<Stencil> f_       = nullptr;
    m_ptr<Prof>    prof_    = nullptr;

    Field* field_y1_ = nullptr;
    Field* field_y2_ = nullptr;

   public:
    explicit RK3_TVD(m_ptr<Grid> grid, m_ptr<Field> state, m_ptr<Stencil> f, m_ptr<Prof> prof);
    virtual ~RK3_TVD();

    void DoDt(const real_t dt, real_t* time);

    real_t ComputeDt(const real_t max_vel = 1.0);
};
#endif