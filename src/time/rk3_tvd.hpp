#ifndef SRC_RK3_HPP
#define SRC_RK3_HPP

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "tools/prof.hpp"
#include "rkfunctor.hpp"

/**
 * @brief provides an implementation of a RK3 -TVD using two registers
 * 
 * For more detail, cfr
 *      Total Variation Diminishing Runge-Kutta Schemes, Gottlieb and Shu, 1998
 *      Strong Stability-Preserving High-Order Time Discretization Methods, Gottlieb et al., 2001
 * 
 */
class RK3_TVD {
    real_t cfl_max_ = 1.0;

    m_ptr<Grid>      grid_    = nullptr;
    m_ptr<Field>     field_u_ = nullptr;
    m_ptr<RKFunctor> f_       = nullptr;
    m_ptr<Prof>      prof_    = nullptr;

    Field* field_y1_ = nullptr;
    Field* field_y2_ = nullptr;

   public:
    explicit RK3_TVD(m_ptr<Grid> grid, m_ptr<Field> state, m_ptr<RKFunctor> f, m_ptr<Prof> prof, const real_t cfl_max = 1.0);
    virtual ~RK3_TVD();

    void DoDt(const real_t dt, real_t* time);

    real_t ComputeDt(m_ptr<const RKFunctor> rhs, m_ptr<const Field> velocity) const;
};
#endif