#ifndef SRC_TIME_RKFUNCTOR_HPP
#define SRC_TIME_RKFUNCTOR_HPP

#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "grid/grid.hpp"

/**
 * @brief Functor that defines the interface used by the RK 
 * 
 */
class RKFunctor {
   public:
    // return stability info
    virtual real_t cfl_rk3() const   = 0;  //!< the CFL constrain
    virtual real_t rdiff() const = 0;  //!< the diffusion constrain

    // do the rhs
    virtual void RhsSet(m_ptr<const Grid> grid, const real_t time, m_ptr<Field> field_u, m_ptr<Field> field_y) = 0;
    virtual void RhsAcc(m_ptr<const Grid> grid, const real_t time, m_ptr<Field> field_u, m_ptr<Field> field_y) = 0;
};

#endif