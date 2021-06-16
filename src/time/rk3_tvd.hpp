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
    Grid*       grid_    = nullptr;
    Field*      field_u_ = nullptr;
    RKFunctor*  f_       = nullptr;
    Prof*       prof_    = nullptr;

    Field* field_y1_ = nullptr;
    Field* field_y2_ = nullptr;

   public:
    explicit RK3_TVD(Grid*  grid, Field*  state, RKFunctor*  f, Prof*  prof);
    virtual ~RK3_TVD();

    void DoDt(const real_t dt, real_t* time);

    real_t ComputeDt(const RKFunctor*  rhs, const Field*  velocity) const;
};
#endif