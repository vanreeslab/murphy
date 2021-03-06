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
    real_t     cfl_max_ = 1.0;
    Grid*      grid_    = nullptr;
    Field*     field_u_ = nullptr;
    RKFunctor* f_       = nullptr;
    Prof*      prof_    = nullptr;

    Field* field_y1_ = nullptr;
    Field* field_y2_ = nullptr;

   public:
    explicit RK3_TVD(Grid* grid, Field* state, RKFunctor* f, Prof* prof, const real_t cfl_max = 1.0);
    virtual ~RK3_TVD();

    void DoDt(const real_t dt, real_t* time) const;

    real_t ComputeDt(const RKFunctor* rhs, const real_t max_ve, const real_t nu = 0.0) const;
    real_t ComputeDt(const RKFunctor* rhs, const Field* velocity, const real_t nu = 0.0) const;
};
#endif