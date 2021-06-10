#ifndef SRC_OPERATOR_ADVECTION_HPP_
#define SRC_OPERATOR_ADVECTION_HPP_

#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "operator/stencil.hpp"
#include "time/rkfunctor.hpp"

typedef enum AdvectionType_t {
    M_WENO_Z
} AdvectionType_t;

/**
 * @brief compute the advection stencil on a given field using centered FD
 *
 * @tparam length the size of the advection stencil. (e.g. 3 for 2nd order and 5 for 4th order)
 */
template <AdvectionType_t type, short_t order>
class Advection : public Stencil, public RKFunctor {
   protected:
    bool               accumulate_ = false;    //!<  used to determine if we accumuate or not the solution
    m_ptr<const Field> u_          = nullptr;  //!< velocity field used for the advection

   public:
    // create the advection, just store the u field
    explicit Advection(m_ptr<const Field> u) : u_(u), Stencil(), RKFunctor(){};
    ~Advection() = default;

    // return the stability conditions
    real_t cfl_rk3() const override { return 1.0; };
    real_t rdiff() const override { return std::numeric_limits<real_t>::max(); };  // no limit, return +inf

    // return the number of ghost points, needed, depend on the stencil etc
    lid_t NGhost() const override { return 0; };

    void RhsSet(m_ptr<const Grid> grid, const real_t time, m_ptr<Field> field_u, m_ptr<Field> field_y) override {
        // -------------------------------------------------------------------------
        accumulate_ = false;
        Stencil::operator()(grid, field_u, field_y);
        accumulate_ = false;  // reset the value to false for other calls
        // -------------------------------------------------------------------------
    };
    void RhsAcc(m_ptr<const Grid> grid, const real_t time, m_ptr<Field> field_u, m_ptr<Field> field_y) override {
        // -------------------------------------------------------------------------
        accumulate_ = true;
        Stencil::operator()(grid, field_u, field_y);
        accumulate_ = false;  // reset the value to false for other calls
        // -------------------------------------------------------------------------
    };

   protected:
    void DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const override {
        m_assert(false, "This combination of type and order does not exist!");
    };
};

// //==============================================================================
// List the specializations
template <>
inline real_t Advection<M_WENO_Z, 3>::cfl_rk3() const {
    // comes from the "worst case scenario" von neumann analysis on the stencil, i.e. the 5th order stencil.
    // we assume that the smoothing will only introduce diffusion and hence will ease the CLF constraint
    return 1.6;
};
template <>
inline lid_t Advection<M_WENO_Z, 3>::NGhost() const {
    // we need one point outside the domain that need 1 ghost point
    return 2;
};
template <>
void Advection<M_WENO_Z, 3>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

template <>
inline real_t Advection<M_WENO_Z, 5>::cfl_rk3() const {
    // comes from the "worst case scenario" von neumann analysis on the stencil, i.e. the 5th order stencil.
    // we assume that the smoothing will only introduce diffusion and hence will ease the CLF constraint
    return 1.4;
};
template <>
inline lid_t Advection<M_WENO_Z, 5>::NGhost() const {
    // we need one point outside the domain that need 2 ghost point
    return 3;
};
template <>
void Advection<M_WENO_Z, 5>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

#endif  // SRC_ADVECTION_DIFFUSION_HPP_
