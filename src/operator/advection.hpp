#ifndef SRC_OPERATOR_ADVECTION_HPP_
#define SRC_OPERATOR_ADVECTION_HPP_

#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "operator/stencil.hpp"
#include "time/rkfunctor.hpp"

using AdvectionType = enum {
    M_ADV_CENTER,
    M_ADV_WENO_VEL
};

/**
 * @brief compute the advection stencil on a given field using centered FD
 *
 * @tparam length the size of the advection stencil. (e.g. 3 for 2nd order and 5 for 4th order)
 */
template <AdvectionType type, short_t order>
class Advection : public Stencil, public RKFunctor {
   protected:
    bool               accumulate_ = false;    //!<  used to determine if we accumuate or not the solution
    m_ptr<const Field> u_          = nullptr;  //!< velocity field used for the advection

   public:
    // create the advection, just store the u field
    explicit Advection(m_ptr<const Field> u) : u_(u), Stencil(), RKFunctor(){};

    // return the stability conditions
    real_t cfl() const override { return 1.0; };
    real_t rdiff() const override { return 100.0; };  // std::numeric_limits<real_t>::max(); };  // no limit, return +inf

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

// only delcare the specialization but define them in the cpp (ease readying and reduce compilation time)
//==============================================================================
// CENTER - 2nd order
template <>
lid_t Advection<M_ADV_CENTER, 2>::NGhost() const;
template <>
void Advection<M_ADV_CENTER, 2>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

//==============================================================================
// CENTER - 4nd order
template <>
lid_t Advection<M_ADV_CENTER, 4>::NGhost() const;
template <>
void Advection<M_ADV_CENTER, 4>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

//==============================================================================
// CENTER - 6nd order
template <>
lid_t Advection<M_ADV_CENTER, 6>::NGhost() const;
template <>
void Advection<M_ADV_CENTER, 6>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

//==============================================================================
// WENO - 3rd order
template <>
lid_t Advection<M_ADV_WENO_VEL, 3>::NGhost() const;
template <>
void Advection<M_ADV_WENO_VEL, 3>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

//==============================================================================
// WENO - 5th order
template <>
lid_t Advection<M_ADV_WENO_VEL, 5>::NGhost() const;
template <>
void Advection<M_ADV_WENO_VEL, 5>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const;

#endif  // SRC_ADVECTION_DIFFUSION_HPP_
