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

//==============================================================================
template <>
inline lid_t Advection<M_ADV_CENTER, 2>::NGhost() const { return 1; };
template <>
inline void Advection<M_ADV_CENTER, 2>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    // -------------------------------------------------------------------------
    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();
    const real_t  oneoh[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    // get the alpha factor
    const real_t alpha = (accumulate_) ? 1.0 : 0.0;

    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);
        const real_t* src = data_src + m_idx(i0, i1, i2);
        const real_t* u   = data_u + m_idx(i0, i1, i2);
        const real_t* v   = data_v + m_idx(i0, i1, i2);
        const real_t* w   = data_w + m_idx(i0, i1, i2);

        // reset the value if needed
        trg[0] *= alpha;

        //advection
        trg[0] -= oneoh[0] * u[0] * 0.5 * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]);
        trg[0] -= oneoh[1] * v[0] * 0.5 * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]);
        trg[0] -= oneoh[2] * w[0] * 0.5 * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]);
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

        for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
        for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
        for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
};

//==============================================================================
template <>
inline lid_t Advection<M_ADV_CENTER, 4>::NGhost() const { return 2; };
template <>
inline void Advection<M_ADV_CENTER, 4>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    // -------------------------------------------------------------------------
    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();
    const real_t  oneoh[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    // get the alpha factor
    const real_t alpha = (accumulate_) ? 1.0 : 0.0;

    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);
        const real_t* src = data_src + m_idx(i0, i1, i2);
        const real_t* u   = data_u + m_idx(i0, i1, i2);
        const real_t* v   = data_v + m_idx(i0, i1, i2);
        const real_t* w   = data_w + m_idx(i0, i1, i2);

        constexpr real_t one_twelve = 1.0 / 12.0;
        constexpr real_t two_third  = 2.0 / 3.0;

        // reset the value if needed
        trg[0] *= alpha;

        //advection
        trg[0] -= oneoh[0] * u[0] * (-one_twelve * (src[m_idx(+2, 0, 0)] - src[m_idx(-2, 0, 0)]) + two_third * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
        trg[0] -= oneoh[1] * v[0] * (-one_twelve * (src[m_idx(0, +2, 0)] - src[m_idx(0, -2, 0)]) + two_third * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
        trg[0] -= oneoh[2] * w[0] * (-one_twelve * (src[m_idx(0, 0, +2)] - src[m_idx(0, 0, -2)]) + two_third * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

        for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
        for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
        for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
};

//==============================================================================
template <>
inline lid_t Advection<M_ADV_CENTER, 6>::NGhost() const { return 3; };
template <>
inline void Advection<M_ADV_CENTER, 6>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    // -------------------------------------------------------------------------
    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();
    const real_t  oneoh[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    // get the alpha factor
    const real_t alpha = (accumulate_) ? 1.0 : 0.0;

    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);
        const real_t* src = data_src + m_idx(i0, i1, i2);
        const real_t* u   = data_u + m_idx(i0, i1, i2);
        const real_t* v   = data_v + m_idx(i0, i1, i2);
        const real_t* w   = data_w + m_idx(i0, i1, i2);

        constexpr real_t one_sixtieth    = 1.0 / 60.0;
        constexpr real_t three_twentieth = 3.0 / 20.0;
        constexpr real_t three_forth     = 3.0 / 4.0;

        // reset the value if needed
        trg[0] *= alpha;

        //advection
        trg[0] -= oneoh[0] * u[0] * (+one_sixtieth * (src[m_idx(+3, 0, 0)] - src[m_idx(-3, 0, 0)]) - three_twentieth * (src[m_idx(+2, 0, 0)] - src[m_idx(-2, 0, 0)]) + three_forth * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
        trg[0] -= oneoh[1] * v[0] * (+one_sixtieth * (src[m_idx(0, +3, 0)] - src[m_idx(0, -3, 0)]) - three_twentieth * (src[m_idx(0, +2, 0)] - src[m_idx(0, -2, 0)]) + three_forth * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
        trg[0] -= oneoh[2] * w[0] * (+one_sixtieth * (src[m_idx(0, 0, +3)] - src[m_idx(0, 0, -3)]) - three_twentieth * (src[m_idx(0, 0, +2)] - src[m_idx(0, 0, -2)]) + three_forth * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

        for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
        for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
        for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
};

//==============================================================================
template <>
inline lid_t Advection<M_ADV_WENO_VEL, 3>::NGhost() const { return 2; };
template <>
inline void Advection<M_ADV_WENO_VEL, 3>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    m_assert(u_->ghost_status(), "the ghost values of the velocity must be known!");
    // -------------------------------------------------------------------------
    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();
    const real_t  oneoh[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    // get the alpha factor
    const real_t alpha = (accumulate_) ? 1.0 : 0.0;

    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);
        const real_t* src = data_src + m_idx(i0, i1, i2);
        const real_t* u   = data_u + m_idx(i0, i1, i2);
        const real_t* v   = data_v + m_idx(i0, i1, i2);
        const real_t* w   = data_w + m_idx(i0, i1, i2);

        // reset the value if needed
        trg[0] *= alpha;

        // compute the fluxes in X, Y and Z
        real_t flux[6];

        constexpr real_t eps = 1.0e-6;
        {  // X
            const real_t* vel     = u;
            const real_t  h_inv   = oneoh[0];
            const real_t  d1m2    = src[m_idx(-1, 0, 0)] * vel[m_idx(-1, 0, 0)] - src[m_idx(-2, 0, 0)] * vel[m_idx(-2, 0, 0)];
            const real_t  d0m1    = src[m_idx(+0, 0, 0)] * vel[m_idx(+0, 0, 0)] - src[m_idx(-1, 0, 0)] * vel[m_idx(-1, 0, 0)];
            const real_t  d1m0    = src[m_idx(+1, 0, 0)] * vel[m_idx(+1, 0, 0)] - src[m_idx(+0, 0, 0)] * vel[m_idx(+0, 0, 0)];
            const real_t  d2m1    = src[m_idx(+2, 0, 0)] * vel[m_idx(+2, 0, 0)] - src[m_idx(+1, 0, 0)] * vel[m_idx(+1, 0, 0)];
            const real_t  r_denom = 1.0 / (eps + pow(d1m0 - d0m1, 2));

            const real_t rm = (eps + pow(d0m1 - d1m2, 2)) * r_denom;
            const real_t wm = 1.0 / (1 + 2.0 * pow(rm, 2));
            flux[0]         = 0.5 * h_inv * ((d0m1 + d1m0) - wm * (d1m2 - 2.0 * d0m1 + d1m0));

            const real_t rp = (eps + pow(d2m1 - d1m0, 2)) * r_denom;
            const real_t wp = 1.0 / (1 + 2.0 * pow(rp, 2));
            flux[1]         = 0.5 * h_inv * ((d0m1 + d1m0) - wp * (d2m1 - 2.0 * d1m0 + d0m1));
        }
        {  // Y
            const real_t* vel     = v;
            const real_t  h_inv   = oneoh[1];
            const real_t  d1m2    = src[m_idx(0, -1, 0)] * vel[m_idx(0, -1, 0)] - src[m_idx(0, -2, 0)] * vel[m_idx(0, -2, 0)];
            const real_t  d0m1    = src[m_idx(0, +0, 0)] * vel[m_idx(0, +0, 0)] - src[m_idx(0, -1, 0)] * vel[m_idx(0, -1, 0)];
            const real_t  d1m0    = src[m_idx(0, +1, 0)] * vel[m_idx(0, +1, 0)] - src[m_idx(0, +0, 0)] * vel[m_idx(0, +0, 0)];
            const real_t  d2m1    = src[m_idx(0, +2, 0)] * vel[m_idx(0, +2, 0)] - src[m_idx(0, +1, 0)] * vel[m_idx(0, +1, 0)];
            const real_t  r_denom = 1.0 / (eps + pow(d1m0 - d0m1, 2));

            const real_t rm = (eps + pow(d0m1 - d1m2, 2)) * r_denom;
            const real_t wm = 1.0 / (1 + 2.0 * pow(rm, 2));
            flux[2]         = 0.5 * h_inv * ((d0m1 + d1m0) - wm * (d1m2 - 2.0 * d0m1 + d1m0));

            const real_t rp = (eps + pow(d2m1 - d1m0, 2)) * r_denom;
            const real_t wp = 1.0 / (1 + 2.0 * pow(rp, 2));
            flux[3]         = 0.5 * h_inv * ((d0m1 + d1m0) - wp * (d2m1 - 2.0 * d1m0 + d0m1));
        }
        {  // Z
            const real_t* vel     = w;
            const real_t  h_inv   = oneoh[2];
            const real_t  d1m2    = src[m_idx(0, 0, -1)] * vel[m_idx(0, 0, -1)] - src[m_idx(0, 0, -2)] * vel[m_idx(0, 0, -2)];
            const real_t  d0m1    = src[m_idx(0, 0, +0)] * vel[m_idx(0, 0, +0)] - src[m_idx(0, 0, -1)] * vel[m_idx(0, 0, -1)];
            const real_t  d1m0    = src[m_idx(0, 0, +1)] * vel[m_idx(0, 0, +1)] - src[m_idx(0, 0, +0)] * vel[m_idx(0, 0, +0)];
            const real_t  d2m1    = src[m_idx(0, 0, +2)] * vel[m_idx(0, 0, +2)] - src[m_idx(0, 0, +1)] * vel[m_idx(0, 0, +1)];
            const real_t  r_denom = 1.0 / (eps + pow(d1m0 - d0m1, 2));

            const real_t rm = (eps + pow(d0m1 - d1m2, 2)) * r_denom;
            const real_t wm = 1.0 / (1 + 2.0 * pow(rm, 2));
            flux[4]         = 0.5 * h_inv * ((d0m1 + d1m0) - wm * (d1m2 - 2.0 * d0m1 + d1m0));

            const real_t rp = (eps + pow(d2m1 - d1m0, 2)) * r_denom;
            const real_t wp = 1.0 / (1 + 2.0 * pow(rp, 2));
            flux[5]         = 0.5 * h_inv * ((d0m1 + d1m0) - wp * (d2m1 - 2.0 * d1m0 + d0m1));
        }
        // we use a simple upwind/downwind approach, also mentionned in Johnsen2006 as we have no pressure
        trg[0] -= 0.5 * (1.0 + m_sign(u[0])) * flux[0] + 0.5 * (1.0 - m_sign(u[0])) * flux[1];
        trg[0] -= 0.5 * (1.0 + m_sign(v[0])) * flux[2] + 0.5 * (1.0 - m_sign(v[0])) * flux[3];
        trg[0] -= 0.5 * (1.0 + m_sign(w[0])) * flux[4] + 0.5 * (1.0 - m_sign(w[0])) * flux[5];
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

        for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
        for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
        for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
};

//==============================================================================
// see Jiang 2000
static real_t PhiWeno(const real_t a, const real_t b, const real_t c, const real_t d) {
    constexpr real_t one_third = 1.0 / 3.0;
    constexpr real_t one_sixth = 1.0 / 6.0;
    constexpr real_t eps       = 1.0e-6;

    const real_t IS0 = 13.0 * pow(a - b, 2) + 3 * pow(a - 3.0 * b, 2);
    const real_t IS1 = 13.0 * pow(b - c, 2) + 3 * pow(b + c, 2);
    const real_t IS2 = 13.0 * pow(c - d, 2) + 3 * pow(3 * c - d, 2);
    const real_t a0  = 1.0 / pow(eps + IS0, 2);
    const real_t a1  = 6.0 / pow(eps + IS1, 2);
    const real_t a2  = 3.0 / pow(eps + IS2, 2);

    const real_t w_denom = (a0 + a1 + a2);
    const real_t w0      = a0 / w_denom;
    const real_t w2      = a2 / w_denom;
    const real_t phi     = one_third * w0 * (a - 2.0 * b + c) + one_sixth * (w2 - 0.5) * (b - 2.0 * c + d);
    return phi;
};

template <>
inline lid_t Advection<M_ADV_WENO_VEL, 5>::NGhost() const { return 3; };

/**
 * @brief WENO fith order implementation of the conservative advection equation, use a velocity based flux-splitting
 * 
 * This implementation relies on Jiang2000
 * 
 */
template <>
inline void Advection<M_ADV_WENO_VEL, 5>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    m_assert(u_->ghost_status(), "the ghost values of the velocity must be known!");
    // -------------------------------------------------------------------------
    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();
    const real_t  oneoh[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    // get the alpha factor
    const real_t alpha = (accumulate_) ? 1.0 : 0.0;

    auto op = [=, &data_trg](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);
        const real_t* src = data_src + m_idx(i0, i1, i2);
        const real_t* u   = data_u + m_idx(i0, i1, i2);
        const real_t* v   = data_v + m_idx(i0, i1, i2);
        const real_t* w   = data_w + m_idx(i0, i1, i2);

        // reset the value if needed
        trg[0] *= alpha;

        // compute the fluxes in X, Y and Z
        real_t flux[6];

        constexpr real_t one_twelve = 1.0 / 12.0;
        {  // X
            const real_t* vel      = u;
            const real_t  h_inv    = oneoh[0];
            const real_t  d2m3     = src[m_idx(-2, 0, 0)] * vel[m_idx(-2, 0, 0)] - src[m_idx(-3, 0, 0)] * vel[m_idx(-3, 0, 0)];
            const real_t  d1m2     = src[m_idx(-1, 0, 0)] * vel[m_idx(-1, 0, 0)] - src[m_idx(-2, 0, 0)] * vel[m_idx(-2, 0, 0)];
            const real_t  d0m1     = src[m_idx(+0, 0, 0)] * vel[m_idx(+0, 0, 0)] - src[m_idx(-1, 0, 0)] * vel[m_idx(-1, 0, 0)];
            const real_t  d1m0     = src[m_idx(+1, 0, 0)] * vel[m_idx(+1, 0, 0)] - src[m_idx(+0, 0, 0)] * vel[m_idx(+0, 0, 0)];
            const real_t  d2m1     = src[m_idx(+2, 0, 0)] * vel[m_idx(+2, 0, 0)] - src[m_idx(+1, 0, 0)] * vel[m_idx(+1, 0, 0)];
            const real_t  d3m2     = src[m_idx(+3, 0, 0)] * vel[m_idx(+3, 0, 0)] - src[m_idx(+2, 0, 0)] * vel[m_idx(+2, 0, 0)];
            const real_t  d_center = one_twelve * h_inv * (-d1m2 + 7.0 * d0m1 + 7 * d1m0 - d2m1);

            flux[0] = d_center - PhiWeno(/* a */ h_inv * (d1m2 - d2m3),
                                         /* b */ h_inv * (d0m1 - d1m2),
                                         /* c */ h_inv * (d1m0 - d0m1),
                                         /* d */ h_inv * (d2m1 - d1m0));

            flux[1] = d_center + PhiWeno(/* a */ h_inv * (d3m2 - d2m1),
                                         /* b */ h_inv * (d2m1 - d1m0),
                                         /* c */ h_inv * (d1m0 - d0m1),
                                         /* d */ h_inv * (d0m1 - d1m2));
        }
        {  // Y
            const real_t* vel      = v;
            const real_t  h_inv    = oneoh[1];
            const real_t  d2m3     = src[m_idx(0, -2, 0)] * vel[m_idx(0, -2, 0)] - src[m_idx(0, -3, 0)] * vel[m_idx(0, -3, 0)];
            const real_t  d1m2     = src[m_idx(0, -1, 0)] * vel[m_idx(0, -1, 0)] - src[m_idx(0, -2, 0)] * vel[m_idx(0, -2, 0)];
            const real_t  d0m1     = src[m_idx(0, +0, 0)] * vel[m_idx(0, +0, 0)] - src[m_idx(0, -1, 0)] * vel[m_idx(0, -1, 0)];
            const real_t  d1m0     = src[m_idx(0, +1, 0)] * vel[m_idx(0, +1, 0)] - src[m_idx(0, +0, 0)] * vel[m_idx(0, +0, 0)];
            const real_t  d2m1     = src[m_idx(0, +2, 0)] * vel[m_idx(0, +2, 0)] - src[m_idx(0, +1, 0)] * vel[m_idx(0, +1, 0)];
            const real_t  d3m2     = src[m_idx(0, +3, 0)] * vel[m_idx(0, +3, 0)] - src[m_idx(0, +2, 0)] * vel[m_idx(0, +2, 0)];
            const real_t  d_center = one_twelve * h_inv * (-d1m2 + 7.0 * d0m1 + 7 * d1m0 - d2m1);

            flux[2] = d_center - PhiWeno(/* a */ h_inv * (d1m2 - d2m3),
                                         /* b */ h_inv * (d0m1 - d1m2),
                                         /* c */ h_inv * (d1m0 - d0m1),
                                         /* d */ h_inv * (d2m1 - d1m0));

            flux[3] = d_center + PhiWeno(/* a */ h_inv * (d3m2 - d2m1),
                                         /* b */ h_inv * (d2m1 - d1m0),
                                         /* c */ h_inv * (d1m0 - d0m1),
                                         /* d */ h_inv * (d0m1 - d1m2));
        }
        {  // Z
            const real_t* vel      = w;
            const real_t  h_inv    = oneoh[2];
            const real_t  d2m3     = src[m_idx(0, 0, -2)] * vel[m_idx(0, 0, -2)] - src[m_idx(0, 0, -3)] * vel[m_idx(0, 0, -3)];
            const real_t  d1m2     = src[m_idx(0, 0, -1)] * vel[m_idx(0, 0, -1)] - src[m_idx(0, 0, -2)] * vel[m_idx(0, 0, -2)];
            const real_t  d0m1     = src[m_idx(0, 0, +0)] * vel[m_idx(0, 0, +0)] - src[m_idx(0, 0, -1)] * vel[m_idx(0, 0, -1)];
            const real_t  d1m0     = src[m_idx(0, 0, +1)] * vel[m_idx(0, 0, +1)] - src[m_idx(0, 0, +0)] * vel[m_idx(0, 0, +0)];
            const real_t  d2m1     = src[m_idx(0, 0, +2)] * vel[m_idx(0, 0, +2)] - src[m_idx(0, 0, +1)] * vel[m_idx(0, 0, +1)];
            const real_t  d3m2     = src[m_idx(0, 0, +3)] * vel[m_idx(0, 0, +3)] - src[m_idx(0, 0, +2)] * vel[m_idx(0, 0, +2)];
            const real_t  d_center = one_twelve * h_inv * (-d1m2 + 7.0 * d0m1 + 7 * d1m0 - d2m1);

            flux[4] = d_center - PhiWeno(/* a */ h_inv * (d1m2 - d2m3),
                                         /* b */ h_inv * (d0m1 - d1m2),
                                         /* c */ h_inv * (d1m0 - d0m1),
                                         /* d */ h_inv * (d2m1 - d1m0));

            flux[5] = d_center + PhiWeno(/* a */ h_inv * (d3m2 - d2m1),
                                         /* b */ h_inv * (d2m1 - d1m0),
                                         /* c */ h_inv * (d1m0 - d0m1),
                                         /* d */ h_inv * (d0m1 - d1m2));
        }
        // we use a simple upwind/downwind approach, also mentionned in Johnsen2006 as we have no pressure
        trg[0] -= 0.5 * (1.0 + m_sign(u[0])) * flux[0] + 0.5 * (1.0 - m_sign(u[0])) * flux[1];
        trg[0] -= 0.5 * (1.0 + m_sign(v[0])) * flux[2] + 0.5 * (1.0 - m_sign(v[0])) * flux[3];
        trg[0] -= 0.5 * (1.0 + m_sign(w[0])) * flux[4] + 0.5 * (1.0 - m_sign(w[0])) * flux[5];
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(&op);
    } else {
        // do the most on the X side to use vectorization
        for_loop<0, M_N, 0, M_N, 0, M_GS>(&op);          // Z-
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(&op);  // Z+

        for_loop<0, M_N, 0, M_GS, M_GS, M_N - M_GS>(&op);          // Y-
        for_loop<0, M_N, M_N - M_GS, M_N, M_GS, M_N - M_GS>(&op);  // Y+

        for_loop<0, M_GS, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);          // X-
        for_loop<M_N - M_GS, M_N, M_GS, M_N - M_GS, M_GS, M_N - M_GS>(&op);  // X+
    }
    // -------------------------------------------------------------------------
};

#endif  // SRC_ADVECTION_DIFFUSION_HPP_
