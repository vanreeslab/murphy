#ifndef SRC_OPERATOR_ADVECTION_HPP_
#define SRC_OPERATOR_ADVECTION_HPP_

#include "core/forloop.hpp"
#include "core/macros.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "operator/stencil.hpp"

/**
 * @brief compute the advection stencil on a given field using centered FD
 * 
 * @tparam length the size of the advection stencil. (e.g. 3 for 2nd order and 5 for 4th order)
 */
template <sid_t length>
class Advection : public Stencil {
   protected:
    m_ptr<const Field> u_;

   public:
    explicit Advection(m_ptr<const Field> u) : u_(u), Stencil(){};
    virtual lid_t NGhost() const override { return (length / 2); };

   protected:
    void DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) override;
};

template <sid_t length>
void Advection<length>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) {
    m_assert(fid_src->lda() == fid_trg->lda(), "the source and target field must have the same dimension ");
    static_assert(length == 2 || length == 4 || length == 6, "the size of the stencil is not supported");
    //-------------------------------------------------------------------------
    const real_t scale_d[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();

    constexpr real_t one_half        = 1.0 / 2.0;
    constexpr real_t two_third       = 2.0 / 3.0;
    constexpr real_t one_twelve      = 1.0 / 12.0;
    constexpr real_t one_sixtieth    = 1.0 / 60.0;
    constexpr real_t three_twentieth = 3.0 / 20.0;
    constexpr real_t three_forth     = 3.0 / 4.0;

    auto op = [=, &data_trg, *this](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writting
        const real_t* src = data_src + m_idx(i0, i1, i2);  // cache line for reading
        const real_t* u   = data_u + m_idx(i0, i1, i2);    // cache line for reading
        const real_t* v   = data_v + m_idx(i0, i1, i2);    // cache line for reading
        const real_t* w   = data_w + m_idx(i0, i1, i2);    // cache line for reading

        //advection
        if constexpr (length == 2) {
            trg[0] -= scale_d[0] * u[0] * (one_half * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
            trg[0] -= scale_d[1] * v[0] * (one_half * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
            trg[0] -= scale_d[2] * w[0] * (one_half * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
        } else if constexpr (length == 4) {
            trg[0] -= scale_d[0] * u[0] * (-one_twelve * (src[m_idx(+2, 0, 0)] - src[m_idx(-2, 0, 0)]) + two_third * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
            trg[0] -= scale_d[1] * v[0] * (-one_twelve * (src[m_idx(0, +2, 0)] - src[m_idx(0, -2, 0)]) + two_third * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
            trg[0] -= scale_d[2] * w[0] * (-one_twelve * (src[m_idx(0, 0, +2)] - src[m_idx(0, 0, -2)]) + two_third * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
        } else if constexpr (length == 6) {
            trg[0] -= scale_d[0] * u[0] * (+one_sixtieth * (src[m_idx(+3, 0, 0)] - src[m_idx(-3, 0, 0)]) - three_twentieth * (src[m_idx(+2, 0, 0)] - src[m_idx(-2, 0, 0)]) + three_forth * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
            trg[0] -= scale_d[1] * v[0] * (+one_sixtieth * (src[m_idx(0, +3, 0)] - src[m_idx(0, -3, 0)]) - three_twentieth * (src[m_idx(0, +2, 0)] - src[m_idx(0, -2, 0)]) + three_forth * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
            trg[0] -= scale_d[2] * w[0] * (+one_sixtieth * (src[m_idx(0, 0, +3)] - src[m_idx(0, 0, -3)]) - three_twentieth * (src[m_idx(0, 0, +2)] - src[m_idx(0, 0, -2)]) + three_forth * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
        }
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
    //-------------------------------------------------------------------------
}

/**
 * @brief compute the advection stencil on a given field using a conservative centered
 * 
 * @tparam length the size of the advection stencil. (e.g. 3 for 2nd order and 5 for 4th order)
 */
template <sid_t length>
class ConsAdvection : public Stencil {
   protected:
    m_ptr<const Field> u_;

   public:
    explicit ConsAdvection(m_ptr<const Field> u) : u_(u), Stencil(){};
    virtual lid_t NGhost() const override { return (length / 2); };

   protected:
    void DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) override;
};

template <sid_t length>
void ConsAdvection<length>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) {
    m_assert(fid_src->lda() == fid_trg->lda(), "the source and target field must have the same dimension ");
    m_assert(u_->ghost_status(), "the velocity must have up to date ghosts");
    static_assert(length == 4, "the size of the stencil is not supported");
    //-------------------------------------------------------------------------
    const real_t scale_d[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();

    constexpr real_t one_half        = 1.0 / 2.0;
    constexpr real_t one_sixth       = 1.0 / 6.0;
    constexpr real_t five_sixth      = 5.0 / 6.0;
    constexpr real_t one_third       = 1.0 / 3.0;
    constexpr real_t two_third       = 2.0 / 3.0;
    constexpr real_t one_twelve      = 1.0 / 12.0;
    constexpr real_t seven_twelve    = 7.0 / 12.0;
    constexpr real_t one_sixtieth    = 1.0 / 60.0;
    constexpr real_t three_twentieth = 3.0 / 20.0;
    constexpr real_t three_forth     = 3.0 / 4.0;

    auto op = [=, &data_trg, *this](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writing
        const real_t* src = data_src + m_idx(i0, i1, i2);  // cache line for reading
        const real_t* u   = data_u + m_idx(i0, i1, i2);    // cache line for reading
        const real_t* v   = data_v + m_idx(i0, i1, i2);    // cache line for reading
        const real_t* w   = data_w + m_idx(i0, i1, i2);    // cache line for reading

        //advection
        real_t flux[6];
        if constexpr (length == 4) {
            // compute the 6 fluxes

            const real_t u_x_m = 0.5 * (u[m_idx(0, 0, 0)] + u[m_idx(-1, 0, 0)]);
            const real_t u_x_p = 0.5 * (u[m_idx(0, 0, 0)] + u[m_idx(+1, 0, 0)]);
            const real_t v_y_m = 0.5 * (v[m_idx(0, 0, 0)] + v[m_idx(0, -1, 0)]);
            const real_t v_y_p = 0.5 * (v[m_idx(0, 0, 0)] + v[m_idx(0, +1, 0)]);
            const real_t w_z_m = 0.5 * (w[m_idx(0, 0, 0)] + w[m_idx(0, 0, -1)]);
            const real_t w_z_p = 0.5 * (w[m_idx(0, 0, 0)] + w[m_idx(0, 0, +1)]);

            flux[0] = (u_x_m > 0.0) * (-one_sixth * (src[m_idx(-2, 0, 0)] * u[m_idx(-2, 0, 0)]) + five_sixth * (src[m_idx(-1, 0, 0)] * u[m_idx(-1, 0, 0)]) + one_third * (src[m_idx(+0, 0, 0)] * u[m_idx(+0, 0, 0)])) +
                      (u_x_m < 0.0) * (+one_third * (src[m_idx(-1, 0, 0)] * u[m_idx(-1, 0, 0)]) + five_sixth * (src[m_idx(+0, 0, 0)] * u[m_idx(+0, 0, 0)]) - one_sixth * (src[m_idx(+1, 0, 0)] * u[m_idx(+1, 0, 0)]));
            flux[1] = (u_x_p > 0.0) * (-one_sixth * (src[m_idx(-1, 0, 0)] * u[m_idx(-1, 0, 0)]) + five_sixth * (src[m_idx(+0, 0, 0)] * u[m_idx(+0, 0, 0)]) + one_third * (src[m_idx(+1, 0, 0)] * u[m_idx(+1, 0, 0)])) +
                      (u_x_p < 0.0) * (+one_third * (src[m_idx(+0, 0, 0)] * u[m_idx(+0, 0, 0)]) + five_sixth * (src[m_idx(+1, 0, 0)] * u[m_idx(+1, 0, 0)]) - one_sixth * (src[m_idx(+2, 0, 0)] * u[m_idx(+2, 0, 0)]));
            flux[2] = (v_y_m > 0.0) * (-one_sixth * (src[m_idx(0, -2, 0)] * v[m_idx(0, -2, 0)]) + five_sixth * (src[m_idx(0, -1, 0)] * v[m_idx(0, -1, 0)]) + one_third * (src[m_idx(0, +0, 0)] * v[m_idx(0, +0, 0)])) +
                      (v_y_m < 0.0) * (+one_third * (src[m_idx(0, -1, 0)] * v[m_idx(0, -1, 0)]) + five_sixth * (src[m_idx(0, +0, 0)] * v[m_idx(0, +0, 0)]) - one_sixth * (src[m_idx(0, +1, 0)] * v[m_idx(0, +1, 0)]));
            flux[3] = (v_y_p > 0.0) * (-one_sixth * (src[m_idx(0, -1, 0)] * v[m_idx(0, -1, 0)]) + five_sixth * (src[m_idx(0, +0, 0)] * v[m_idx(0, +0, 0)]) + one_third * (src[m_idx(0, +1, 0)] * v[m_idx(0, +1, 0)])) +
                      (v_y_p < 0.0) * (+one_third * (src[m_idx(0, +0, 0)] * v[m_idx(0, +0, 0)]) + five_sixth * (src[m_idx(0, +1, 0)] * v[m_idx(0, +1, 0)]) - one_sixth * (src[m_idx(0, +2, 0)] * v[m_idx(0, +2, 0)]));
            flux[4] = (w_z_m > 0.0) * (-one_sixth * (src[m_idx(0, 0, -2)] * w[m_idx(0, 0, -2)]) + five_sixth * (src[m_idx(0, 0, -1)] * w[m_idx(0, 0, -1)]) + one_third * (src[m_idx(0, 0, +0)] * w[m_idx(0, 0, +0)])) +
                      (w_z_m < 0.0) * (+one_third * (src[m_idx(0, 0, -1)] * w[m_idx(0, 0, -1)]) + five_sixth * (src[m_idx(0, 0, +0)] * w[m_idx(0, 0, +0)]) - one_sixth * (src[m_idx(0, 0, +1)] * w[m_idx(0, 0, +1)]));
            flux[5] = (w_z_p > 0.0) * (-one_sixth * (src[m_idx(0, 0, -1)] * w[m_idx(0, 0, -1)]) + five_sixth * (src[m_idx(0, 0, +0)] * w[m_idx(0, 0, +0)]) + one_third * (src[m_idx(0, 0, +1)] * w[m_idx(0, 0, +1)])) +
                      (w_z_p < 0.0) * (+one_third * (src[m_idx(0, 0, +0)] * w[m_idx(0, 0, +0)]) + five_sixth * (src[m_idx(0, 0, +1)] * w[m_idx(0, 0, +1)]) - one_sixth * (src[m_idx(0, 0, +2)] * w[m_idx(0, 0, +2)]));

            // // this is the main contribution
            // trg[0] -= scale_d[0] * 0.5 * m_max(u[0] + u[m_idx(+1, 0, 0)],0.0) * (-one_twelve * (src[m_idx(+2, 0, 0)] + src[m_idx(-1, 0, 0)]) + seven_twelve * (src[m_idx(+1, 0, 0)] - src[m_idx(0, 0, 0)]));
            // trg[0] -= scale_d[1] * 0.5 * m_max(v[0] + v[m_idx(0, +1, 0)],0.0) * (-one_twelve * (src[m_idx(0, +2, 0)] + src[m_idx(0, -1, 0)]) + seven_twelve * (src[m_idx(0, +1, 0)] - src[m_idx(0, 0, 0)]));
            // trg[0] -= scale_d[2] * 0.5 * m_max(w[0] + w[m_idx(0, 0, +1)],0.0) * (-one_twelve * (src[m_idx(0, 0, +2)] + src[m_idx(0, 0, -1)]) + seven_twelve * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, 0)]));

            // // add the side contributions, the + side first
            // trg[0] -= scale_d[0] * 0.5 * u[m_idx(+1, 0, 0)] * (-one_twelve * (src[m_idx(+2, 0, 0)] + src[m_idx(-2, 0, 0)]) + two_third * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
            // trg[0] -= scale_d[1] * 0.5 * v[m_idx(0, +1, 0)] * (-one_twelve * (src[m_idx(0, +2, 0)] + src[m_idx(0, -2, 0)]) + two_third * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
            // trg[0] -= scale_d[2] * 0.5 * w[m_idx(0, 0, +1)] * (-one_twelve * (src[m_idx(0, 0, +2)] + src[m_idx(0, 0, -2)]) + two_third * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
        }

        // add the fluxes
        trg[0] -= scale_d[0] * (flux[1] - flux[0]);
        trg[0] -= scale_d[1] * (flux[3] - flux[2]);
        trg[0] -= scale_d[2] * (flux[5] - flux[4]);
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
    //-------------------------------------------------------------------------
}

/**
 * @brief implement the weno advection in the sense of Jiang 2000
 * 
 * The choice of the flux at the end is quite basic, maybe to be improved
 * 
 * @tparam length 
 */
template <sid_t length>
class WenoAdvection : public Stencil {
   protected:
    m_ptr<const Field> u_;

   public:
    explicit WenoAdvection(m_ptr<const Field> u) : u_(u), Stencil(){};
    virtual lid_t NGhost() const override { return (length / 2 + 1); };

   protected:
    void DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) override;
};

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

template <sid_t length>
void WenoAdvection<length>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) {
    m_assert(fid_src->lda() == fid_trg->lda(), "the source and target field must have the same dimension ");
    m_assert(u_->ghost_status(), "the velocity must have up to date ghosts");
    static_assert(length == 3 || length == 5, "only the weighted third order ENO has been coded");
    //-------------------------------------------------------------------------
    const real_t scale_d[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    const real_t* data_u   = block->data(u_, 0).Read();
    const real_t* data_v   = block->data(u_, 1).Read();
    const real_t* data_w   = block->data(u_, 2).Read();
    const real_t* data_src = block->data(fid_src, ida_).Read();
    real_t*       data_trg = block->data(fid_trg, ida_).Write();

    auto op = [=, &data_trg, *this](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_t*       trg = data_trg + m_idx(i0, i1, i2);  // cache line for writing
        const real_t* src = data_src + m_idx(i0, i1, i2);  // cache line for reading
        const real_t* u   = data_u + m_idx(i0, i1, i2);    // cache line for reading
        const real_t* v   = data_v + m_idx(i0, i1, i2);    // cache line for reading
        const real_t* w   = data_w + m_idx(i0, i1, i2);    // cache line for reading

        //advection
        real_t flux[6];
        if constexpr (length == 3) {  // third order
            constexpr real_t eps = 1.0e-6;
            {  // X
                const real_t* vel     = u;
                const real_t  h_inv   = scale_d[0];
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
                const real_t  h_inv   = scale_d[0];
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
                const real_t  h_inv   = scale_d[0];
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
        }
        if constexpr (length == 5) {  // fifth order
            constexpr real_t one_twelve = 1.0 / 12.0;
            {  // X
                const real_t* vel      = u;
                const real_t  h_inv    = scale_d[0];
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
                const real_t  h_inv    = scale_d[1];
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
                const real_t  h_inv    = scale_d[2];
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

            // const real_t u_x_m = 0.5 * (u[m_idx(0, 0, 0)] + u[m_idx(-1, 0, 0)]);
            // const real_t u_x_p = 0.5 * (u[m_idx(0, 0, 0)] + u[m_idx(+1, 0, 0)]);
            // const real_t v_y_m = 0.5 * (v[m_idx(0, 0, 0)] + v[m_idx(0, -1, 0)]);
            // const real_t v_y_p = 0.5 * (v[m_idx(0, 0, 0)] + v[m_idx(0, +1, 0)]);
            // const real_t w_z_m = 0.5 * (w[m_idx(0, 0, 0)] + w[m_idx(0, 0, -1)]);
            // const real_t w_z_p = 0.5 * (w[m_idx(0, 0, 0)] + w[m_idx(0, 0, +1)]);

            // // this is the main contribution
            // trg[0] -= scale_d[0] * 0.5 * m_max(u[0] + u[m_idx(+1, 0, 0)],0.0) * (-one_twelve * (src[m_idx(+2, 0, 0)] + src[m_idx(-1, 0, 0)]) + seven_twelve * (src[m_idx(+1, 0, 0)] - src[m_idx(0, 0, 0)]));
            // trg[0] -= scale_d[1] * 0.5 * m_max(v[0] + v[m_idx(0, +1, 0)],0.0) * (-one_twelve * (src[m_idx(0, +2, 0)] + src[m_idx(0, -1, 0)]) + seven_twelve * (src[m_idx(0, +1, 0)] - src[m_idx(0, 0, 0)]));
            // trg[0] -= scale_d[2] * 0.5 * m_max(w[0] + w[m_idx(0, 0, +1)],0.0) * (-one_twelve * (src[m_idx(0, 0, +2)] + src[m_idx(0, 0, -1)]) + seven_twelve * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, 0)]));

            // // add the side contributions, the + side first
            // trg[0] -= scale_d[0] * 0.5 * u[m_idx(+1, 0, 0)] * (-one_twelve * (src[m_idx(+2, 0, 0)] + src[m_idx(-2, 0, 0)]) + two_third * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]));
            // trg[0] -= scale_d[1] * 0.5 * v[m_idx(0, +1, 0)] * (-one_twelve * (src[m_idx(0, +2, 0)] + src[m_idx(0, -2, 0)]) + two_third * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]));
            // trg[0] -= scale_d[2] * 0.5 * w[m_idx(0, 0, +1)] * (-one_twelve * (src[m_idx(0, 0, +2)] + src[m_idx(0, 0, -2)]) + two_third * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]));
        }

        // add the fluxes -> godunov type, we take the max of the fluxes given the sign of the velocity
        trg[0] -= (u[0] > 0.0) ? flux[0] : flux[1];
        trg[0] -= (v[0] > 0.0) ? flux[2] : flux[3];
        trg[0] -= (w[0] > 0.0) ? flux[4] : flux[5];
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
    //-------------------------------------------------------------------------
}

#endif  // SRC_ADVECTION_DIFFUSION_HPP_