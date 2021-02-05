#ifndef SRC_ADVECTION_DIFFUSION_HPP_
#define SRC_ADVECTION_DIFFUSION_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "forloop.hpp"
#include "stencil.hpp"

/**
 * @brief compute the advection + diffusion stencil on a given field
 * 
 * @tparam length_advection the size of the advection stencil. (e.g. 3 for 2nd order and 5 for 4th order)
 * @tparam length_diffusion the size of the diffusion stencil. (e.g. 3 for 2nd order and 5 for 4th order)
 */
template <sid_t length_advection, sid_t length_diffusion>
class AdvectionDiffusion : public Stencil {
   protected:
    real_t nu_          = 0.0;
    real_t u_stream_[3] = {0.0, 0.0, 1.0};

   public:
    explicit AdvectionDiffusion(const real_t nu, const real_t u_stream[3]) : Stencil() {
        //-------------------------------------------------------------------------
        nu_ = nu;
        for (lda_t i = 0; i < 3; ++i) {
            u_stream_[i] = u_stream[i];
        }
        //-------------------------------------------------------------------------
    }

    virtual lid_t NGhost() const override { return m_max(length_advection / 2, length_diffusion / 2); };

   protected:
    void DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const override;
};

template <sid_t length_advection, sid_t length_diffusion>
void AdvectionDiffusion<length_advection, length_diffusion>::DoMagic(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, const bool is_outer, m_ptr<const Field> fid_src, m_ptr<Field> fid_trg) const {
    m_assert((fid_src->lda() == 3) && (fid_trg->lda() == 3), "the source must be 3 times smaller than the target ");
    static_assert(length_advection == 3 || length_advection == 5, "the size of the stencil is not supported");
    static_assert(length_diffusion == 3 || length_diffusion == 5, "the size of the stencil is not supported");
    //-------------------------------------------------------------------------
    const real_t scale_d[3]  = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};
    const real_t scale_d2[3] = {nu_ / pow(block->hgrid(0), 2), nu_ / pow(block->hgrid(1), 2), nu_ / pow(block->hgrid(2), 2)};

    real_p data_src = block->data(fid_src, ida_);
    real_p data_trg = block->data(fid_trg, ida_);
    m_assume_aligned(data_trg);
    m_assume_aligned(data_src);

    constexpr real_t two_third  = 2.0 / 3.0;
    constexpr real_t one_twelve = 1.0 / 12.0;
    constexpr real_t four_third = 4.0 / 3.0;
    constexpr real_t five_half  = 5.0 / 2.0;

    auto op = [=, &data_trg, *this](const lid_t i0, const lid_t i1, const lid_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_p       trg  = data_trg + m_idx(i0, i1, i2);                // cache line for writting
        const real_p src  = data_src + m_idx(i0, i1, i2);                // cache line for reading
        const real_t u[3] = {u_stream_[0], u_stream_[1], u_stream_[2]};  // todo add the real velocity here

        //advection
        if constexpr (length_advection == 3) {
            (*trg) -= u[0] * scale_d[0] * 0.5 * (src[m_idx(+1, 0, 0)] - src[m_idx(-1, 0, 0)]);
            (*trg) -= u[1] * scale_d[1] * 0.5 * (src[m_idx(0, +1, 0)] - src[m_idx(0, -1, 0)]);
            (*trg) -= u[2] * scale_d[2] * 0.5 * (src[m_idx(0, 0, +1)] - src[m_idx(0, 0, -1)]);
        } else if constexpr (length_advection == 5) {
            (*trg) -= u[0] * scale_d[0] * (-one_twelve * src[m_idx(+2, 0, 0)] + two_third * src[m_idx(+1, 0, 0)] - two_third * src[m_idx(-1, 0, 0)] + one_twelve * src[m_idx(-2, 0, 0)]);
            (*trg) -= u[1] * scale_d[1] * (-one_twelve * src[m_idx(0, +2, 0)] + two_third * src[m_idx(0, +1, 0)] - two_third * src[m_idx(0, -1, 0)] + one_twelve * src[m_idx(0, -2, 0)]);
            (*trg) -= u[2] * scale_d[2] * (-one_twelve * src[m_idx(0, 0, +2)] + two_third * src[m_idx(0, 0, +1)] - two_third * src[m_idx(0, 0, -1)] + one_twelve * src[m_idx(0, 0, -2)]);
        }
        // diffusion
        if constexpr (length_diffusion == 3) {
            (*trg) += scale_d2[0] * (src[m_idx(+1, 0, 0)] - 2.0 * src[m_idx(0, 0, 0)] + src[m_idx(-1, 0, 0)]);
            (*trg) += scale_d2[1] * (src[m_idx(0, +1, 0)] - 2.0 * src[m_idx(0, 0, 0)] + src[m_idx(0, -1, 0)]);
            (*trg) += scale_d2[2] * (src[m_idx(0, 0, +1)] - 2.0 * src[m_idx(0, 0, 0)] + src[m_idx(0, 0, -1)]);
        } else if constexpr (length_diffusion == 5) {
            (*trg) += scale_d2[0] * (-one_twelve * src[m_idx(+2, 0, 0)] + four_third * src[m_idx(+1, 0, 0)] - five_half * src[m_idx(0, 0, 0)] + four_third * src[m_idx(-1, 0, 0)] - one_twelve * src[m_idx(-2, 0, 0)]);
            (*trg) += scale_d2[1] * (-one_twelve * src[m_idx(0, +2, 0)] + four_third * src[m_idx(0, +1, 0)] - five_half * src[m_idx(0, 0, 0)] + four_third * src[m_idx(0, -1, 0)] - one_twelve * src[m_idx(0, -2, 0)]);
            (*trg) += scale_d2[2] * (-one_twelve * src[m_idx(0, 0, +2)] + four_third * src[m_idx(0, 0, +1)] - five_half * src[m_idx(0, 0, 0)] + four_third * src[m_idx(0, 0, -1)] - one_twelve * src[m_idx(0, 0, -2)]);
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

#endif  // SRC_ADVECTION_DIFFUSION_HPP_