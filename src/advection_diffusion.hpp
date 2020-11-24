#ifndef SRC_ADVECTION_DIFFUSION_HPP_
#define SRC_ADVECTION_DIFFUSION_HPP_

#include "forloop.hpp"
#include "defs.hpp"
#include "stencil.hpp"

template <sid_t length>
class AdvectionDiffusion : public Stencil {
   protected:
    real_t nu_          = 0.0;
    real_t u_stream_[3] = {0.0, 0.0, 1.0};

    real_t coef_d_[length];   //!< coefficients for the first derivative
    real_t coef_d2_[length];  //!< coefficients for the second derivative

   public:
    explicit AdvectionDiffusion(const real_t nu, const real_t u_stream[3]) : Stencil() {
        //-------------------------------------------------------------------------
        nu_ = nu;
        for (lda_t i = 0; i < 3; ++i) {
            u_stream_[i] = u_stream[i];
        }

        // get the coefficients
        if constexpr (length == 3) {
            coef_d_[0]  = -1.0 / 2.0;
            coef_d_[1]  = 0.0;
            coef_d_[2]  = +1.0 / 2.0;
            coef_d2_[0] = +1.0;
            coef_d2_[1] = -2.0;
            coef_d2_[2] = +1.0;
        } else {
            m_assert(false, "not coded yet");
        }
        //-------------------------------------------------------------------------
    }

    virtual void ApplyStencilInner(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) override;
    virtual void ApplyStencilOuter(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) override;

   protected:
    void DoMagic(const bool is_outer, const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg);
};

template <sid_t length>
void AdvectionDiffusion<length>::ApplyStencilInner(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) {
    DoMagic(false, qid, block, fid_src, fid_trg);
};

template <sid_t length>
void AdvectionDiffusion<length>::ApplyStencilOuter(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) {
    DoMagic(true, qid, block, fid_src, fid_trg);
};

template <sid_t length>
void AdvectionDiffusion<length>::DoMagic(const bool is_outer, const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) {
    m_assert((fid_src->lda() == 3) && (fid_trg->lda() == 3), "the source must be 3 times smaller than the target ");
    //-------------------------------------------------------------------------
    const real_t* coef_d      = coef_d_ + (length / 2);
    const real_t* coef_d2     = coef_d2_ + (length / 2);
    const real_t  scale_d[3]  = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};
    const real_t  scale_d2[3] = {nu_ / pow(block->hgrid(0), 2), nu_ / pow(block->hgrid(1), 2), nu_ / pow(block->hgrid(2), 2)};

    real_p data_src = block->data(fid_src, ida_);
    real_p data_trg = block->data(fid_trg, ida_);
    m_assume_aligned(data_trg);
    m_assume_aligned(data_src);

    // default = capture by value = read-only
    auto op = [=, &data_trg, *this](const lid_t i0, const lid_t i1, const lid_t i2) -> void {
        // get the data pointer in front of the row for every cache line
        real_p       trg  = data_trg + m_idx(i0, i1, i2);                // cache line for writting
        const real_p src  = data_src + m_idx(i0, i1, i2);                // cache line for reading
        const real_t u[3] = {u_stream_[0], u_stream_[1], u_stream_[2]};  // todo add the real velocity here

        // loop on the stencil block
        (*trg) = 0.0;
        for (int is = -(length / 2); is <= (length / 2); is++) {
            // advection: ux * d/dx + uy * d/dy + ux * d/dz
            (*trg) -= u[0] * coef_d[is] * scale_d[0] * src[m_idx(is, 0, 0)];
            (*trg) -= u[1] * coef_d[is] * scale_d[1] * src[m_idx(0, is, 0)];
            (*trg) -= u[2] * coef_d[is] * scale_d[2] * src[m_idx(0, 0, is)];

            // diffusion nu * d/dx + nu* * d/dy + nu * d/dz (nu is in the scale_d2 factors)
            (*trg) += coef_d2[is] * scale_d2[0] * src[m_idx(is, 0, 0)];
            (*trg) += coef_d2[is] * scale_d2[1] * src[m_idx(0, is, 0)];
            (*trg) += coef_d2[is] * scale_d2[2] * src[m_idx(0, 0, is)];
        }
    };

    if (!is_outer) {
        for_loop<M_GS, M_N - M_GS>(op);
    } else {
        for_loop<0, M_GS, 0, M_N, 0, M_N>(op);          // X-
        for_loop<0, M_N, 0, M_GS, 0, M_N>(op);          // Y-
        for_loop<0, M_N, 0, M_N, 0, M_GS>(op);          // Z-
        for_loop<M_N - M_GS, M_N, 0, M_N, 0, M_N>(op);  // X+
        for_loop<0, M_N, M_N - M_GS, M_N, 0, M_N>(op);  // Y+
        for_loop<0, M_N, 0, M_N, M_N - M_GS, M_N>(op);  // Z+
    }
    //-------------------------------------------------------------------------
}

#endif  // SRC_ADVECTION_DIFFUSION_HPP_