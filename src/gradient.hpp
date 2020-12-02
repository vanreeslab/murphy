#ifndef SRC_GRADIENT_HPP_
#define SRC_GRADIENT_HPP_

#include "defs.hpp"
#include "stencil.hpp"

/**
 * @brief Implements a stencil laplacian in cross of order (length+1)/2
 * 
 * e.g. if the length is 3, the order is 2
 * 
 * @tparam length 
 */
template <sid_t length>
class Gradient : public Stencil {
   protected:
    real_t coef_[length];  //!< coefficients of the laplacian to apply

   public:
    explicit Gradient() : Stencil() {
        // get the stencil
        if constexpr (length == 3) {
            coef_[0] = -1.0 / 2.0;
            coef_[1] = 0.0;
            coef_[2] = +1.0 / 2.0;
        } else if (length == 5) {
            coef_[0] = +1.0 / 12.0;
            coef_[1] = -2.0 / 3.0;
            coef_[2] = 0.0;
            coef_[3] = +2.0 / 3.0;
            coef_[4] = -1.0 / 12.0;
        } else {
            m_assert(false, "not coded yet");
        }
    }

    virtual lid_t NGhost() const override { return length / 2; }

    virtual void ApplyStencilInner(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) override;
    virtual void ApplyStencilOuter(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) override;
};

/**
 * @brief computes the inner lapgradientlacian contribution of the dimension Stencil::ida_ of the source field on the traget field
 * 
 * @tparam length
 * @param qid
 * @param block 
 * @param fid_src the source field
 * @param fid_trg the target field
 */
template <sid_t length>
void Gradient<length>::ApplyStencilInner(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) {
    m_assert((fid_src->lda() == 1) && (fid_trg->lda() == 3), "the source must be 3 times smaller than the target ");
    //-------------------------------------------------------------------------
    const real_t* coef     = coef_ + (length / 2);
    const real_t  scale[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    real_p data_src   = block->data(fid_src, ida_);
    real_p data_trg_0 = block->data(fid_trg, 0);
    real_p data_trg_1 = block->data(fid_trg, 1);
    real_p data_trg_2 = block->data(fid_trg, 2);
    m_assume_aligned(data_trg_0);
    m_assume_aligned(data_trg_1);
    m_assume_aligned(data_trg_2);
    m_assume_aligned(data_src);

    for (lid_t i2 = M_GS; i2 < M_N - M_GS; ++i2) {
        for (lid_t i1 = M_GS; i1 < M_N - M_GS; ++i1) {
            for (lid_t i0 = M_GS; i0 < M_N - M_GS; ++i0) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

/**
 * @brief computes the 6 different outer laplacian contribution of the dimension Stencil::ida_ of the source field's ghost values on the traget field
 * 
 * @tparam length 
 * @param qid 
 * @param block 
 * @param fid_src 
 * @param fid_trg 
 */
template <sid_t length>
void Gradient<length>::ApplyStencilOuter(const qid_t* qid, GridBlock* block, Field* fid_src, Field* fid_trg) {
    m_assert(fid_src->lda() == (fid_trg->lda() / 3), "the source must be 3 times smaller than the target ");
    //-------------------------------------------------------------------------
    const real_t* coef     = coef_ + (length / 2);
    const real_t  scale[3] = {1.0 / block->hgrid(0), 1.0 / block->hgrid(1), 1.0 / block->hgrid(2)};

    real_p data_src   = block->data(fid_src, ida_);
    real_p data_trg_0 = block->data(fid_trg, 0);
    real_p data_trg_1 = block->data(fid_trg, 1);
    real_p data_trg_2 = block->data(fid_trg, 2);
    m_assume_aligned(data_trg_0);
    m_assume_aligned(data_trg_1);
    m_assume_aligned(data_trg_2);
    m_assume_aligned(data_src);

    //------ X -
    for (lid_t i2 = 0; i2 < M_N ; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_GS; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ X +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = M_N - M_GS; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Y -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_GS; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Y +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = M_N - M_GS; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Z -
    for (lid_t i2 = 0; i2 < M_GS; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    //------ Z +
    for (lid_t i2 = M_N - M_GS; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg_0 = data_trg_0 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_1 = data_trg_1 + m_idx(i0, i1, i2);  // cache line for writting
                real_p       trg_2 = data_trg_2 + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src   = data_src + m_idx(i0, i1, i2);    // cache line for reading
                // loop on the stencil block
                (*trg_0) = 0.0;
                (*trg_1) = 0.0;
                (*trg_2) = 0.0;
                for (int is = -(length / 2); is <= (length / 2); is++) {
                    (*trg_0) += coef[is] * scale[0] * src[m_idx(is, 0, 0)];
                    (*trg_1) += coef[is] * scale[1] * src[m_idx(0, is, 0)];
                    (*trg_2) += coef[is] * scale[2] * src[m_idx(0, 0, is)];
                }
            }
        }
    }
    // -------------------------------------------------------------------------
}


#endif  // SRC_GRADIENT_HPP_
