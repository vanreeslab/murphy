#ifndef SRC_GAUSSSEIDEL_IPP_
#define SRC_GAUSSSEIDEL_IPP_

#include "gaussseidel.hpp"

template <sid_t length>
void GaussSeidel<length>::ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_sol, Field* fid_rhs) {
    //-------------------------------------------------------------------------
    // take the stencil coef
    real_t* coef     = this->coef_ + (length / 2);
    real_t  scale[3] = {1.0 / pow(block->hgrid(0), 2), 1.0 / pow(block->hgrid(1), 2), 1.0 / pow(block->hgrid(2), 2)};
    // compute GS coefficients
    real_t coef_gs[length];
    real_t denom    = -coef[0] * (scale[0] + scale[1] + scale[2]);
    real_t coef_rhs = (alpha_ / denom);
    for (sid_t is = -(length / 2); is <= (length / 2); is++) {
        coef_gs[(length / 2) + is] = coef[is] * coef_rhs;
    }
    coef    = coef_gs + (length / 2);
    coef[0] = 0.0;

    // set the coef array to the correct position
    real_p data_rhs = block->data(fid_rhs, this->ida_);
    real_p data_sol = block->data(fid_sol, this->ida_);
    m_assume_aligned(data_rhs);
    m_assume_aligned(data_sol);

    // we ignore the dependencies, it will add some shuffle, which is fine
    // #pragma ivdep
    for (lid_t i2 = M_GS; i2 < M_N - M_GS; i2++) {
        for (lid_t i1 = M_GS; i1 < M_N - M_GS; i1++) {
            for (lid_t i0 = M_GS; i0 < M_N - M_GS; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    //-------------------------------------------------------------------------
}

template <sid_t length>
void GaussSeidel<length>::ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_sol, Field* fid_rhs) {
    //-------------------------------------------------------------------------
    // take the stencil coef
    real_t* coef     = this->coef_ + (length / 2);
    real_t  scale[3] = {1.0 / pow(block->hgrid(0), 2), 1.0 / pow(block->hgrid(1), 2), 1.0 / pow(block->hgrid(2), 2)};
    // compute GS coefficients
    real_t coef_gs[length];
    real_t denom    = -coef[0] * (scale[0] + scale[1] + scale[2]);
    real_t coef_rhs = (alpha_ / denom);
    for (sid_t is = -(length / 2); is <= (length / 2); is++) {
        coef_gs[(length / 2) + is] = coef[is] * coef_rhs;
    }
    coef    = coef_gs + (length / 2);
    coef[0] = 0.0;

    real_p data_rhs = block->data(fid_rhs, this->ida_);
    real_p data_sol = block->data(fid_sol, this->ida_);
    m_assume_aligned(data_sol);
    m_assume_aligned(data_rhs);

    //------ X -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_GS; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);

                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    //------ X +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = M_N - M_GS; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);

                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    //------ Y -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_GS; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);

                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    //------ Y +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = M_N - M_GS; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);

                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    //------ Z -
    for (lid_t i2 = 0; i2 < M_GS; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);

                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    //------ Z +
    for (lid_t i2 = M_N - M_GS; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);
                const real_p src = data_rhs + m_idx(i0, i1, i2);

                real_t value = -coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * scale[0] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * scale[1] * trg[m_idx(0, is, 0)];
                    value += coef[is] * scale[2] * trg[m_idx(0, 0, is)];
                }
                trg[0] = (1.0 - alpha_) * trg[0] + value;
            }
        }
    }
    // -------------------------------------------------------------------------
}

#endif  // SRC_GAUSSSEIDEL_IPP_