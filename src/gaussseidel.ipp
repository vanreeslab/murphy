#ifndef SRC_GAUSSSEIDEL_IPP_
#define SRC_GAUSSSEIDEL_IPP_

#include "gaussseidel.hpp"

template <sid_t length>
void GaussSeidel<length>::ApplyOpDerivInner(const qid_t* qid, GridBlock* block, const Field* fid_rhs, Field* fid_sol) {
    //-------------------------------------------------------------------------
    // we incorporate the alpha parameter to the scale coefficient
    real_t coef_tmp[length];
    real_t denom = -this->coef_[length / 2] / pow(block->hgrid(0), 2) - this->coef_[length / 2] / pow(block->hgrid(1), 2) - this->coef_[length / 2] / pow(block->hgrid(2), 2);
    for (sid_t is = -(length / 2); is <= (length / 2); is++) {
        coef_tmp[is] = (this->coef_[length / 2 + is]) * alpha_ / denom;
    }
    real_t coef_rhs = alpha_ / denom;
    // we reset the alpha coefficient so that the sum = 1 - alpha
    coef_tmp[length / 2] = (1.0 - alpha_) / 3.0;
    // set the coef array to the correct position
    real_t* coef = coef_tmp + (length / 2);

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
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    //-------------------------------------------------------------------------
}

template <sid_t length>
void GaussSeidel<length>::ApplyOpDerivOuter(const qid_t* qid, GridBlock* block, const Field* fid_rhs, Field* fid_sol) {
    //-------------------------------------------------------------------------
    // we incorporate the alpha parameter to the scale coefficient
    real_t coef_tmp[length];
    real_t denom = -this->coef_[length / 2] / pow(block->hgrid(0), 2) - this->coef_[length / 2] / pow(block->hgrid(1), 2) - this->coef_[length / 2] / pow(block->hgrid(2), 2);
    for (sid_t is = -(length / 2); is <= (length / 2); is++) {
        coef_tmp[is] = (this->coef_[length / 2 + is]) * alpha_ / denom;
    }
    real_t coef_rhs = alpha_ / denom;
    // we reset the alpha coefficient so that the sum = 1 - alpha
    coef_tmp[length / 2] = (1.0 - alpha_) / 3.0;
    // set the coef array to the correct position
    real_t* coef = coef_tmp + (length / 2);


    real_p data_rhs = block->data(fid_rhs, this->ida_);
    real_p data_sol = block->data(fid_sol, this->ida_);
    m_assume_aligned(data_sol);
    m_assume_aligned(data_rhs);

    //------ X -
    for (lid_t i2 = 0; i2 < M_N ; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_GS; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    //------ X +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = M_N - M_GS; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    //------ Y -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_GS; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    //------ Y +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = M_N - M_GS; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    //------ Z -
    for (lid_t i2 = 0; i2 < M_GS; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    //------ Z +
    for (lid_t i2 = M_N - M_GS; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the data pointer in front of the row for every cache line
                real_p       trg = data_sol + m_idx(i0, i1, i2);  // cache line for writting
                const real_p src = data_rhs + m_idx(i0, i1, i2);  // cache line for reading
                // loop on the stencil block
                real_t value = coef_rhs * src[0];
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    value += coef[is] * trg[m_idx(is, 0, 0)];
                    value += coef[is] * trg[m_idx(0, is, 0)];
                    value += coef[is] * trg[m_idx(0, 0, is)];
                }
                trg[0] = value;
            }
        }
    }
    // -------------------------------------------------------------------------
}


#endif  // SRC_GAUSSSEIDEL_IPP_