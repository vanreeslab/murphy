#ifndef SRC_JACOBI_INC_
#define SRC_JACOBI_INC_

#include <cstring>

#include "defs.hpp"
#include "jacobi.hpp"

using std::memcpy;

template <sid_t length>
void Jacobi<length>::IterativeSolverPrep(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp) {
    m_assert(fid_sol != nullptr, "solution cannot be null");
    m_assert(fid_tmp != nullptr, "temporay cannot be null");
    //-------------------------------------------------------------------------
    // set the coef array to the correct position
    real_p data_sol = block->data(fid_sol, this->ida_);
    real_p data_tmp = block->data(fid_tmp, this->ida_);
    m_assume_aligned(data_sol);
    m_assume_aligned(data_tmp);
    // copy the current solution to the temp
    memcpy(data_tmp - m_zeroidx(0, block), data_sol - m_zeroidx(0, block), sizeof(real_t) * m_blockmemsize(1));
    //-------------------------------------------------------------------------
};

template <sid_t length>
void Jacobi<length>::IterativeSolverInner(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp) {
    m_assert(fid_sol != nullptr, "solution cannot be null");
    m_assert(fid_rhs != nullptr, "right hande side cannot be null");
    m_assert(fid_tmp != nullptr, "temporay cannot be null");
    //-------------------------------------------------------------------------
    // get the laplacian stencil and scale factor
    real_t* coef_lapla = coef_lapla_ + (length / 2);
    real_t  scale[3]   = {1.0 / pow(block->hgrid(0), 2), 1.0 / pow(block->hgrid(1), 2), 1.0 / pow(block->hgrid(2), 2)};
    real_t  denom      = coef_lapla[0] * (scale[0] + scale[1] + scale[2]);
    real_t  coef_mult  = alpha_ / denom;
    // set the coef array to the correct position
    real_p data_rhs = block->data(fid_rhs, this->ida_);
    real_p data_sol = block->data(fid_sol, this->ida_);
    real_p data_tmp = block->data(fid_tmp, this->ida_);
    m_assume_aligned(data_rhs);
    m_assume_aligned(data_sol);
    m_assume_aligned(data_tmp);

    // we ignore the dependencies, it will add some shuffle, which is fine
    for (lid_t i2 = M_GS; i2 < M_N - M_GS; i2++) {
        for (lid_t i1 = M_GS; i1 < M_N - M_GS; i1++) {
            for (lid_t i0 = M_GS; i0 < M_N - M_GS; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //-------------------------------------------------------------------------
};

template <sid_t length>
void Jacobi<length>::IterativeSolverOuter(const qid_t* qid, GridBlock* block, Field* fid_sol, Field* fid_rhs, Field* fid_tmp) {
    //-------------------------------------------------------------------------
    // get the laplacian stencil and scale factor
    real_t* coef_lapla = coef_lapla_ + (length / 2);
    real_t  scale[3]   = {1.0 / pow(block->hgrid(0), 2), 1.0 / pow(block->hgrid(1), 2), 1.0 / pow(block->hgrid(2), 2)};
    real_t  denom      = coef_lapla[0] * (scale[0] + scale[1] + scale[2]);
    real_t  coef_mult  = alpha_ / denom;

    // set the coef array to the correct position
    real_p data_rhs = block->data(fid_rhs, this->ida_);
    real_p data_sol = block->data(fid_sol, this->ida_);
    real_p data_tmp = block->data(fid_tmp, this->ida_);
    m_assume_aligned(data_rhs);
    m_assume_aligned(data_sol);
    m_assume_aligned(data_tmp);

    //--------- X -
    for (lid_t i2 = 0; i2 < M_GS; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //--------- X +
    for (lid_t i2 = M_N - M_GS; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //--------- Y -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_GS; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //--------- Y +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = M_N - M_GS; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_N; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //--------- Z -
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = 0; i0 < M_GS; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //--------- Z +
    for (lid_t i2 = 0; i2 < M_N; i2++) {
        for (lid_t i1 = 0; i1 < M_N; i1++) {
            for (lid_t i0 = M_N - M_GS; i0 < M_N; i0++) {
                // get the local data pointer
                real_p       lsol = data_sol + m_idx(i0, i1, i2);
                const real_p lrhs = data_rhs + m_idx(i0, i1, i2);
                const real_p ltmp = data_tmp + m_idx(i0, i1, i2);
                // loop on the stencil block
                real_t res = lrhs[0] ;
                for (sid_t is = -(length / 2); is <= (length / 2); is++) {
                    res -= coef_lapla[is] * scale[0] * ltmp[m_idx(is, 0, 0)];
                    res -= coef_lapla[is] * scale[1] * ltmp[m_idx(0, is, 0)];
                    res -= coef_lapla[is] * scale[2] * ltmp[m_idx(0, 0, is)];
                }
                lsol[0] = ltmp[0] + coef_mult * res;
            }
        }
    }
    //-------------------------------------------------------------------------
};

#endif  // SRC_JACOBI_INC_