#include "blas.hpp"

Dcopy::Dcopy() : BlockOperator(nullptr){};
Dcopy::Dcopy(const Wavelet* interp) : BlockOperator(interp){};

void Dcopy::ComputeDcopyGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y) {
    m_assert(fid_x->lda() == fid_y->lda(), "the dimensions must match");
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        const data_ptr data_x = block->data(fid_x, ida);
        data_ptr       data_y = block->data(fid_y, ida);
        m_assume_aligned(data_x);
        m_assume_aligned(data_y);
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    const size_t idx = m_idx(i0, i1, i2);
                    data_y[idx]      = data_x[idx];
                }
            }
        }
    }
}

Daxpy::Daxpy(real_t alpha) : BlockOperator(nullptr) {
    //-------------------------------------------------------------------------
    alpha_ = alpha;
    //-------------------------------------------------------------------------
}

void Daxpy::ComputeDaxpyGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z) {
    m_assert(fid_x->lda() == fid_y->lda(), "the dimensions must match");
    m_assert(fid_y->lda() == fid_z->lda(), "the dimensions must match");
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        const data_ptr data_x = block->data(fid_x, ida);
        const data_ptr data_y = block->data(fid_y, ida);
        data_ptr       data_z = block->data(fid_z, ida);
        m_assume_aligned(data_x);
        m_assume_aligned(data_y);
        m_assume_aligned(data_z);
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = 0; i2 < M_N; i2++) {
            for (lid_t i1 = 0; i1 < M_N; i1++) {
                for (lid_t i0 = 0; i0 < M_N; i0++) {
                    const size_t idx = m_idx(i0, i1, i2);
                    data_z[idx]      = alpha_ * data_x[idx] + data_y[idx];
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}

Scale::Scale(real_t alpha) : BlockOperator(nullptr) {
    //-------------------------------------------------------------------------
    alpha_ = alpha;
    //-------------------------------------------------------------------------
}

void Scale::ComputeScaleGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x) {
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        const data_ptr data_x = block->data(fid_x, ida);
        m_assume_aligned(data_x);
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = 0; i2 < M_N; i2++) {
            for (lid_t i1 = 0; i1 < M_N; i1++) {
                for (lid_t i0 = 0; i0 < M_N; i0++) {
                    const size_t idx = m_idx(i0, i1, i2);
                    data_x[idx] *= alpha_;
                }
            }
        }
    }
    //-------------------------------------------------------------------------
}
