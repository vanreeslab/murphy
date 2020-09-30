#include "daxpy.hpp"

Daxpy::Daxpy(real_t alpha) : BlockOperator(nullptr) {
    alpha_ = alpha;
}

void Daxpy::ComputeDaxpyGridBlock(const qid_t* qid, GridBlock* block, Field* fid_x, Field* fid_y, Field* fid_z) {
    m_assert(fid_x->lda() == fid_y->lda(), "the dimensions must match");
    m_assert(fid_y->lda() == fid_z->lda(), "the dimensions must match");
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        real_p data_x = block->data(fid_x, ida);
        real_p data_y = block->data(fid_y, ida);
        real_p data_z = block->data(fid_z, ida);
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
