#include "blas.hpp"

//-----------------------------------------------------------------------------
Dset::Dset() : BlockOperator(nullptr){};
Dset::Dset(m_ptr<const Wavelet> interp) : BlockOperator(interp){};

void Dset::operator()(m_ptr<const ForestGrid> grid, const real_t value, m_ptr<Field> fid_x) {
    m_begin;
    //-------------------------------------------------------------------------
    value_ = value;

    DoOpMesh(this, &Dset::ComputeDsetGridBlock, grid, fid_x);
    // update the ghost
    fid_x->ghost_status(this->do_ghost());
    //-------------------------------------------------------------------------
    m_end;
}

void Dset::ComputeDsetGridBlock(m_ptr<const qid_t> qid, m_ptr<const GridBlock> block, m_ptr<Field> fid_x) {
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        real_t* data_x = block->data(fid_x, ida).Write();
        // m_assume_aligned(data_x);
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    const size_t idx = m_idx(i0, i1, i2);
                    data_x[idx]      = value_;
                }
            }
        }
    }
}

//-----------------------------------------------------------------------------
Dcopy::Dcopy() : BlockOperator(nullptr){};
Dcopy::Dcopy(m_ptr<const Wavelet> interp) : BlockOperator(interp){};

void Dcopy::operator()(m_ptr<const ForestGrid> grid, m_ptr<const Field> fid_x, m_ptr<Field> fid_y) {
    m_begin;
    //-------------------------------------------------------------------------
    DoOpMesh(this, &Dcopy::ComputeDcopyGridBlock, grid, fid_x, fid_y);
    // update the ghost
    fid_y->ghost_status(this->do_ghost());
    //-------------------------------------------------------------------------
    m_end;
}

void Dcopy::ComputeDcopyGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid_x, m_ptr<Field> fid_y) {
    m_assert(fid_x->lda() == fid_y->lda(), "the dimensions must match");
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        const real_t* data_x = block->data(fid_x, ida).Read();
        real_t*       data_y = block->data(fid_y, ida).Write();

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

//-----------------------------------------------------------------------------
Daxpy::Daxpy() : BlockOperator(nullptr){};
Daxpy::Daxpy(m_ptr<const Wavelet> interp) : BlockOperator(interp){};

void Daxpy::operator()(m_ptr<const ForestGrid> grid, const real_t alpha, m_ptr<const Field> fid_x, m_ptr<const Field> fid_y, m_ptr<Field> fid_z) {
    m_begin;
    //-------------------------------------------------------------------------
    alpha_ = alpha;
    DoOpMesh(this, &Daxpy::ComputeDaxpyGridBlock, grid, fid_x, fid_y, fid_z);

    // update the ghost
    fid_z->ghost_status(this->do_ghost());
    //-------------------------------------------------------------------------
    m_end;
}

void Daxpy::ComputeDaxpyGridBlock(m_ptr<const qid_t> qid, m_ptr<const GridBlock> block, m_ptr<const Field> fid_x, m_ptr<const Field> fid_y, m_ptr<Field> fid_z) {
    m_assert(fid_x->lda() == fid_y->lda(), "the dimensions must match");
    m_assert(fid_y->lda() == fid_z->lda(), "the dimensions must match");
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        const real_t* data_x = block->data(fid_x, ida).Read();
        const real_t* data_y = block->data(fid_y, ida).Read();
        real_t*       data_z = block->data(fid_z, ida).Write();
        // m_assume_aligned(data_x);
        // m_assume_aligned(data_y);
        // m_assume_aligned(data_z);
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

//-----------------------------------------------------------------------------
Dscale::Dscale() : BlockOperator(nullptr){};
Dscale::Dscale(m_ptr<const Wavelet> interp) : BlockOperator(interp){};

void Dscale::operator()(const ForestGrid* grid, const real_t alpha, m_ptr<Field> fid_x) {
    m_begin;
    //-------------------------------------------------------------------------
    alpha_ = alpha;
    DoOpMesh(this, &Dscale::ComputeDscaleGridBlock, grid, fid_x);
    // update the ghost
    fid_x->ghost_status(this->do_ghost());
    //-------------------------------------------------------------------------
    m_end;
}

void Dscale::ComputeDscaleGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<Field> fid_x) {
    //-------------------------------------------------------------------------
    const sid_t lda = fid_x->lda();
    for (sid_t ida = 0; ida < lda; ida++) {
        // get the data pointers
        real_t* data_x = block->data(fid_x, ida).Write();
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
