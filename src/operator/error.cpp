#include "error.hpp"

// declare the specialization, implement them in the cpp
template <>
void Error::ErrorOnGridBlock<Field>(const qid_t* qid, GridBlock* block, const Field* fid, const Field* sol) {
    //-------------------------------------------------------------------------
    const_data_ptr ptr_field = block->data(fid);
    const_data_ptr ptr_sol   = block->data(sol);

    real_t e2 = 0.0;
    real_t ei = 0.0;
    for (sid_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const real_t* data_field = ptr_field.Read(ida, block);
        const real_t* data_sol   = ptr_sol.Read(ida, block);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // we need to discard the physical BC for the edges
            const real_t error = data_field[m_idx(i0, i1, i2)] - data_sol[m_idx(i0, i1, i2)];
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field[m_idx(i0, i1, i2)]);
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, start_, end_);
    }

    // add the result
    const real_t* hgrid = block->hgrid();

    //#pragma omp critical
    {  // no max atomic in OpenMP
        error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
};
template <>
void Error::ErrorFieldOnGridBlock<Field>(const qid_t* qid, GridBlock* block, const Field* fid, const Field* sol, const Field* error) {
    //-------------------------------------------------------------------------
    const_data_ptr ptr_field = block->data(fid);
    const_data_ptr ptr_sol   = block->data(sol);
    data_ptr       ptr_error = block->data(error);

    real_t e2 = 0.0;
    real_t ei = 0.0;
    for (sid_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const real_t* data_field = ptr_field.Read(ida, block);
        const real_t* data_sol   = ptr_sol.Read(ida, block);
        real_t*       data_error = ptr_error.Write(ida, block);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the local error
            real_t error = data_field[m_idx(i0, i1, i2)] - data_sol[m_idx(i0, i1, i2)];
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field[m_idx(i0, i1, i2)]);
            // store it
            data_error[m_idx(i0, i1, i2)] = error;
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, start_, end_);
    }
    // add the result
    const real_t* hgrid = block->hgrid();

    //#pragma omp critical
    {  // no max atomic in OpenMP
        error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
};

template <>
void Error::ErrorOnGridBlock<lambda_error_t>(const qid_t* qid, GridBlock* block, const Field* fid, const lambda_error_t* sol) {
    //-------------------------------------------------------------------------
    const_data_ptr ptr_field = block->data(fid);

    real_t e2 = 0.0;
    real_t ei = 0.0;
    for (sid_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const real_t* data_field = ptr_field.Read(ida, block);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // we need to discard the physical BC for the edges
            const real_t sol_val  = (*sol)(i0, i1, i2, block);
            const real_t data_val = data_field[m_idx(i0, i1, i2)];
            const real_t error    = data_val - sol_val;
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field[m_idx(i0, i1, i2)]);
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, start_, end_);
    }

    // add the result
    const real_t* hgrid = block->hgrid();

    //#pragma omp critical
    {  // no max atomic in OpenMP
        error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
};

template <>
void Error::ErrorFieldOnGridBlock<lambda_error_t>(const qid_t* qid, GridBlock* block, const Field* fid, const lambda_error_t* sol, const Field* error) {
    //-------------------------------------------------------------------------
    const_data_ptr ptr_field = block->data(fid);
    data_ptr       ptr_error = block->data(error);

    real_t e2 = 0.0;
    real_t ei = 0.0;
    for (sid_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const real_t* data_field = ptr_field.Read(ida, block);
        real_t*       data_error = ptr_error.Write(ida, block);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the local error
            real_t error = data_field[m_idx(i0, i1, i2)] - (*sol)(i0, i1, i2, block);
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field[m_idx(i0, i1, i2)]);
            // store it
            data_error[m_idx(i0, i1, i2)] = error;
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, start_, end_);
    }
    // add the result
    const real_t* hgrid = block->hgrid();

    //#pragma omp critical
    {  // no max atomic in OpenMP
        error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
};
