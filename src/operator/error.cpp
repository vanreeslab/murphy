#include "error.hpp"

// declare the specialization, implement them in the cpp
template <>
void Error::ErrorOnGridBlock<Field>(const qid_t* qid, const CartBlock* block, const Field* fid, const Field* sol) {
    //-------------------------------------------------------------------------
    real_t e2 = 0.0;
    real_t ei = 0.0;

    for (lda_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const ConstMemData data_field = block->data(fid, ida);
        const ConstMemData data_sol   = block->data(sol, ida);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // we need to discard the physical BC for the edges
            const real_t error = data_field(i0, i1, i2) - data_sol(i0, i1, i2);
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2));
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, span_);
    }

    //#pragma omp critical
    {  // no max atomic in OpenMP
        const real_t* hgrid = block->hgrid();
        error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
};
template <>
void Error::ErrorFieldOnGridBlock<Field>(const qid_t* qid, const CartBlock* block, const Field* fid, const Field* sol, const Field* error) {
    //-------------------------------------------------------------------------
    real_t e2 = 0.0;
    real_t ei = 0.0;

    for (lda_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const ConstMemData data_field = block->data(fid, ida);
        const ConstMemData data_sol   = block->data(sol, ida);
        const MemData      data_error = block->data(error, ida);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the local error
            real_t error = data_field(i0, i1, i2) - data_sol(i0, i1, i2);
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2));
            // store it
            data_error(i0, i1, i2) = error;
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, span_);
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
void Error::ErrorOnGridBlock<lambda_error_t>(const qid_t* qid, const CartBlock* block, const Field* fid, const lambda_error_t* sol) {
    //-------------------------------------------------------------------------
    real_t e2 = 0.0;
    real_t ei = 0.0;

    for (lda_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const ConstMemData data_field = block->data(fid, ida);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // we need to discard the physical BC for the edges
            const real_t sol_val  = (*sol)(i0, i1, i2, block);
            const real_t data_val = data_field(i0, i1, i2);
            const real_t error    = data_val - sol_val;
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f vs %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2), sol_val);
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, span_);
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
void Error::ErrorFieldOnGridBlock<lambda_error_t>(const qid_t* qid, const CartBlock* block, const Field* fid, const lambda_error_t* sol, const Field* error) {
    //-------------------------------------------------------------------------
    real_t e2 = 0.0;
    real_t ei = 0.0;

    for (lda_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers
        const ConstMemData data_field = block->data(fid, ida);
        const MemData      data_error = block->data(error, ida);

        auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
            // get the local error
            real_t error = data_field(i0, i1, i2) - (*sol)(i0, i1, i2, block);
            m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2));
            // store it
            data_error(i0, i1, i2) = error;
            // update the block errors
            e2 += error * error;
            ei = m_max(std::fabs(error), ei);
        };
        for_loop(&op, span_);
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
