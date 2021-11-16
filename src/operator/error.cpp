#include "error.hpp"

#include "core/forloop.hpp"

template <>
void Error::ErrorMagic<Field>(const ForestGrid* grid, const level_t level, const Field* field, const Field* sol, Field* error, real_t* norm_2, real_t* norm_i) const {
    m_assert(IsGhostValid(field), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
    m_assert(IsGhostValid(sol), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
    //--------------------------------------------------------------------------
    for (lda_t ida = 0; ida < field->lda(); ++ida) {
        const bool no_store = (error == nullptr);
        // define the operator
        real_t error_i = 0.0;
        auto   op      = [=, &error_i](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* block) -> real_t {
            // get the data pointers
            const ConstMemData data_field = block->data(field, ida);
            const ConstMemData data_sol   = block->data(sol, ida);

            // we need to discard the physical BC for the edges
            const real_t value = data_field(i0, i1, i2) - data_sol(i0, i1, i2);
            m_assert(std::isfinite(value), "the error cannot be nan: @ %d %d %d: %f", i0, i1, i2, data_field(i0, i1, i2));

            if (!no_store) {
                const MemData data_error = block->data(error, ida);
                data_error(i0, i1, i2)   = value;
            }

            // ei is camputed and updated
            const bool is_erri_valid = (span_.start[0] <= i0) && (i0 < span_.end[0]) &&
                                       (span_.start[1] <= i1) && (i1 < span_.end[1]) &&
                                       (span_.start[2] <= i2) && (i2 < span_.end[2]);
            error_i = m_max(std::fabs(value)*is_erri_valid, error_i);
            /// error 2 is returned to be integrated
            const real_t error2 = value * value;
            return error2;
        };
        // integrate error 2 and get the errori
        real_t error_2 = 0.0;
        if (level >= 0) {
            ComputeIntegral(grid, level, op, &error_2);
        } else {
            ComputeIntegral(grid, op, &error_2);
        }

        // finalize the computation
        if (norm_i != nullptr) {
            MPI_Allreduce(&error_i,norm_i+ida,1,M_MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
        }
        if (norm_2 != nullptr) {
            norm_2[ida] = sqrt(error_2);
        }
    }
    //--------------------------------------------------------------------------
}

template <>
void Error::ErrorMagic<lambda_error_t>(const ForestGrid* grid, const level_t level, const Field* field, const lambda_error_t* sol, Field* error, real_t* norm_2, real_t* norm_i) const{
    m_assert(IsGhostValid(field), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
    //--------------------------------------------------------------------------
    for (lda_t ida = 0; ida < field->lda(); ++ida) {
        // define the operator
        real_t error_i = 0.0;
        auto   op      = [=, &error_i](const bidx_t i0, const bidx_t i1, const bidx_t i2, const CartBlock* block) -> real_t {
            const bool no_store = (error == nullptr);
            // get the data pointers
            const ConstMemData data_field = block->data(field, ida);
            const real_t       data_sol   = (*sol)(i0, i1, i2, block);

            // we need to discard the physical BC for the edges
            const real_t value = data_field(i0, i1, i2) - data_sol;
            m_assert(std::isfinite(value), "the error cannot be nan: @ %d %d %d: %f", i0, i1, i2, data_field(i0, i1, i2));

            if (!no_store) {
                const MemData data_error = block->data(error, ida);
                data_error(i0, i1, i2)   = value;
            }

            // ei is camputed and updated
            const bool is_erri_valid = (span_.start[0] <= i0) && (i0 < span_.end[0]) &&
                                       (span_.start[1] <= i1) && (i1 < span_.end[1]) &&
                                       (span_.start[2] <= i2) && (i2 < span_.end[2]);
            error_i = m_max(std::fabs(value)*is_erri_valid, error_i);
            // error 2 is returned to be integrated
            const real_t error2 = value * value;
            return  error2 ;
        };
        // integrate error 2 and get the errori
        real_t error_2 = 0.0;
        if (level >= 0) {
            ComputeIntegral(grid, level, op, &error_2);
        } else {
            ComputeIntegral(grid, op, &error_2);
        }

        // finalize the computation
        if (norm_i != nullptr) {
            MPI_Allreduce(&error_i,norm_i+ida,1,M_MPI_REAL,MPI_SUM,MPI_COMM_WORLD);
        }
        if (norm_2 != nullptr) {
            norm_2[ida] = sqrt(error_2);
        }
    }
    //--------------------------------------------------------------------------
}

// // declare the specialization, implement them in the cpp
// template <>
// void Error::ErrorOnGridBlock<Field>(const qid_t* qid, const CartBlock* block, const Field* fid, const Field* sol) {
//     //-------------------------------------------------------------------------
//     real_t e2 = 0.0;
//     real_t ei = 0.0;

//     for (lda_t ida = 0; ida < fid->lda(); ++ida) {
//         // get the data pointers
//         const ConstMemData data_field = block->data(fid, ida);
//         const ConstMemData data_sol   = block->data(sol, ida);

//         auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // we need to discard the physical BC for the edges
//             const real_t error = data_field(i0, i1, i2) - data_sol(i0, i1, i2);
//             m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2));
//             // update the block errors
//             e2 += error * error;
//             ei = m_max(std::fabs(error), ei);
//         };
//         for_loop(&op, span_);
//     }

//     //#pragma omp critical
//     {  // no max atomic in OpenMP
//         const real_t* hgrid = block->hgrid();
//         error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
//         error_i_ = m_max(error_i_, ei);
//     }
//     //-------------------------------------------------------------------------
// };
// template <>
// void Error::ErrorFieldOnGridBlock<Field>(const qid_t* qid, const CartBlock* block, const Field* fid, const Field* sol, const Field* error) {
//     //-------------------------------------------------------------------------
//     real_t e2 = 0.0;
//     real_t ei = 0.0;

//     for (lda_t ida = 0; ida < fid->lda(); ++ida) {
//         // get the data pointers
//         const ConstMemData data_field = block->data(fid, ida);
//         const ConstMemData data_sol   = block->data(sol, ida);
//         const MemData      data_error = block->data(error, ida);

//         auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // get the local error
//             real_t error = data_field(i0, i1, i2) - data_sol(i0, i1, i2);
//             m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2));
//             // store it
//             data_error(i0, i1, i2) = error;
//             // update the block errors
//             e2 += error * error;
//             ei = m_max(std::fabs(error), ei);
//         };
//         for_loop(&op, span_);
//     }
//     // add the result
//     const real_t* hgrid = block->hgrid();

//     //#pragma omp critical
//     {  // no max atomic in OpenMP
//         error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
//         error_i_ = m_max(error_i_, ei);
//     }
//     //-------------------------------------------------------------------------
// };

// template <>
// void Error::ErrorOnGridBlock<lambda_error_t>(const qid_t* qid, const CartBlock* block, const Field* fid, const lambda_error_t* sol) {
//     //-------------------------------------------------------------------------
//     real_t e2 = 0.0;
//     real_t ei = 0.0;

//     for (lda_t ida = 0; ida < fid->lda(); ++ida) {
//         // get the data pointers
//         const ConstMemData data_field = block->data(fid, ida);

//         auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // we need to discard the physical BC for the edges
//             const real_t sol_val  = (*sol)(i0, i1, i2, block);
//             const real_t data_val = data_field(i0, i1, i2);
//             const real_t error    = data_val - sol_val;
//             m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f vs %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2), sol_val);
//             // update the block errors
//             e2 += error * error;
//             ei = m_max(std::fabs(error), ei);
//         };
//         for_loop(&op, span_);
//     }

//     // add the result
//     const real_t* hgrid = block->hgrid();

//     //#pragma omp critical
//     {  // no max atomic in OpenMP
//         error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
//         error_i_ = m_max(error_i_, ei);
//     }
//     //-------------------------------------------------------------------------
// };

// template <>
// void Error::ErrorFieldOnGridBlock<lambda_error_t>(const qid_t* qid, const CartBlock* block, const Field* fid, const lambda_error_t* sol, const Field* error) {
//     //-------------------------------------------------------------------------
//     real_t e2 = 0.0;
//     real_t ei = 0.0;

//     for (lda_t ida = 0; ida < fid->lda(); ++ida) {
//         // get the data pointers
//         const ConstMemData data_field = block->data(fid, ida);
//         const MemData      data_error = block->data(error, ida);

//         auto op = [=, &e2, &ei](const bidx_t i0, const bidx_t i1, const bidx_t i2) -> void {
//             // get the local error
//             real_t error = data_field(i0, i1, i2) - (*sol)(i0, i1, i2, block);
//             m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field(i0, i1, i2));
//             // store it
//             data_error(i0, i1, i2) = error;
//             // update the block errors
//             e2 += error * error;
//             ei = m_max(std::fabs(error), ei);
//         };
//         for_loop(&op, span_);
//     }
//     // add the result
//     const real_t* hgrid = block->hgrid();

//     //#pragma omp critical
//     {  // no max atomic in OpenMP
//         error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);
//         error_i_ = m_max(error_i_, ei);
//     }
//     //-------------------------------------------------------------------------
// };
