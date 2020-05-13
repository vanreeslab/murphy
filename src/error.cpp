#include "error.hpp"

#include "gridblock.hpp"

/**
 * @brief returns the infinite norm of the error, see @ref ErrorCalculator::Norms()
 * 
 * @param grid the grid
 * @param field the field
 * @param sol the analytical solution field
 * @param norm_i the infinite norm
 */
void ErrorCalculator::Normi(Grid* grid, Field* field, Field* sol, real_t* norm_i) {
    Norms(grid, field, sol, nullptr, norm_i);
}

/**
 * @brief returns the 2norm of the error, see @ref ErrorCalculator::Norms()
 * 
 * @param grid the grid
 * @param field the field
 * @param sol the analytical solution
 * @param norm_2 the 2-norm
 */
void ErrorCalculator::Norm2(Grid* grid, Field* field, Field* sol, real_t* norm_2) {
    Norms(grid, field, sol, norm_2, nullptr);
}

/**
 * @brief computes both the 2-norm and the infinite norm
 * 
 * If one of the norm is not wanted, set the pointer to `nullptr`.
 * 
 * The error is computed on every dimension of the field, i.e. if `error[i]` is the error in the ith dimension and
 * `local_volume` is the local cell volume, associated to the grid spacing on each block:
 * the 2 norm is defined as `sqrt( (error[0]^2 + error[1]^2 + error[2]^2) * local_volume)` and
 * the infinite norm is defined as `max(error[0], error[1], error[2])`
 * 
 * @param grid the grid
 * @param field the field
 * @param sol the analytical solution
 * @param norm_2 the 2-norm, can be `nullptr`
 * @param norm_i the infinite norm, can be `nullptr`
 */
void ErrorCalculator::Norms(Grid* grid, Field* field, Field* sol, real_t* norm_2, real_t* norm_i) {
    m_begin;
    m_assert(field->lda() == sol->lda(), "the two fields must have the same dimension");
    //-------------------------------------------------------------------------
    error_2_ = 0.0;
    error_i_ = 0.0;
    // apply
    ConstOperatorFF::operator()(grid, field, sol);
    // do the gathering into
    if (norm_2 != nullptr) {
        MPI_Allreduce(&error_2_, norm_2, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        *norm_2 = std::sqrt(norm_2[0]);
    }
    if (norm_i != nullptr) {
        MPI_Allreduce(&error_i_, norm_i, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    }
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief computes the 2 norm and the infinite norm between two field for a given block.
 * 
 * The error is computed on every dimension of the field
 * 
 * @param qid the id of the block, see @ref qid_t
 * @param block the block
 * @param fid the field
 * @param sol the analytical solution
 */
void ErrorCalculator::ApplyConstOpFF(const qid_t* qid, GridBlock* block, const Field* fid, const Field* sol) {
    //-------------------------------------------------------------------------
    const real_t* hgrid   = block->hgrid();

    real_t e2 = 0.0;
    real_t ei = 0.0;

    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        // get the data pointers
        real_p data_field = block->data(fid, ida);
        real_p data_sol   = block->data(sol, ida);
        m_assume_aligned(data_field);
        m_assume_aligned(data_sol);
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = 0; i2 < M_N; i2++) {
            for (lid_t i1 = 0; i1 < M_N; i1++) {
                for (lid_t i0 = 0; i0 < M_N; i0++) {
                    real_t error = data_field[m_idx(i0, i1, i2)] - data_sol[m_idx(i0, i1, i2)];
                    e2 += error * error;
                    ei = m_max(std::fabs(error), ei);
                }
            }
        }
    }
    // add the result
#pragma omp atomic update
    error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);

#pragma omp critical
    { // no max atomic in OpenMP
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
}
