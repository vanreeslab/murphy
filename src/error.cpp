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
    const sid_t lda      = field->lda();
    const lid_t nthreads = omp_get_max_threads();
    // init the errors to compute both of them
    error_2_ = (real_t*)m_calloc(sizeof(real_t) * lda * nthreads);
    error_i_ = (real_t*)m_calloc(sizeof(real_t) * lda * nthreads);

    // apply
    ConstOperatorFF::operator()(grid, field, sol);

    // do the gathering into
    if(norm_2 != nullptr){
        (*norm_2) = 0.0;
        for (sid_t ida = 0; ida < lda; ida++) {
            for (lid_t it = 0; it < nthreads; it++) {
                (*norm_2) += error_2_[it * lda + ida];
            }
        }
        (*norm_2) = std::sqrt((*norm_2));
    }
    if(norm_i != nullptr){
        (*norm_i) = 0.0;
        for (sid_t ida = 0; ida < lda; ida++) {
            for (lid_t it = 0; it < nthreads; it++) {
                (*norm_i) = m_max(error_i_[it * lda + ida], (*norm_i));
            }
        }
    }

    // free the memory
    m_free(error_2_);
    m_free(error_i_);
    
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
    const lid_t   ithread = omp_get_thread_num();
    const real_t* hgrid   = block->hgrid();
    const real_t  vol     = hgrid[0] * hgrid[1] * hgrid[2];
    
    for (sid_t ida = 0; ida < fid->lda(); ida++) {
        // get the data pointers
        real_p data_field = block->data(fid, ida);
        real_p data_sol   = block->data(sol, ida);
        // get the correct place given the current thread and the dimension
        real_t* e2 = error_2_ + ithread * fid->lda() + ida;
        real_t* ei = error_i_ + ithread * fid->lda() + ida;
        // set the errors to 0
        (*e2) = 0.0;
        (*ei) = 0.0;

        for (lid_t i2 = 0; i2 < M_N; i2++) {
            for (lid_t i1 = 0; i1 < M_N; i1++) {
                for (lid_t i0 = 0; i0 < M_N; i0++) {
                    real_t error = data_field[m_idx(i0, i1, i2)] - data_sol[m_idx(i0, i1, i2)];

                    (*e2) += error * error;
                    (*ei) = m_max(std::fabs(error), (*ei));
                }
            }
        }
        (*e2) = (*e2) * vol;
    }
    //-------------------------------------------------------------------------
}
