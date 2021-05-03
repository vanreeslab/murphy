#include "error.hpp"

#include "gridblock.hpp"

/**
 * @brief Construct a new Error Calculator, if the grid is not nullptr, compute the error on the ghost points as well
 * 
 * @param interp the Wavelet to use to get the number of actual ghost points, see BlockOperator::BlockOperator()
 */
ErrorCalculator::ErrorCalculator(m_ptr<const Wavelet> interp) : BlockOperator(interp) {}
ErrorCalculator::ErrorCalculator() : BlockOperator(nullptr) {}

/**
 * @brief returns the infinite norm of the error, see @ref ErrorCalculator::Norms()
 * 
 * @param grid the grid
 * @param field the field
 * @param sol the analytical solution field
 * @param norm_i the infinite norm
 */
void ErrorCalculator::Normi(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_i) {
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
void ErrorCalculator::Norm2(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2) {
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
void ErrorCalculator::Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
    m_begin;
    //-------------------------------------------------------------------------
    Norms(grid, field, sol, nullptr, norm_2, norm_i);
    //-------------------------------------------------------------------------
    m_end;
}

void ErrorCalculator::Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<Field> error, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
    m_begin;
    m_assert(field->lda() == sol->lda(), "the two fields must have the same dimension");
    //-------------------------------------------------------------------------
    error_2_ = 0.0;
    error_i_ = 0.0;

    // check that if we do the ghost, the ghost are updated
    m_assert(!(do_ghost_ && (!field->ghost_status())), "we cannot compute the ghost");
    // call the operator
    if (error.IsEmpty()) {
        DoOpMesh(this, &ErrorCalculator::ErrorOnGridBlock, grid, field, sol);
    } else {
        m_ptr<const Field> error_cst = error;
        DoOpMesh(this, &ErrorCalculator::ErrorFieldOnGridBlock, grid, field, sol, error_cst);
    }
    if (!error.IsEmpty()) {
        error->ghost_status(do_ghost());
    }

    // do the gathering into
    if (!norm_2.IsEmpty()) {
        (*norm_2()) = 0.0;
        MPI_Allreduce(&error_2_, norm_2(), 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        (*norm_2()) = std::sqrt(norm_2()[0]);
    }
    if (!norm_i.IsEmpty()) {
        (*norm_i()) = 0.0;
        MPI_Allreduce(&error_i_, norm_i(), 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
    }
}

/**
 * @brief same as Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) but operates on the given level only
 * 
 * @param grid 
 * @param level the operation level
 * @param field 
 * @param sol 
 * @param norm_2 
 * @param norm_i 
 */
void ErrorCalculator::Norms(m_ptr<const Grid> grid, const level_t level, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
    m_begin;
    m_assert(field->lda() == sol->lda(), "the two fields must have the same dimension");
    //-------------------------------------------------------------------------
    error_2_ = 0.0;
    error_i_ = 0.0;

    // check that if we do the ghost, the ghost are updated
    m_assert(!(do_ghost_ && (!field->ghost_status())), "we cannot compute the ghost");
    // call the operator
    DoOpMeshLevel(this, &ErrorCalculator::ErrorOnGridBlock, grid, level, field, sol);

    // do the gathering into
    if (!norm_2.IsEmpty()) {
        (*norm_2()) = 0.0;
        MPI_Allreduce(&error_2_, norm_2(), 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        (*norm_2()) = std::sqrt(norm_2()[0]);
    }
    if (!norm_i.IsEmpty()) {
        (*norm_i()) = 0.0;
        MPI_Allreduce(&error_i_, norm_i(), 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
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
void ErrorCalculator::ErrorOnGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, m_ptr<const Field> sol) {
    //-------------------------------------------------------------------------
    const real_t* hgrid = block->hgrid();

    real_t e2 = 0.0;
    real_t ei = 0.0;

    // m_log("compute error %d to %d for block @ %f %f %f", start_, end_, block->xyz(0), block->xyz(1), block->xyz(2));

    const_data_ptr block_field = block->data(fid);
    const_data_ptr block_sol   = block->data(sol);

    for (sid_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers

        const real_t* data_field = block_field.Read(ida, block());
        const real_t* data_sol   = block_sol.Read(ida, block());

        // m_assume_aligned(data_field);
        // m_assume_aligned(data_sol);
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    real_t pos[3];
                    m_pos(pos, i0, i1, i2, hgrid, block->xyz());
                    // we need to discard the physical BC for the edges

                    real_t error = data_field[m_idx(i0, i1, i2)] - data_sol[m_idx(i0, i1, i2)];

                    m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field[m_idx(i0, i1, i2)]);
                    e2 += error * error;
                    ei = m_max(std::fabs(error), ei);

                    // if (fabs(error) > 7.8) {
                    //     m_log("block %d @ %f %f %f (ida=%d) , @ %d %d %d field = %e, sol = %e, error = %e", qid->cid, block->xyz(0), block->xyz(1), block->xyz(2), ida, i0, i1, i2, data_field[m_idx(i0, i1, i2)], data_sol[m_idx(i0, i1, i2)], error);
                    // }
                }
            }
        }
    }
    // add the result
    //#pragma omp atomic update
    error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);

    //#pragma omp critical
    {  // no max atomic in OpenMP
        error_i_ = m_max(error_i_, ei);
    }

    // m_log("compute error done: err2 = %e, erri = %e (local: %e %e)", error_2_, error_i_, e2, ei);

    //-------------------------------------------------------------------------
}

/**
 * @brief computes the 2 norm and the infinite norm between two field for a given block and fill the error field with the error
 * 
 * The error is computed on every dimension of the field
 * 
 * @param qid the id of the block, see @ref qid_t
 * @param block the block
 * @param fid the field
 * @param sol the analytical solution
 * @param error the error field to fill
 */
void ErrorCalculator::ErrorFieldOnGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, m_ptr<const Field> sol, m_ptr<const Field> error) {
    //-------------------------------------------------------------------------
    const real_t* hgrid = block->hgrid();

    real_t e2 = 0.0;
    real_t ei = 0.0;

    const_data_ptr ptr_field = block->data(fid);
    const_data_ptr ptr_sol   = block->data(sol);
    data_ptr       ptr_error = block->data(error);

    for (sid_t ida = 0; ida < fid->lda(); ++ida) {
        // get the data pointers

        const real_t* data_field = ptr_field.Read(ida, block());
        const real_t* data_sol   = ptr_sol.Read(ida, block());
        real_t*       data_error = ptr_error.Write(ida, block());
        // get the correct place given the current thread and the dimension
        for (lid_t i2 = start_; i2 < end_; i2++) {
            for (lid_t i1 = start_; i1 < end_; i1++) {
                for (lid_t i0 = start_; i0 < end_; i0++) {
                    real_t pos[3];
                    m_pos(pos, i0, i1, i2, hgrid, block->xyz());
                    // we need to discard the physical BC for the edges

                    real_t error = data_field[m_idx(i0, i1, i2)] - data_sol[m_idx(i0, i1, i2)];

                    data_error[m_idx(i0, i1, i2)] = error;

                    m_assert(error == error, "the error cannot be nan: tree %d block %d @ %d %d %d: %f", qid->tid, qid->qid, i0, i1, i2, data_field[m_idx(i0, i1, i2)]);
                    e2 += error * error;
                    ei = m_max(std::fabs(error), ei);
                }
            }
        }
    }
    // add the result
    //#pragma omp atomic update
    error_2_ += e2 * (hgrid[0] * hgrid[1] * hgrid[2]);

    //#pragma omp critical
    {  // no max atomic in OpenMP
        error_i_ = m_max(error_i_, ei);
    }
    //-------------------------------------------------------------------------
}
