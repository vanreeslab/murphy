#ifndef SRC_ERROR_HPP_
#define SRC_ERROR_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "operator/blockoperator.hpp"
#include "wavelet/wavelet.hpp"

/**
 * @brief Computes the error of a field against the given solution.
 * 
 * @tparam Sol the solution to evaluate: typically a field or a lambda fct to evaluate
 */
class Error : public BlockOperator {
    real_t error_2_ = 0.0;  //!< the 2 norm of the error on the grid
    real_t error_i_ = 0.0;  //!< the infinite norm of the error on the grid

   public:
    explicit Error() : BlockOperator(nullptr){};
    explicit Error(m_ptr<const Wavelet> interp) : BlockOperator(interp){};

    template <class Sol>
    void Normi(m_ptr<const Grid> grid, m_ptr<const Field> field, Sol sol, m_ptr<real_t> norm_i) {
        Norms(grid, field, sol, nullptr, norm_i);
    };

    template <class Sol>
    void Norm2(m_ptr<const Grid> grid, m_ptr<const Field> field, Sol sol, m_ptr<real_t> norm_2) {
        Norms(grid, field, sol, norm_2, nullptr);
    };

    template <class Sol>
    void Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, Sol sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
        Norms(grid, -1, field, sol, nullptr, norm_2, norm_i);
    };

    template <class Sol>
    void Norms(m_ptr<const Grid> grid, const level_t level, m_ptr<const Field> field, Sol sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
        Norms(grid, level, field, sol, nullptr, norm_2, norm_i);
    };

    template <class Sol>
    void Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, Sol sol, m_ptr<Field> error, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
        Norms(grid, -1, field, sol, error, norm_2, norm_i);
    };

    template <class Sol>
    void Norms(m_ptr<const Grid> grid, const level_t level, m_ptr<const Field> field, Sol sol, m_ptr<Field> error, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i) {
        m_begin;
        m_assert(!(do_ghost_ && (!field->ghost_status())), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
        //-------------------------------------------------------------------------
        error_2_ = 0.0;
        error_i_ = 0.0;

        const bool no_level = (level == (-1));
        const bool no_error = (error.IsEmpty());

        // call the operator
        if (no_error && no_level) {
            DoOpMesh(this, &Error::ErrorOnGridBlock<Sol>, grid, field, sol);
        } else if (!no_error && no_level) {
            m_ptr<const Field> error_cst = error;
            DoOpMesh(this, &Error::ErrorFieldOnGridBlock<Sol>, grid, field, sol, error_cst);
        } else if (no_error && !no_level) {
            m_assert(level >= 0 && level < P8EST_MAXLEVEL, "the level = %d must be >= 0 and < %d", level, P8EST_MAXLEVEL);
            DoOpMeshLevel(this, &Error::ErrorOnGridBlock<Sol>, grid, level, field, sol);
        } else {  // error && level
            m_assert(level >= 0 && level < P8EST_MAXLEVEL, "the level = %d must be >= 0 and < %d", level, P8EST_MAXLEVEL);
            m_ptr<const Field> error_cst = error;
            DoOpMeshLevel(this, &Error::ErrorFieldOnGridBlock<Sol>, grid, level, field, sol, error_cst);
        }

        if (!no_error) {
            error->ghost_status(do_ghost());
        }

        // do the gathering into
        if (!norm_2.IsEmpty()) {
            real_t recv_buf = 0.0;
            MPI_Allreduce(&error_2_, &recv_buf, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
            norm_2[0] = std::sqrt(recv_buf);
        }
        if (!norm_i.IsEmpty()) {
            real_t recv_buf = 0.0;
            MPI_Allreduce(&error_i_, &recv_buf, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
            norm_i[0] = recv_buf;
        }
        //-------------------------------------------------------------------------
        m_end;
    };

    template <class Sol>
    void ErrorOnGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, Sol sol) {
        //-------------------------------------------------------------------------
        m_assert(false, "Function needs to be specialized: sol");
        //-------------------------------------------------------------------------
    };

    template <class Sol>
    void ErrorFieldOnGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, Sol sol, m_ptr<const Field> error) {
        //-------------------------------------------------------------------------
        m_assert(false, "Function needs to be specialized");
        //-------------------------------------------------------------------------
    };
};

// declare the specialization, implement them in the cpp
template <>
void Error::ErrorOnGridBlock<m_ptr<const Field> >(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, m_ptr<const Field> sol);
template <>
void Error::ErrorFieldOnGridBlock<m_ptr<const Field> >(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, m_ptr<const Field> sol, m_ptr<const Field> error);

template <>
void Error::ErrorOnGridBlock<lambda_i3block_t*>(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, lambda_i3block_t* sol);
template <>
void Error::ErrorFieldOnGridBlock<lambda_i3block_t*>(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, lambda_i3block_t* sol, m_ptr<const Field> error);

#endif  // SRC_ERROR_HPP
