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
    explicit Error(const Wavelet*  interp) : BlockOperator(interp){};

    template <class Sol>
    void Normi(const Grid*  grid, const Field*  field, Sol sol, real_t*  norm_i) {
        Norms(grid, field, sol, nullptr, norm_i);
    };

    template <class Sol>
    void Norm2(const Grid*  grid, const Field*  field, Sol sol, real_t*  norm_2) {
        Norms(grid, field, sol, norm_2, nullptr);
    };

    template <class Sol>
    void Norms(const Grid*  grid, const Field*  field, Sol sol, real_t*  norm_2, real_t*  norm_i) {
        Norms(grid, -1, field, sol, nullptr, norm_2, norm_i);
    };

    template <class Sol>
    void Norms(const Grid*  grid, const level_t level, const Field*  field, Sol sol, real_t*  norm_2, real_t*  norm_i) {
        Norms(grid, level, field, sol, nullptr, norm_2, norm_i);
    };

    template <class Sol>
    void Norms(const Grid*  grid, const Field*  field, Sol sol, Field*  error, real_t*  norm_2, real_t*  norm_i) {
        Norms(grid, -1, field, sol, error, norm_2, norm_i);
    };

    template <class Sol>
    void Norms(const Grid*  grid, const level_t level, const Field*  field, Sol sol, Field*  error, real_t*  norm_2, real_t*  norm_i) {
        m_begin;
        m_assert(!(do_ghost_ && (!field->ghost_status())), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
        //-------------------------------------------------------------------------
        error_2_ = 0.0;
        error_i_ = 0.0;

        const bool no_level = level == (-1);
        const bool no_error = (error == nullptr);

        // call the operator
        if (no_error && no_level) {
            DoOpMesh(this, &Error::ErrorOnGridBlock<Sol>, grid, field, sol);
        } else if (!no_error && no_level) {
            const Field* error_cst = error;
            DoOpMesh(this, &Error::ErrorFieldOnGridBlock<Sol>, grid, field, sol, error_cst);
        } else if (no_error && !no_level) {
            m_assert(level >= 0 && level < P8EST_MAXLEVEL, "the level = %d must be >= 0 and < %d", level, P8EST_MAXLEVEL);
            DoOpMeshLevel(this, &Error::ErrorOnGridBlock<Sol>, grid, level, field, sol);
        } else {  // error && level
            m_assert(level >= 0 && level < P8EST_MAXLEVEL, "the level = %d must be >= 0 and < %d", level, P8EST_MAXLEVEL);
            const Field* error_cst = error;
            DoOpMeshLevel(this, &Error::ErrorFieldOnGridBlock<Sol>, grid, level, field, sol, error_cst);
        }

        if (!no_error) {
            error->ghost_status(do_ghost());
        }

        // do the gathering into
        if (!(norm_2 == nullptr)) {
            norm_2[0] = 0.0;
            MPI_Allreduce(&error_2_, norm_2, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
            norm_2[0] = std::sqrt(norm_2[0]);
        }
        if (!(norm_i == nullptr)) {
            norm_i[0] = 0.0;
            MPI_Allreduce(&error_i_, norm_i, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
        }
        //-------------------------------------------------------------------------
        m_end;
    };

    template <class Sol>
    void ErrorOnGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid, Sol sol) {
        //-------------------------------------------------------------------------
        m_assert(false, "Function needs to be specialized: sol");
        //-------------------------------------------------------------------------
    };

    template <class Sol>
    void ErrorFieldOnGridBlock(const qid_t*  qid, GridBlock*  block, const Field*  fid, Sol sol, const Field*  error) {
        //-------------------------------------------------------------------------
        m_assert(false, "Function needs to be specialized");
        //-------------------------------------------------------------------------
    };
};

// declare the specialization, implement them in the cpp
template <>
void Error::ErrorOnGridBlock<const Field*  >(const qid_t*  qid, GridBlock*  block, const Field*  fid, const Field*  sol);
template <>
void Error::ErrorFieldOnGridBlock<const Field*  >(const qid_t*  qid, GridBlock*  block, const Field*  fid, const Field*  sol, const Field*  error);

template <>
void Error::ErrorOnGridBlock<lambda_i3block_t*>(const qid_t*  qid, GridBlock*  block, const Field*  fid, lambda_i3block_t* sol);
template <>
void Error::ErrorFieldOnGridBlock<lambda_i3block_t*>(const qid_t*  qid, GridBlock*  block, const Field*  fid, lambda_i3block_t* sol, const Field*  error);

#endif  // SRC_ERROR_HPP
