#ifndef SRC_ERROR_HPP_
#define SRC_ERROR_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/field.hpp"
#include "grid/forestgrid.hpp"
#include "operator/integral.hpp"
#include "wavelet/wavelet.hpp"

using lambda_error_t = lambda_i3_t<real_t, const CartBlock*>;

/**
 * @brief Computes the error of a field against the given solution.
 * 
 * @tparam Sol the solution to evaluate: typically a field or a lambda fct to evaluate
 */
class Error : public Integral<2> {
   public:
    Error() noexcept : Integral<2>(){};
    Error(const bidx_t* ghost_len) noexcept : Integral<2>(ghost_len){};

    // no error storage - norm i only
    template <class Sol>
    void Normi(const ForestGrid* grid, const Field* field, const Sol* sol, real_t* norm_i) const {
        Norms(grid, field, sol, nullptr, norm_i);
    };
    // no error storage - norm2 only
    template <class Sol>
    void Norm2(const ForestGrid* grid, const Field* field, const Sol* sol, real_t* norm_2) const {
        Norms(grid, field, sol, norm_2, nullptr);
    };

    // no error storage - all norms
    template <class Sol>
    void Norms(const ForestGrid* grid, const Field* field, const Sol* sol, real_t* norm_2, real_t* norm_i) const {
        ErrorMagic(grid, -1, field, sol, nullptr, norm_2, norm_i);
    };

    // no error storage - all norms but only on a given level
    template <class Sol>
    void Norms(const ForestGrid* grid, const level_t level, const Field* field, const Sol* sol, real_t* norm_2, real_t* norm_i) const {
        ErrorMagic(grid, level, field, sol, nullptr, norm_2, norm_i);
    };

    // store the error field
    template <class Sol>
    void Norms(const ForestGrid* grid, const Field* field, const Sol* sol, Field* error, real_t* norm_2, real_t* norm_i) const {
        ErrorMagic(grid, -1, field, sol, error, norm_2, norm_i);
    };

   protected:
    template <class Sol>
    void ErrorMagic(const ForestGrid* grid, const level_t level, const Field* field, const Sol* sol, Field* error, real_t* norm_2, real_t* norm_i) const {
        m_begin;
        m_assert(IsGhostValid(field), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
        m_assert(IsGhostValid(sol), "we cannot compute the ghost, please get the ghost before for field <%s>", field->name().c_str());
        //-------------------------------------------------------------------------
        m_assert(false, "No overloading for the type of solution you provided");
        // // call the operator
        // if (no_error && no_level) {
        //     DoOpMesh(this, &Error::ErrorOnGridBlock<Sol>, grid, field, sol);
        // } else if (!no_error && no_level) {
        //     DoOpMesh(this, &Error::ErrorFieldOnGridBlock<Sol>, grid, field, sol, error);
        // } else if (no_error && !no_level) {
        //     m_assert(level >= 0 && level < P8EST_MAXLEVEL, "the level = %d must be >= 0 and < %d", level, P8EST_MAXLEVEL);
        //     DoOpMeshLevel(this, &Error::ErrorOnGridBlock<Sol>, grid, level, field, sol);
        // } else {  // error && level
        //     m_assert(level >= 0 && level < P8EST_MAXLEVEL, "the level = %d must be >= 0 and < %d", level, P8EST_MAXLEVEL);
        //     DoOpMeshLevel(this, &Error::ErrorFieldOnGridBlock<Sol>, grid, level, field, sol, error);
        // }

        // if (!no_error) {
        //     error->ghost_len(ghost_len_res_);
        // }

        // // do the gathering into
        // if (!(norm_2 == nullptr)) {
        //     norm_2[0] = 0.0;
        //     MPI_Allreduce(&error_2_, norm_2, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        //     norm_2[0] = std::sqrt(norm_2[0]);
        // }
        // if (!(norm_i == nullptr)) {
        //     norm_i[0] = 0.0;
        //     MPI_Allreduce(&error_i_, norm_i, 1, M_MPI_REAL, MPI_MAX, MPI_COMM_WORLD);
        // }
        //-------------------------------------------------------------------------
        m_end;
    };

    // template <class Sol>
    // void ErrorOnGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid, const Sol* sol) {
    //     //-------------------------------------------------------------------------
    //     m_assert(false, "Function needs to be specialized: sol");
    //     //-------------------------------------------------------------------------
    // };

    // template <class Sol>
    // void ErrorFieldOnGridBlock(const qid_t* qid, const CartBlock* block, const Field* fid, const Sol* sol, const Field* error) {
    //     //-------------------------------------------------------------------------
    //     m_assert(false, "Function needs to be specialized");
    //     //-------------------------------------------------------------------------
    // };
};

template <>
void Error::ErrorMagic<Field>(const ForestGrid* grid, const level_t level, const Field* field, const Field* sol, Field* error, real_t* norm_2, real_t* norm_i) const;
template <>
void Error::ErrorMagic<lambda_error_t>(const ForestGrid* grid, const level_t level, const Field* field, const lambda_error_t* sol, Field* error, real_t* norm_2, real_t* norm_i) const;

// declare the specialization, implement them in the cpp
// template <>
// void Error::ErrorOnGridBlock<Field>(const qid_t* qid, const CartBlock* block, const Field* fid, const Field* sol);
// template <>
// void Error::ErrorFieldOnGridBlock<Field>(const qid_t* qid, const CartBlock* block, const Field* fid, const Field* sol, const Field* error);

// template <>
// void Error::ErrorOnGridBlock<lambda_error_t>(const qid_t* qid, const CartBlock* block, const Field* fid, const lambda_error_t* sol);
// template <>
// void Error::ErrorFieldOnGridBlock<lambda_error_t>(const qid_t* qid, const CartBlock* block, const Field* fid, const lambda_error_t* sol, const Field* error);

#endif  // SRC_ERROR_HPP
