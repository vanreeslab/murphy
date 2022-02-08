#ifndef SRC_OPERATOR_INTEGRAL_HPP_
#define SRC_OPERATOR_INTEGRAL_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/cartblock.hpp"
#include "grid/forestgrid.hpp"
#include "operator/blockoperator.hpp"
#include "core/doop.hpp"

using integrand_t = lambda_i3_t<real_t, const CartBlock*>;

/**
 * @brief Implements the integration of an integrant on the grid.
 * 
 * The integration is done using a local reconstruction of a polynomial of order N
 * - ORDER = 0 sum of values
 * - ORDER = 2 linear interpolation = trapZ
 * - ORDER = 4 third order interpolation and integration
 * 
 * 
 * @tparam O the order of the integration method
 */
template <short_t ORDER>
class Integral : public BlockOperator {
   public:
    Integral() : BlockOperator(nullptr) {
        //----------------------------------------------------------------------
        // defines the ghost point size
        ghost_len_need_[0] = ghost_len_res_[0] + m_max(ORDER / 2 - 1, 0);
        ghost_len_need_[1] = ghost_len_res_[1] + m_max(ORDER / 2, 0);
        //----------------------------------------------------------------------
    };
    Integral(const bidx_t* ghost_len) noexcept : BlockOperator(ghost_len) {
        //----------------------------------------------------------------------
        // defines the ghost point size
        ghost_len_need_[0] = ghost_len_res_[0] + ORDER / 2 - 1;
        ghost_len_need_[1] = ghost_len_res_[1] + ORDER / 2;
        //----------------------------------------------------------------------
    };

    /**
     * @brief Loop on each block and integrate the integrant and put the result in the result array
     * 
     * @warning the user is responsible to check that the integrant is actually valid when evaluated
     * @warning the user is responsible to check the completion of the request
     */
    void ComputeIntegral(const ForestGrid* grid, const integrand_t op, real_t* result) const {
        m_begin;
        //----------------------------------------------------------------------
        // allocate a small array
        real_t local_integral = 0.0;
        // go on the blocks
        DoOpMesh(this, &Integral<ORDER>::ComputeIntegralBlock, grid, op, &local_integral);
        // get the sum over the whole grid
        MPI_Allreduce(&local_integral, result, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD); 
        //--------------------------------------------------------------------------
        m_end;
    };

    /**
     * @brief Loop on each block at level @ref il and integrate the integrant and put the result in the result array
     * 
     * @warning the user is responsible to check that the integrant is actually valid when evaluated
     * @warning the user is responsible to check the completion of the request
     */
    void ComputeIntegral(const ForestGrid* grid, const level_t il, const integrand_t op, real_t* result) const {
        m_begin;
        m_assert(0 <= il, "the level you asked is negative: %d", il);
        //----------------------------------------------------------------------
        // allocate a small array
        real_t local_integral = 0.0;
        // go on the blocks
        DoOpMeshLevel(this, &Integral<ORDER>::ComputeIntegralBlock, grid, il, op, &local_integral);
        // get the sum over the whole grid
        MPI_Allreduce(&local_integral, result, 1, M_MPI_REAL, MPI_SUM, MPI_COMM_WORLD);
        //--------------------------------------------------------------------------
        m_end;
    };

   protected:
    /**
     * @brief Perform the sum on the 
     * 
     * @param qid 
     * @param block 
     * @param field 
     * @param ida 
     * @param op 
     */
    void ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const integrand_t op, real_t* result) const {
        //----------------------------------------------------------------------
        m_assert(false, "The integral at order %d is not defined", ORDER);
        result[0] = 0.0;
        //----------------------------------------------------------------------
    };
};

template <>
void Integral<0>::ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const integrand_t op, real_t* result) const;
template <>
void Integral<2>::ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const integrand_t op, real_t* result) const;
template <>
void Integral<4>::ComputeIntegralBlock(const qid_t* qid, const CartBlock* block, const integrand_t op, real_t* result) const;

#endif