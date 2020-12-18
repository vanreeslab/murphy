#ifndef SRC_ERROR_HPP_
#define SRC_ERROR_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "operator/blockoperator.hpp"
#include "wavelet/wavelet.hpp"

class ErrorCalculator : public BlockOperator {
    real_t error_2_ = 0.0;  //!< the 2 norm of the error on the grid
    real_t error_i_ = 0.0;  //!< the infinite norm of the error on the grid

   public:
    explicit ErrorCalculator();
    explicit ErrorCalculator(m_ptr<const Wavelet> interp);

    void Normi(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_i);
    void Norm2(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2);
    void Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i);
    void Norms(m_ptr<const Grid> grid, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<Field> error, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i);
    void Norms(m_ptr<const Grid> grid, const level_t level, m_ptr<const Field> field, m_ptr<const Field> sol, m_ptr<real_t> norm_2, m_ptr<real_t> norm_i);

    void ErrorOnGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, m_ptr<const Field> sol);
    void ErrorFieldOnGridBlock(m_ptr<const qid_t> qid, m_ptr<GridBlock> block, m_ptr<const Field> fid, m_ptr<const Field> sol, m_ptr<const Field> error);
};

#endif  // SRC_ERROR_HPP
