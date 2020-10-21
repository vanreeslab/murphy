#ifndef SRC_ERROR_HPP_
#define SRC_ERROR_HPP_

#include "blockoperator.hpp"
#include "field.hpp"
#include "grid.hpp"
#include "interpolator.hpp"
#include "murphy.hpp"

class ErrorCalculator : public BlockOperator{
    real_t error_2_ = 0.0;  //!< the 2 norm of the error on the grid
    real_t error_i_ = 0.0;  //!< the infinite norm of the error on the grid

   public:
    explicit ErrorCalculator();
    explicit ErrorCalculator(const InterpolatingWavelet* interp);

    void Normi(Grid* grid, const Field* field, const Field* sol, real_t* norm_i);
    void Norm2(Grid* grid, const Field* field, const Field* sol, real_t* norm_2);
    void Norms(Grid* grid, const Field* field, const Field* sol, real_t* norm_2, real_t* norm_i);
    void Norms(Grid* grid, const level_t level, const Field* field, const Field* sol, real_t* norm_2, real_t* norm_i);

    void ErrorOnGridBlock(const qid_t* qid, GridBlock* block, const Field* fid, const Field* sol);
};

#endif  // SRC_ERROR_HPP
