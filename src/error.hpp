#ifndef SRC_ERROR_HPP_
#define SRC_ERROR_HPP_

#include "grid.hpp"
#include "murphy.hpp"
#include "operator.hpp"

class ErrorCalculator : public ConstOperatorFF {
    lid_t  start_;
    lid_t  end_;
    real_t error_2_;
    real_t error_i_;

   public:
    explicit ErrorCalculator();
    explicit ErrorCalculator(const Grid* grid);

    void Normi(Grid* grid, Field* field, Field* sol, real_t* norm_i);
    void Norm2(Grid* grid, Field* field, Field* sol, real_t* norm_2);
    void Norms(Grid* grid, Field* field, Field* sol, real_t* norm_2, real_t* norm_i);

    // override the ConstOperatorFF function
    void ApplyConstOpFF(const qid_t* qid, GridBlock* block, const Field* fid, const Field* sol) override;
};

#endif  // SRC_ERROR_HPP
