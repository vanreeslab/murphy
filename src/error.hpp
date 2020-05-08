#ifndef SRC_ERROR_HPP_
#define SRC_ERROR_HPP_

#include "murphy.hpp"
#include "operator.hpp"
#include "grid.hpp"

class ErrorCalculator : public ConstOperatorFF {
    real_t error_2_;
    real_t error_i_;

   public:
    void Normi(Grid* grid, Field* field, Field* sol, real_t* norm_i);
    void Norm2(Grid* grid, Field* field, Field* sol, real_t* norm_2);
    void Norms(Grid* grid, Field* field, Field* sol, real_t* norm_2, real_t* norm_i);
    
    void ApplyConstOpFF(const qid_t* qid, GridBlock* block, const Field* fid, const Field* sol) override;
};

#endif  // SRC_ERROR_HPP