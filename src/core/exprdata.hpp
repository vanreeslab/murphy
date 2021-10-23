#ifndef SRC_CORE_EXPRDATA_HPP_
#define SRC_CORE_EXPRDATA_HPP_

#include "core/data.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"

class ExprData : public Data<const real_t> {
    //--------------------------------------------------------------------------
    const lda_t                       ida_ = 0;
    const lambda_i3_t<real_t, lda_t>* expr_;

   public:
    explicit ExprData() = delete;
    explicit ExprData(const lambda_i3_t<real_t, lda_t>* expr, const lda_t ida) : expr_(expr), ida_(ida){};

    //--------------------------------------------------------------------------
    /**
     * @brief return a copy of data located at (i0,i1,i2) position
     * 
     */
    virtual M_INLINE const real_t at(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept override {
        return (*expr_)(i0,  i1, i2, ida_);;
    };
};

#endif