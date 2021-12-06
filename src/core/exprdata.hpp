#ifndef SRC_CORE_EXPRDATA_HPP_
#define SRC_CORE_EXPRDATA_HPP_

#include "core/data.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"

class ExprData : public Data<const real_t> {
    //--------------------------------------------------------------------------
    const lda_t         ida_ = 0;
    const lambda_expr_t* expr_;  //!< analytical expression (x,y,z,ida)
    const real_t        offset_pos_[3];
    const real_t        h_[3];

   public:
    explicit ExprData() = delete;
    explicit ExprData(const lambda_expr_t* expr, const real_t offset_pos[3], const real_t h[3], const lda_t ida) : ida_(ida),
                                                                                                                   expr_(expr),
                                                                                                                   offset_pos_{offset_pos[0], offset_pos[1], offset_pos[2]},
                                                                                                                   h_{h[0], h[1], h[2]} {};

    /**
     * @brief return a copy of data located at (i0,i1,i2) position
     * 
     */
    virtual M_INLINE const real_t at(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept override {
        //--------------------------------------------------------------------------
        const real_t x = i0 * h_[0] + offset_pos_[0];
        const real_t y = i1 * h_[1] + offset_pos_[1];
        const real_t z = i2 * h_[2] + offset_pos_[2];

        // return the value of the analytical expresssion
        return (*expr_)(x, y, z, ida_);
        //--------------------------------------------------------------------------
    };
};

#endif