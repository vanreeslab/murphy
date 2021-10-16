#ifndef SRC_OPERATOR_DIAGNOSTICS_HPP
#define SRC_OPERATOR_DIAGNOSTICS_HPP

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/field.hpp"
#include "grid/grid.hpp"
#include "grid/gridblock.hpp"
#include "operator/blockoperator.hpp"
#include "wavelet/wavelet.hpp"
#include "operator/error.hpp"

/**
 * @brief Performs the diagnostics of the detail vs error
 * 
 * @note we could have heritated from Stencil, however we need all the GP at onces for the wavelets.
 * it's then easier to take it from BlockOperator
 * 
 */
class DetailVsError : public BlockOperator {
    bool do_error_ = true;
    iter_t n_cat_   = 48;
    real_t max_cat_ = 1.0;
    real_t min_cat_ = 1.0e-16;

    bidx_t* n_blocks_loc_  = nullptr;
    bidx_t* n_blocks_glob_ = nullptr;

   public:
    explicit DetailVsError() = delete;
    explicit DetailVsError(const Wavelet* interp) noexcept;

    void operator()(const iter_t id, const std::string folder, const std::string suffix, const Grid* grid, Field* field, const lambda_error_t* sol);
    void DoMagic(const qid_t* qid, GridBlock* block, const Wavelet* interp, const Field* field, const lambda_error_t* sol);
};

#endif