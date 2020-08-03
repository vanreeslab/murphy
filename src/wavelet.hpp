#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

/**
 * @brief Interpolating wavelet computation
 * 
 * The wavelet application considered here is based on the lifting/dual lifting scheme.
 * More details on the wavelet approach can be found in the wavelet tutorial (https://github.com/van-Rees-Lab/wavelet_tutorial).
 * 
 * For any couple N,Nt (with Nt <= N) the filter coefficients are computed by (https://github.com/van-Rees-Lab/wavelet_tutorial/blob/master/code/wavelet_xy.py)
 * 
 * The following filters are used:
 * - ha = compute the scaling from a given level -> coarsening
 * - ga = compute the details from a given level -> detail computation
 * - gs = reconstruction from the scaling and 0 detail of the scaling coef -> refinement
 * 
 * 
 * @warning we MUST satisfy Nt <= N
 *  
 * @tparam N the number of vanishing moments of the dual wavelet -> interpolates polynomial of N-1
 * @tparam Nt the number of vanishing moments of the primal wavelet -> preserves Nt vanishing moments in the interpolation
 */
class Wavelet : public Interpolator {
   public:
#if (M_WAVELET_N == 2) && (M_WAVELET_NT == 2)
    static constexpr lid_t  len_ha_ = 5;
    static constexpr lid_t  len_ga_ = 3;
    static constexpr lid_t  len_gs_ = 2;
    static constexpr real_t ha_[]   = {-0.125, 0.25, 0.75, 0.25, -0.125};
    static constexpr real_t ga_[]   = {-0.5, 1.0, -0.5};
    static constexpr real_t gs_[]   = {0.5, 0.5};
#elif (M_WAVELET_N == 4) && (M_WAVELET_NT == 0)
    static constexpr lid_t  len_ha_ = 1;
    static constexpr lid_t  len_ga_ = 7;
    static constexpr lid_t  len_gs_ = 4;
    static constexpr real_t ha_[]   = {1.0};
    static constexpr real_t ga_[]   = {0.0625, 0.0, -0.5625, 1.0, -0.5625, 0.0, 0.0625};
    static constexpr real_t gs_[]   = {-0.0625, 0.5625, 0.5625, -0.0625};
#elif (M_WAVELET_N == 4) && (M_WAVELET_NT == 2)
    static constexpr lid_t  len_ha_ = 9;
    static constexpr lid_t  len_ga_ = 7;
    static constexpr lid_t  len_gs_ = 4;
    static constexpr real_t ha_[]   = {0.015625, 0.0, -0.125, 0.25, 0.71875, 0.25, -0.125, 0.0, 0.015625};
    static constexpr real_t ga_[]   = {0.0625, 0.0, -0.5625, 1.0, -0.5625, 0.0, 0.0625};
    static constexpr real_t gs_[]   = {-0.0625, 0.5625, 0.5625, -0.0625};
#else
    static constexpr lid_t  len_ha_ = 0;
    static constexpr lid_t  len_ga_ = 0;
    static constexpr lid_t  len_gs_ = 0;
    static constexpr real_t ha_[]   = nullptr;
    static constexpr real_t ga_[]   = nullptr;
    static constexpr real_t gs_[]   = nullptr;
#endif

    // front ghost length
    lid_t ncoarsen_front() const override { return len_ha_ / 2; }
    lid_t nrefine_front() const override { return len_gs_ / 2 - 1; }
    lid_t ncriterion_front() const override { return len_ga_ / 2 - 1; }
    // back ghost length
    lid_t nrefine_back() const override { return len_gs_ / 2; }
    lid_t ncoarsen_back() const override { return len_ha_ / 2 - 1; }
    lid_t ncriterion_back() const override { return len_ga_ / 2; }
    /*
    * @name function overriding Interpolator class
    * @{
    */
    string Identity() const override { return "interpolating wavelet " + std::to_string(M_WAVELET_N) + "." + std::to_string(M_WAVELET_NT); }

    real_t Criterion(MemLayout* block, real_p data) override;

   protected:
    void Coarsen_(const interp_ctx_t* ctx) override;
    void Refine_(const interp_ctx_t* ctx) override;
    /* @}*/

    /*
    * @name implementation specific function
    * @{
    */
   public:
    void Details(MemLayout* block, real_p data, real_t* details_max);

   protected:
    void Detail_(const interp_ctx_t* ctx, real_t* details_max);
    /* @}*/
};

#endif  // SRC_WAVELET_HPP_