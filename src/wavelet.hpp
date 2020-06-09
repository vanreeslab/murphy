#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

/**
 * @brief Interpolating wavelet computation
 * 
 * More details on the wavelet coefficients can be found in the wavelet tutorial (https://github.com/van-Rees-Lab/wavelet_tutorial).
 * 
 * 
 * The application of the wavelet is based on Mallat 2009, as a succession of split, dual lifting and primal lifting.
 * It is said that this approach reduces by a factor of 2 the number of computations (hence 8 in 3D?).
 * 
 * The coarsening operation is driven by N and Nt, the templates param:
 * ```
 * ----c-------c-------c-------c-------c----    even-odd split
 *     |       |       |       |       |
 *   1 |     1 |     1 |     1 |     1 |
 *     V       V       V       V       V
 * ----c-------d-------c-------d-------c----    predict = dual lifting -> dictates N
 *     |\______|______/|\______|______/|            N = 2 -> -1/2, 1 ,-1/2
 *     |       |       |       |       |            N = 4 -> 1/16, -9/16, 1, -9/16, 1/16
 *     V       V       V       V       V
 * ----c-------d-------c-------d-------c----    update = lifting -> dictates Nt
 *     |______/ \______|______/ \______|            Nt = 2 -> +1/4, 1, +1/4
 *     |               |               |            Nt = 4 -> -1/32,  9/32, 9/32, -1/32
 *     V               V               V
 * ----c-------d-------c-------d-------c----
 * ```
 * We can check that the combined application of N=4 and Nt=2, we obtain the same results as table 2 in Sweldens 1996.
 * 
 * Numbef of needed ghost points, i.e. the number of finer info on the side of the last interessing point
 * every coarse point needs Nt/2 detail information to be updated. Those lifting info are found within the 2*(Nt/2)-1 first ghost points.
 * Each of these details will in turn need N/2 scaling points to be up to date.
 * Hence, the last detail coefficient to be updates will need 2*(N/2)-1 points in the GP
 * 
 * The whole process needs then 2*(Nt/2)-1 + 2*(N/2)-1 = (N+Nt)-2 ghost points
 * 
 * 
 * Refinement: we proceed to the operations in the reverse order
 * ```
 * ----c---------------c---------------c----    inverse update = lifting -> dictates Nt
 *     |_____     _____|_____     _____|
 *     |     \   /     |     \   /     |
 *     V      V V      V      V V      V
 * ----c-------d-------c-------d-------c----    inverse predict = dual-lifting -> dictates N
 *     |  _____|_____  |  _____|_____  |
 *     | /     |     \ | /     |     \ |
 *     VV      V      VVV      V      VV
 * ----c-------d-------c-------d-------c----    inverse split
 *     |       |       |       |       |
 *   1 |     1 |     1 |     1 |     1 |
 *     V       V       V       V       V
 * ----c-------c-------c-------c-------c----
 * ```
 * 
 * Number of needed ghost points, i.e. the number of coarser info on the side of the last interessing point
 * every fine info will need N/2 detail coefficient to be updated, which are among the (2*N/2)-1 first fine ghost points.
 * Then, each detail needs Nt/2 scaling coef to be updated, found among the (Nt/2-1) fine information on its left.
 * Finally, we then need (N+Nt)-2 fine information, which is equivalent to (N+Nt)/2-1 coarse information
 *  
 * 
 * @tparam N the number of vanishing moments of the dual wavelet -> interpolates polynomial of N-1
 * @tparam Nt the number of vanishing moments of the primal wavelet -> preserves Nt vanishing moments in the interpolation
 */
template <int N, int Nt>
class Wavelet : public Interpolator {
   protected:
    real_t* s_lift_;
    real_t* s_dual_;

   public:
    Wavelet() {
        s_dual_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * (2 * N - 1)));
        s_dual_ = s_dual_ + N - 1;
        s_lift_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * (2 * Nt - 1)));
        s_lift_ = s_lift_ + Nt - 1;
        // get the correct lifting filter given N
        if (N == 2) {
            // lifting coefficients
            s_dual_[-1] = 1.0 / 2.0;
            s_dual_[0]  = 1.0;
            s_dual_[1]  = 1.0 / 2.0;
        } else if (N == 4) {
            // lifting coefficients
            s_dual_[-3] = +1.0 / 16.0;
            s_dual_[-2] = 0.0;
            s_dual_[-1] = -9.0 / 16.0;
            s_dual_[0]  = 1.0;
            s_dual_[1]  = -9.0 / 16.0;
            s_dual_[2]  = 0.0;
            s_dual_[3]  = +1.0 / 16.0;
        } else {
            m_assert(false, "wavelet N=%d not implemented yet", N);
        }
        // get the correct dual lifting filter given Nt
        if (Nt == 2) {
            // lifting coefficients
            s_lift_[-1] = 1.0 / 4.0;
            s_lift_[0]  = 1.0;
            s_lift_[1]  = 1.0 / 4.0;
        } else if (Nt == 4) {
            // lifting coefficients
            s_lift_[-3] = +1.0 / 32.0;
            s_lift_[-2] = 0.0;
            s_lift_[-1] = -9.0 / 32.0;
            s_lift_[0]  = 1.0;
            s_lift_[1]  = -9.0 / 32.0;
            s_lift_[2]  = 0.0;
            s_lift_[3]  = +1.0 / 32.0;
        } else {
            m_assert(false, "wavelet Nt=%d not implemented yet", Nt);
        }
    }

    ~Wavelet() {
        m_begin;
        //-------------------------------------------------------------------------

        //-------------------------------------------------------------------------
        m_end;
    }

    real_t Criterion(MemLayout* block, real_p data, MemPool* mem_pool) override;
    void   Details(MemLayout* block, real_p data, real_t* criterion);

    string Identity() override { return "interpolating wavelet" + std::to_string(N) + "." + std::to_string(Nt); }
    lid_t  NGhostCoarse() const override { return (N + Nt) / 2 - 1; }
    lid_t  NGhostFine() const override { return N + Nt - 2; }

   protected:
    void Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const override;
    void Refine_(const interp_ctx_t* ctx) const override;
    void Copy_(const interp_ctx_t* ctx) const override;
    void Detail_(const interp_ctx_t* ctx, real_t* details_inf_norm) const;
};

#endif  // SRC_WAVELET_HPP_

#include "wavelet.inc"