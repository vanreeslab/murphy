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
 *     |       |       |       |       |            N = 4 -> 1/16, 0, -9/16, 1, -9/16, 0, 1/16
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
 * every coarse point needs Nt/2 detail information to be updated. Those lifting info are found within the 2*(Nt/2)-1 = Nt -1 first ghost points.
 * Each of these details will in turn need N/2 scaling points to be up to date.
 * Hence, the last detail coefficient to be updates will need 2*(N/2)-1 = N - 1 points in the GP
 * 
 * The whole process needs then Nt-1 + N-1 = (N+Nt)-2 ghost points on each side
 * 
 * 
 * Refinement: we proceed to the operations in the reverse order
 * ```
 * ----c------ 0 ------c------ 0 ------c----    inverse update = lifting -> dictates Nt
 *     |  _____|_____  |  _____|_____  |        in our case, we have 0 as detail, so nothing is done
 *     | /     |     \ | /     |     \ |
 *     VV      V      VVV      V      VV
 * ----c-------d-------c-------d-------c----    inverse predict = dual-lifting -> dictates N
 *     |_____  |  _____|_____  |  _____|
 *     |     \ | /     |     \ | /     |
 *     V      VVV      V      VVV      V
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
 * @warning To ease the implementation, we decided, for the moment, to consider stupid filters to apply.
 * This increases (max factor 2) the number of operations but greatly simplifies the memory access.
 * The filter coefficients are computed using this script: https://github.com/van-Rees-Lab/wavelet_tutorial/blob/master/code/interp_coef.py
 * 
 * @warning we MUST satisfy Nt <= N
 *  
 * 
 * @tparam N the number of vanishing moments of the dual wavelet -> interpolates polynomial of N-1
 * @tparam Nt the number of vanishing moments of the primal wavelet -> preserves Nt vanishing moments in the interpolation
 */
template <int N, int Nt>
class Wavelet : public Interpolator {
   protected:
    sid_t len_ha_   = 0;
    sid_t len_ha_2_ = 0;
    sid_t len_ga_   = 0;
    sid_t len_gs_   = 0;

    real_p ha_   = nullptr;  //!< scaling analysis: 1 level coarsening
    real_p ha_2_ = nullptr;  //!< scaling analysis: 2 level coarsening
    real_p ga_   = nullptr;  //!< detail analysis: 1 level coarsening
    real_p gs_   = nullptr;  //!< detail synthesis: 1 level refinement

   public:
    Wavelet() {
        m_assert(Nt <= N, "we do not support the case of Nt=%d > N=%d", Nt, N);
        sid_t lift_len = 1 + 2 * m_max(Nt - 1, 0);  // length of the lifting filter
        sid_t dual_len = 1 + 2 * m_max(N - 1, 0);   // length of the dual filter

        len_ha_   = 1 + ((Nt == 0) ? 0 : 2 * (lift_len / 2 + dual_len / 2));
        len_ga_   = dual_len;  // details
        len_gs_   = N;
        len_ha_2_ = 1 + ((Nt == 0) ? 0 : 6 * (lift_len / 2 + dual_len / 2));//6 * m_max(N - 1 + Nt - 1, 0);  //(Nt == 0) ? 0 : 6 * ((N + Nt) - 2);

        ha_   = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_ha_));
        ha_2_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_ha_2_));
        ga_   = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_ga_));
        gs_   = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_gs_));
        ha_   = ha_ + len_ha_ / 2;
        ga_   = ga_ + len_ga_ / 2;
        gs_   = gs_ + len_gs_ / 2 - 1;
        ha_2_ = ha_2_ + len_ha_2_ / 2;

        // get the correct lifting filter given N
        if (N == 2 && Nt == 2) {
            // ha
            ha_[-2] = -1.0 / 8.0;
            ha_[-1] = +1.0 / 4.0;
            ha_[0]  = +3.0 / 4.0;
            ha_[1]  = +1.0 / 4.0;
            ha_[2]  = -1.0 / 8.0;
            // ga
            ga_[-1] = -1.0 / 2.0;
            ga_[0]  = +1.0;
            ga_[1]  = -1.0 / 2.0;
            // gs
            gs_[0] = +1.0 / 2.0;
            gs_[1] = +1.0 / 2.0;
            //ha_2
            ha_2_[-6] = 1.0 / 64.0;
            ha_2_[-5] = -1.0 / 32.0;
            ha_2_[-4] = -1.0 / 8.0;
            ha_2_[-3] = 1.0 / 32.0;
            ha_2_[-2] = 7.0 / 64.0;
            ha_2_[-1] = 1.0 / 4.0;
            ha_2_[0]  = 1.0 / 2.0;
            ha_2_[1]  = 1.0 / 4.0;
            ha_2_[2]  = 7.0 / 64.0;
            ha_2_[3]  = 1.0 / 32.0;
            ha_2_[4]  = -1.0 / 8.0;
            ha_2_[5]  = -1.0 / 32.0;
            ha_2_[6]  = 1.0 / 64.0;
        } else if (N == 4 && Nt == 0) {
            // ha
            ha_[0] = 1.0;
            // ga
            ga_[-3] = 1.0 / 16.0;
            ga_[-2] = 0.0;
            ga_[-1] = -9.0 / 16.0;
            ga_[0]  = +1.0;
            ga_[1]  = -9.0 / 16.0;
            ga_[2]  = 0.0;
            ga_[3]  = 1.0 / 16.0;
            // gs
            gs_[-1] = -1.0 / 16.0;
            gs_[0]  = 9.0 / 16.0;
            gs_[1]  = 9.0 / 16.0;
            gs_[2]  = -1.0 / 16.0;
            //ha_2
            ha_2_[0] = 1.0;
        } else if (N == 4 && Nt == 2) {
            // ha
            ha_[-4] = +1.0 / 64.0;
            ha_[-3] = 0.0;
            ha_[-2] = -1.0 / 8.0;
            ha_[-1] = +1.0 / 4.0;
            ha_[0]  = +23.0 / 32.0;
            ha_[1]  = +1.0 / 4.0;
            ha_[2]  = -1.0 / 8.0;
            ha_[3]  = 0.0;
            ha_[4]  = +1.0 / 64.0;
            // ga
            ga_[-3] = 1.0 / 16.0;
            ga_[-2] = 0.0;
            ga_[-1] = -9.0 / 16.0;
            ga_[0]  = +1.0;
            ga_[1]  = -9.0 / 16.0;
            ga_[2]  = 0.0;
            ga_[3]  = 1.0 / 16.0;
            // gs
            gs_[-1] = -1.0 / 16.0;
            gs_[0]  = 9.0 / 16.0;
            gs_[1]  = 9.0 / 16.0;
            gs_[2]  = -1.0 / 16.0;
            //ha_2_
            ha_2_[-12] = 1.0 / 4096.0;
            ha_2_[-11] = 0.0;
            ha_2_[-10] = -1.0 / 512.0;
            ha_2_[-9]  = 1.0 / 256.0;
            ha_2_[-8]  = 19.0 / 2048.0;
            ha_2_[-7]  = 1.0 / 256.0;
            ha_2_[-6]  = 9.0 / 512.0;
            ha_2_[-5]  = -1.0 / 32.0;
            ha_2_[-4]  = -449.0 / 4096.0;
            ha_2_[-3]  = 1.0 / 32.0;
            ha_2_[-2]  = 7.0 / 64.0;
            ha_2_[-1]  = 31.0 / 128.0;
            ha_2_[0]   = 461.0 / 1024.0;
            ha_2_[1]   = 31.0 / 128.0;
            ha_2_[2]   = 7.0 / 64.0;
            ha_2_[3]   = 1.0 / 32.0;
            ha_2_[4]   = -449.0 / 4096.0;
            ha_2_[5]   = -1.0 / 32.0;
            ha_2_[6]   = 9.0 / 512.0;
            ha_2_[7]   = 1.0 / 256.0;
            ha_2_[8]   = 19.0 / 2048.0;
            ha_2_[9]   = 1.0 / 256.0;
            ha_2_[10]  = -1.0 / 512.0;
            ha_2_[11]  = 0.0;
            ha_2_[12]  = 1.0 / 4096.0;
        } else {
            m_assert(false, "wavelet N=%d.Nt=%d not implemented yet", N, Nt);
        }
        m_log("Wavelet %d.%d with ga[%d], ha[%d],gs[%d],ha_2[%d]",N,Nt,len_ga_,len_ha_,len_gs_,len_ha_2_);
    }

    ~Wavelet() {
        m_begin;
        //-------------------------------------------------------------------------
        m_free(ha_ - (len_ha_ / 2));
        m_free(ha_2_ - (len_ha_2_ / 2));
        m_free(ga_ - (len_ga_ / 2));
        m_free(gs_ - (len_gs_ / 2 - 1));
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