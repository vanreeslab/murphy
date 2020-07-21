#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include "interpolator.hpp"
#include "murphy.hpp"

#define M_WFRONT 0
#define M_WBACK 1
#define M_WCOAR 0
#define M_WFINE 1

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
 * ----c-------0-------c-------0-------c----
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
 * This increases (about a factor 2) the number of operations but greatly simplifies the memory access.
 * The filter coefficients, lengths, etc are given by this script: https://github.com/van-Rees-Lab/wavelet_tutorial/blob/master/code/interp_coef.py
 * 
 * @warning we MUST satisfy Nt <= N
 *  
 * @tparam N the number of vanishing moments of the dual wavelet -> interpolates polynomial of N-1
 * @tparam Nt the number of vanishing moments of the primal wavelet -> preserves Nt vanishing moments in the interpolation
 */
template <int N, int Nt>
class Wavelet : public Interpolator {
   protected:
    sid_t len_ha_ = 0;
    sid_t len_ga_ = 0;
    sid_t len_gs_ = 0;

    // sid_t n_ghost_[2]   = {0, 0}; //!< number of ghost needed by the wavelet computation on  for a block: i.e. [front: normal = -1, back: normal=+1]
    // sid_t n_info_[2] = {{0, 0}, {0, 0}};  //!< for each block side, the number of coarse points and fine points needed
    // standard filters
    real_t* ha_ = nullptr;  //!< scaling analysis: 1 level coarsening -> coarsening computation
    real_t* ga_ = nullptr;  //!< detail analysis: 1 level coarsening -> detail computation
    real_t* gs_ = nullptr;  //!< detail synthesis: 1 level refinement -> refinement computation
    // // modified filter
    // real_t** gs_g_[2][2] = {{nullptr, nullptr}, {nullptr, nullptr}};  //!< modified gs filter for the ghosting[front,back][coarse,fine][ighost]

    // real_p ha_2_ = nullptr;  //!< scaling analysis: 2 level coarsening
    // sid_t len_ha_2_ = 0;

    real_t** temp_;

   public:
    Wavelet() {
        m_assert(Nt <= N, "we do not support the case of Nt=%d > N=%d", Nt, N);
        //-------------------------------------------------------------------------

        // get the size of the filters and the number of ghost points
        if (N == 2 && Nt == 2) {
            // length of the filters
            len_ha_ = 5;
            len_ga_ = 3;
            len_gs_ = 2;
            // number of ghosts
            // n_ghost_[M_WFRONT]         = 2;  // total number of ghost @ front
            // n_ghost_[M_WBACK]          = 1;  // total number of ghost @ back
            // n_info_[M_WFRONT][M_WCOAR] = 1;  // front, n_coarse
            // n_info_[M_WFRONT][M_WFINE] = 3;  // front, n_fine
            // n_info_[M_WBACK][M_WCOAR]  = 0;  // back, n_coarse
            // n_info_[M_WBACK][M_WFINE]  = 0;  // back, n_fine
        } else if (N == 4 && Nt == 0) {
            // length of the filters
            len_ha_ = 1;
            len_ga_ = 7;
            len_gs_ = 4;
            // number of ghosts
            // n_ghost_[M_WFRONT]         = 2;
            // n_ghost_[M_WBACK]          = 3;
            // n_info_[M_WFRONT][M_WCOAR] = 2;  // front, n_coarse
            // n_info_[M_WFRONT][M_WFINE] = 3;  // front, n_fine
            // n_info_[M_WBACK][M_WCOAR]  = 3;  // back, n_coarse
            // n_info_[M_WBACK][M_WFINE]  = 2;  // back, n_fine
        } else if (N == 4 && Nt == 2) {
            // length of the filters
            len_ha_ = 9;
            len_ga_ = 7;
            len_gs_ = 4;
            // number of ghosts
            // n_ghost_[0]   = 4;
            // n_ghost_[1]   = 3;
            // n_info_[0][0] = 3;  // front, n_coarse
            // n_info_[0][1] = 7;  // front, n_fine
            // n_info_[1][0] = 3;  // back, n_coarse
            // n_info_[1][1] = 6;  // back, n_fine
        } else {
            m_assert(false, "wavelet N=%d.Nt=%d not implemented yet", N, Nt);
        }

        // allocate the memory for the filters
        ha_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_ha_));
        ga_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_ga_));
        gs_ = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * len_gs_));
        // ha_   = ha_ + len_ha_ / 2;
        // ga_   = ga_ + len_ga_ / 2;
        // gs_   = gs_ + len_gs_ / 2 - 1;
        // ha_2_ = ha_2_ + len_ha_2_ / 2;

        // ghost modified filters, only used for n_ghost/2 detail points
        // sid_t n_ghost_det = n_ghost_[M_WFRONT] / 2;
        // if (n_ghost_det > 0) {
        //     gs_g_[M_WFRONT][M_WCOAR] = reinterpret_cast<real_t**>(m_calloc(sizeof(real_t*) * n_ghost_det));
        //     gs_g_[M_WFRONT][M_WFINE] = reinterpret_cast<real_t**>(m_calloc(sizeof(real_t*) * n_ghost_det));
        //     for (sid_t ig = 0; ig < n_ghost_det; ig++) {
        //         gs_g_[M_WFRONT][M_WCOAR][ig] = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * n_info_[M_WFRONT][M_WCOAR]));
        //         gs_g_[M_WFRONT][M_WFINE][ig] = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * n_info_[M_WFRONT][M_WFINE]));
        //     }
        // }
        // n_ghost_det = n_ghost_[M_WBACK] / 2;
        // if (n_ghost_det > 0) {
        //     gs_g_[M_WBACK][M_WCOAR] = reinterpret_cast<real_t**>(m_calloc(sizeof(real_t*) * n_ghost_det));
        //     gs_g_[M_WBACK][M_WFINE] = reinterpret_cast<real_t**>(m_calloc(sizeof(real_t*) * n_ghost_det));
        //     for (sid_t ig = 0; ig < n_ghost_det; ig++) {
        //         gs_g_[M_WBACK][M_WCOAR][ig] = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * n_info_[M_WBACK][M_WCOAR]));
        //         gs_g_[M_WBACK][M_WFINE][ig] = reinterpret_cast<real_t*>(m_calloc(sizeof(real_t) * n_info_[M_WBACK][M_WFINE]));
        //     }
        // }

        // get the correct lifting filter given N
        if (N == 2 && Nt == 2) {
            // ha
            ha_[0] = -1.0 / 8.0;
            ha_[1] = +1.0 / 4.0;
            ha_[2] = +3.0 / 4.0;
            ha_[3] = +1.0 / 4.0;
            ha_[4] = -1.0 / 8.0;
            // ga
            ga_[0] = -1.0 / 2.0;
            ga_[1] = +1.0;
            ga_[2] = -1.0 / 2.0;
            // gs
            gs_[0] = +1.0 / 2.0;
            gs_[1] = +1.0 / 2.0;
            // modified filters - detail ghost 0
            // sid_t ig = 0;
            // front coarse
            // gs_g_[M_WFRONT][M_WCOAR][ig][0] = 1.0 / 2.0;
            // // front fine
            // gs_g_[M_WFRONT][M_WFINE][ig][0] = 3.0 / 7.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][1] = 1.0 / 7.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][2] = -1.0 / 14.0;
            // // back coarse
            // gs_g_[M_WBACK][M_WCOAR] = nullptr;
            // // back fine
            // gs_g_[M_WBACK][M_WFINE] = nullptr;
        } else if (N == 4 && Nt == 0) {
            // ha
            ha_[0] = 1.0;
            // ga
            ga_[0] = 1.0 / 16.0;
            ga_[1] = 0.0;
            ga_[2] = -9.0 / 16.0;
            ga_[3] = +1.0;
            ga_[4] = -9.0 / 16.0;
            ga_[5] = 0.0;
            ga_[6] = 1.0 / 16.0;
            // gs
            gs_[0] = -1.0 / 16.0;
            gs_[1] = 9.0 / 16.0;
            gs_[2] = 9.0 / 16.0;
            gs_[3] = -1.0 / 16.0;
            // // modified filters
            // sid_t ig = 0;
            // // front coarse
            // gs_g_[M_WFRONT][M_WCOAR][ig][0] = -1.0 / 16.0;
            // gs_g_[M_WFRONT][M_WCOAR][ig][1] = 9.0 / 16.0;
            // // front fine
            // gs_g_[M_WFRONT][M_WFINE][ig][0] = 9.0 / 16.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][1] = 0.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][2] = -1.0 / 16.0;
            // // back coarse
            // gs_g_[M_WBACK][M_WCOAR][ig][0] = 9.0 / 16.0;
            // gs_g_[M_WBACK][M_WCOAR][ig][1] = 9.0 / 16.0;
            // gs_g_[M_WBACK][M_WCOAR][ig][1] = -1.0 / 16.0;
            // // back fine
            // gs_g_[M_WBACK][M_WFINE][ig][0] = -1.0 / 16.0;
            // gs_g_[M_WBACK][M_WFINE][ig][1] = 0.0;
        } else if (N == 4 && Nt == 2) {
            // ha
            ha_[0] = +1.0 / 64.0;
            ha_[1] = 0.0;
            ha_[2] = -1.0 / 8.0;
            ha_[3] = +1.0 / 4.0;
            ha_[4] = +23.0 / 32.0;
            ha_[5] = +1.0 / 4.0;
            ha_[6] = -1.0 / 8.0;
            ha_[7] = 0.0;
            ha_[8] = +1.0 / 64.0;
            // ga
            ga_[0] = 1.0 / 16.0;
            ga_[1] = 0.0;
            ga_[2] = -9.0 / 16.0;
            ga_[3] = +1.0;
            ga_[4] = -9.0 / 16.0;
            ga_[5] = 0.0;
            ga_[6] = 1.0 / 16.0;
            // gs
            gs_[0] = -1.0 / 16.0;
            gs_[1] = 9.0 / 16.0;
            gs_[2] = 9.0 / 16.0;
            gs_[3] = -1.0 / 16.0;
            // // modified filters
            // // ghost 0
            // sid_t ig = 0;
            // // front coarse
            // gs_g_[M_WFRONT][M_WCOAR][ig][0] = -1.0 / 16.0;
            // gs_g_[M_WFRONT][M_WCOAR][ig][1] = 9.0 / 16.0;
            // gs_g_[M_WFRONT][M_WCOAR][ig][2] = 439.0 / 782.0;  // from python: 31617/56320
            // // front fine
            // gs_g_[M_WFRONT][M_WFINE][ig][0] = -369.0 / 7040.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][1] = -63.0 / 3520.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][2] = 44.0 / 4441.0;  // from python 279/28160
            // gs_g_[M_WFRONT][M_WFINE][ig][3] = 1.0 / 3520.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][4] = -9.0 / 7040.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][5] = 0.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][6] = 1.0 / 56320.0;
            // // ghost 1
            // ig = 1;
            // // front coarse
            // gs_g_[M_WFRONT][M_WCOAR][ig][0] = 0.0;
            // gs_g_[M_WFRONT][M_WCOAR][ig][1] = -1.0 / 16.0;
            // gs_g_[M_WFRONT][M_WCOAR][ig][2] = 503.0 / 880.0;
            // // front fine
            // gs_g_[M_WFRONT][M_WFINE][ig][0] = 211 / 440;
            // gs_g_[M_WFRONT][M_WFINE][ig][1] = 8 / 55;
            // gs_g_[M_WFRONT][M_WFINE][ig][2] = -59.0 / 440.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][3] = -1.0 / 55.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][4] = 17.0 / 880.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][5] = 0.0;
            // gs_g_[M_WFRONT][M_WFINE][ig][6] = -1.0 / 880.0;
            // // ghost 0
            // ig = 0;
            // // back coarse
            // gs_g_[M_WBACK][M_WCOAR][ig][0] = 73.0 / 128.0;
            // gs_g_[M_WBACK][M_WCOAR][ig][1] = 575.0 / 1024.0;
            // gs_g_[M_WBACK][M_WCOAR][ig][1] = -1.0 / 16.0;
            // // back fine
            // gs_g_[M_WBACK][M_WFINE][ig][0] = -1.0 / 1024.0;
            // gs_g_[M_WBACK][M_WFINE][ig][1] = 0.0;
            // gs_g_[M_WBACK][M_WFINE][ig][2] = 1.0 / 128.0;
            // gs_g_[M_WBACK][M_WFINE][ig][3] = -1.0 / 64.0;
            // gs_g_[M_WBACK][M_WFINE][ig][4] = -23.0 / 512.0;
            // gs_g_[M_WBACK][M_WFINE][ig][5] = -1.0 / 64.0;
        } else {
            m_assert(false, "wavelet N=%d.Nt=%d not implemented yet", N, Nt);
        }

        int nthread = omp_get_max_threads();
        temp_       = reinterpret_cast<real_t**>(m_calloc(nthread * sizeof(real_t*)));
        for (sid_t it = 0; it < nthread; it++) {
            temp_[it] = reinterpret_cast<real_t*>(m_calloc(len_gs_ * len_gs_ * len_gs_ * sizeof(real_t)));
        }

        m_log("Wavelet %d.%d with ga[%d], ha[%d], gs[%d], (%d,%d) ghosts needed", N, Nt, len_ga_, len_ha_, len_gs_, nghost_front(), nghost_back());
        //-------------------------------------------------------------------------
    }

    ~Wavelet() {
        m_begin;
        //-------------------------------------------------------------------------
        m_free(ha_);
        m_free(ga_);
        m_free(gs_);

        int nthread = omp_get_max_threads();
        for (sid_t it = 0; it < nthread; it++) {
            m_free(temp_[it]);
        }
        m_free(temp_);
        //-------------------------------------------------------------------------
        m_end;
    }

    /*
    * @name function overriding Interpolator class
    * @{
    */

   public:
    string Identity() const override { return "interpolating wavelet " + std::to_string(N) + "." + std::to_string(Nt); }

    // lid_t NGhostCoarseFront() const override { return n_info_[M_WFRONT][M_WCOAR]; }
    // lid_t NGhostCoarseBack() const override { return n_info_[M_WBACK][M_WCOAR]; }
    // front length
    lid_t ncoarsen_front() const override { return len_ha_ / 2; }
    lid_t nrefine_front() const override { return len_gs_ / 2 - 1; }
    lid_t ncriterion_front() const override { return len_ga_ / 2 - 1; }
    // back length
    lid_t ncoarsen_back() const override { return len_ha_ / 2 - 1; }
    lid_t nrefine_back() const override { return len_gs_ / 2; }
    lid_t ncriterion_back() const override { return len_ga_ / 2; }

    real_t Criterion(MemLayout* block, real_p data) override;

   protected:
    void Coarsen_(const interp_ctx_t* ctx) override;
    void Refine_(const interp_ctx_t* ctx) override;
    // void Copy_(const sid_t dlvl, const interp_ctx_t* ctx) override;
    /* @}*/

    /*
    * @name implementation specific function
    * @{
    */
   public:
    void Details(MemLayout* block, real_p data, real_t* criterion);

   protected:
    void Detail_(const interp_ctx_t* ctx, real_t* details_inf_norm);

    /* @}*/

    // void RefineGhost_(const interp_ctx_t* ctx)  override;
};

#endif  // SRC_WAVELET_HPP_

#include "wavelet.inc"