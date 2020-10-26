
#ifndef SRC_INTERPOLATE_HPP_
#define SRC_INTERPOLATE_HPP_

#include <string>

#include "field.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"
#include "p8est.h"

// //................................................
// // WAVELET 2.2
// const static lid_t  len_ha_2_2 = 5;
// const static lid_t  len_gs_2_2 = 2;
// const static real_t ha_2_2[5]  = {-0.125, 0.25, 0.75, 0.25, -0.125};
// const static real_t gs_2_2[2]  = {0.5, 0.5};
// //................................................
// // WAVELET 4.0
// const static lid_t  len_ha_4_0 = 1;
// const static lid_t  len_gs_4_0 = 4;
// const static real_t ha_4_0[1]  = {1.0};
// const static real_t gs_4_0[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};
// //................................................
// // WAVELET 4.2
// const static lid_t  len_ha_4_0 = 9;
// const static lid_t  len_gs_4_0 = 4;
// const static real_t ha_4_0[9]  = {0.015625, 0.0, -0.125, 0.25, 0.71875, 0.25, -0.125, 0.0, 0.015625};
// const static real_t gs_4_0[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};
// //................................................
// // WAVELET 4.4
// const static lid_t  len_ha_4_4 = 13;
// const static lid_t  len_gs_4_4 = 4;
// const static real_t ha_4_4[13] = {-1.0 / 512.0, 0.0, 9.0 / 256.0, -1.0 / 32.0, -63.0 / 512.0, 9.0 / 32.0, 87.0 / 128.0, 9.0 / 32.0, -63.0 / 512.0, -1.0 / 32.0, 9.0 / 256.0, 0.0, -1.0 / 512.0};
// const static real_t gs_4_4[4]  = {-0.0625, 0.5625, 0.5625, -0.0625};

// //................................................
// // WAVELET 6.0


// //................................................
// // WAVELET 6.2


// //................................................
// WAVELET 6.4

/**
 * @brief Defines all the required information to perform an interpolation on a given block
 * 
 * Since the interpolation is done within threads, those values cannot belong to the object 
 * and must be created each time an interpolation is needed
 */
typedef struct interp_ctx_t {
    lid_t srcstr;       //!< the source stride
    lid_t trgstr;       //!< the target stride
    lid_t trgstart[3];  //!< first index needed in the target memory
    lid_t trgend[3];    //!< last index needed in the target memory
#ifndef NDEBUG
    // for debug only
    lid_t srcstart[3];  //!< first index available in the source memory
    lid_t srcend[3];    //!< last index available in the source memory
#endif
    real_t alpha = 0.0;  //!< the constant multiplication factor: target = alpha * constant + interpolation(source)

    /**
     * @name position pointers
     * 
     * They both refer the position (0,0,0) of the target, hence the ghostsize is assumed to be zero
     * @{
     */
    data_ptr sdata;  //!< refers the (0,0,0) location of the target memory, in the source memory layout
    data_ptr cdata;  //!< refers the (0,0,0) location of the target memory, in the constant memory layout
    data_ptr tdata;  //!< refers the (0,0,0) location of the target memory
    /** @} */
} interp_ctx_t;

/**
 * @brief defines the most basic Interpolating wavelet class
 * 
 * The target field is computed as alpha * constant field + interpolation(source field).
 * The interpolation procedure is one of the following:
 * - a copy
 * - a get/put RMA operation on an already activated window
 * - a refinement (to be provided by the child class)
 * - a coarsening (to be provided by the child class)
 * 
 */
// template< lda_t N, lda_t NT>
class InterpolatingWavelet {
   public:
    //................................................
    // need for empty constructor/destructor to call the virtual ones
    explicit InterpolatingWavelet(){};
    virtual ~InterpolatingWavelet(){};

    //................................................
    /**
    * @name basic implemented interpolating functions
    * @{
    */
   public:
    virtual void Copy(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg) const;
    virtual void Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg) const;
    virtual void Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg, const real_t alpha, data_ptr data_cst) const;
    virtual void GetRma(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, MPI_Aint disp_src, const MemLayout* block_trg, data_ptr data_trg, rank_t src_rank, MPI_Win win) const;
    virtual void PutRma(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr ptr_src, const MemLayout* block_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win) const;

    virtual real_t Criterion(MemLayout* block, data_ptr data, mem_ptr data_tmp) const;
    void           Details(MemLayout* block, data_ptr data_block, mem_ptr data_tmp, real_t* details_max) const;

    //................................................
   public:
    virtual const sid_t N() const  = 0;
    virtual const sid_t Nt() const = 0;
    std::string         Identity() const { return "interpolating wavelet " + std::to_string(N()) + "." + std::to_string(Nt()); }

    /** @} */

    //................................................
    /**
     * @name filter management and various lengths
     * @{
     */
   public:
    virtual const sid_t   len_ha() const = 0;
    virtual const sid_t   len_ga() const = 0;
    virtual const sid_t   len_gs() const = 0;
    virtual const real_t* ha() const     = 0;  //!< ha = compute the scaling from a given level -> coarsening
    virtual const real_t* ga() const     = 0;  //!< ga = compute the details from a given level -> detail computation: dual lifting only :-)
    virtual const real_t* gs() const     = 0;  //!< gs = reconstruction from the scaling and 0 detail of the scaling coef -> refinement: should be -ga!!

    // shift for the details
    const lid_t shift_front() const { return m_max((2 * Nt() - 1) / 2, 0); };     //!< return the num of detail to take into account outside the block, in front
    const lid_t shift_back() const { return m_max((2 * Nt() - 1) / 2 - 1, 0); };  //!< return the num of detail to take into account outside the block, in the back

    // nghosts
    const lid_t ncoarsen_front() const { return (len_ha() / 2); };                 //!< returns the number of gp needed for the coarsening operation, in front
    const lid_t ncoarsen_back() const { return (len_ha() / 2 - 1); };              //!< returns the number of gp needed for the coarsening operation, in the back
    const lid_t nrefine_front() const { return (len_gs() / 2 - 1); };              //!< returns the number of gp needed for the refinement operation, in front
    const lid_t nrefine_back() const { return (len_gs() / 2); };                   //!< returns the number of gp needed for the refinement operation, in the back
    const lid_t ncriterion_front() const { return shift_front() + len_gs() - 1; }  //!< returns the number of gp needed for the detail operation, in front
    const lid_t ncriterion_back() const { return shift_back() + len_gs() - 1; }    //!< returns the number of gp needed for the detail operation, in the back

    // half limits
    const sid_t ha_half_lim() const { return (len_ha() / 2); };
    const sid_t ga_half_lim() const { return (len_ga() / 2); };
    const sid_t gs_half_lim() const { return (len_gs() / 2) - 1; };

    // filters:
    const real_t* ha_filter() const { return ha() + ha_half_lim(); };
    const real_t* ga_filter() const { return ga() + ga_half_lim(); };
    const real_t* gs_filter() const { return gs() + gs_half_lim(); };

    // nghosts
    const lid_t nghost_front() const { return m_max(ncoarsen_front(), m_max(ncriterion_front(), nrefine_front())); }
    const lid_t nghost_back() const { return m_max(ncoarsen_back(), m_max(ncriterion_back(), nrefine_back())); }

    /** @} */

    inline lid_t CoarseNGhostFront() const {
        // we need the max between the number of scaling in front
        //      = (nghost_front() / 2)
        // and the number of details
        //      = (nghost_front() / 2) + nrefine_front()
        const lid_t gp = (nghost_front() / 2) + nrefine_front();
        return gp;
    };
    inline lid_t CoarseNGhostBack() const {
        // we need the max between the number of scaling
        //      = (interp->nghost_back()+1) / 2
        const lid_t gp_scaling = ((nghost_back() + 1) / 2);
        // and the number of details
        //      = (interp->nghost_back()) / 2 +  + interp->nrefine_back()
        const lid_t gp_detail = ((nghost_back()) / 2) + nrefine_back();
        return m_max(gp_scaling, gp_detail);
    };
    inline size_t CoarseStride() const {
        const lid_t gp_front = CoarseNGhostFront();
        const lid_t gp_back  = CoarseNGhostBack();
        return gp_front + M_HN + gp_back;
    };
    inline size_t CoarseMemSize() const {
        return CoarseStride() * CoarseStride() * CoarseStride() * sizeof(real_t);
    };
    /**
    * @brief given a starting a ghost range for a block, transform it into the needed range for its coarse representation
    * if the point is inside the block, we divide its id by 2, if the point is in the ghost points,
    * we take ALL the coarse ghost points
    * 
    * @warning this is magic...
    * 
    * we compute:
    * b = a + M_N, such that
    *      if a is in the ghost points, b < M_N
    *      if a is in the center points, M_N <= b < 2* M_N
    *      if a is in the ghost points, 2* M_N <= b
    * 
    * c is 0,1,2,3:
    *      = 0 if a is in the negative ghost points
    *      = 1 if a is in the center points, including 0, we scale it by two but preserve the odd numbers
    *      = 2 if a is M_N
    *      = 3 if a is in the negative GP
    * 
    * the correct ID is returned based on the value of c
    * 
    * @param a the id in the block
    * @param interp the interpolator used
    * @return lid_t 
    */
    inline lid_t CoarseFromBlock(const lid_t a) const {
        const lid_t gp_front = CoarseNGhostFront();
        const lid_t gp_back  = CoarseNGhostBack();
        const lid_t b        = (a + M_N);
        const lid_t c        = (b / M_N) + (a > M_N);
        const lid_t res[4]   = {-gp_front, (a / 2) + (a % 2), M_HN, M_HN + gp_back};
        // return the correct choice
        return res[c];
    }
    inline void CoarseFromFine(const SubBlock* block_fine, SubBlock* block_coarse) const {
        lid_t coarse_start[3], coarse_end[3];
        for (lda_t id = 0; id < 3; ++id) {
            coarse_start[id] = CoarseFromBlock(block_fine->start(id));
            coarse_end[id]   = CoarseFromBlock(block_fine->end(id));
        }
        block_coarse->Reset(CoarseNGhostFront(), CoarseStride(), coarse_start, coarse_end);
    }
    inline void CoarseFromFine(const lid_t id_fine[3], lid_t* id_coarse) const {
        for (lda_t id = 0; id < 3; ++id) {
            id_coarse[id] = CoarseFromBlock(id_fine[id]);
        }
    }

   protected:
    // call Copy_, Coarsen_ or Refine_
    virtual void DoMagic_(const level_t dlvl, const bool force_copy, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg, const real_t alpha, const data_ptr data_cst) const;

    /**
     * @name Interpolation functions, to be implemented
     * @{
     */
    virtual void Coarsen_(const interp_ctx_t* ctx) const;
    virtual void Refine_(const interp_ctx_t* ctx) const;
    virtual void Detail_(const interp_ctx_t* ctx, real_t* details_max) const;
    /** @} */

    // defined function -- might be overriden
    virtual void Copy_(const level_t dlvl, const interp_ctx_t* ctx) const;
};




#endif  // SRC_INTERPOLATE_HPP_
