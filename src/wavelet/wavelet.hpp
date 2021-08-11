#ifndef SRC_WAVELET_HPP_
#define SRC_WAVELET_HPP_

#include <string>

#include "core/macros.hpp"
#include "core/memlayout.hpp"
#include "core/pointers.hpp"
#include "core/types.hpp"
#include "p8est.h"
#include "subblock.hpp"

/**
 * @brief Defines all the required information to perform an interpolation on a given block
 * 
 * Since the interpolation is done within threads, those values cannot belong to the object 
 * and must be created each time an interpolation is needed
 * 
 * @warning the exact meaning of each field depends on the wavelet function called! please check the documentation.
 */
typedef struct interp_ctx_t {
    bidx_t trgstart[3];  //!< first index needed in the target memory
    bidx_t trgend[3];    //!< last index needed in the target memory
#ifndef NDEBUG
    // for debug only
    bidx_t srcstart[3];  //!< first index available in the source memory
    bidx_t srcend[3];    //!< last index available in the source memory
#endif

    real_t alpha = 0.0;  //!< a factor used in different ways depending on the context

    /**
     * @name position pointers
     * @{
     */
    const_data_ptr sdata;  //!< refers the (0,0,0) location of the source memory, in the source memory layout
    data_ptr       tdata;  //!< refers the (0,0,0) location of the target memory
    // data_ptr       temp;   //!< refers the (0,0,0) location of the temporary memory (has the layour of trg!)
    /** @} */

    /** @brief Default, parameterless constructor of the structure */
    interp_ctx_t(){};

#ifndef NDEBUG
    /** @brief Constructor of the structure to avoid the default initialisation of the const_data_ptr and of the data_ptr */
    interp_ctx_t(const_data_ptr sdata_in, bidx_t srcstr_in, bidx_t srcstart_in[3],bidx_t srcend_in[3], data_ptr tdata_in, bidx_t trgstr_in, bidx_t trgstart_in[3],bidx_t trgend_in[3], real_t alpha_in):sdata(sdata_in, srcstr_in), tdata(tdata_in, trgstr_in){
        for(bidx_t i = 0; i < 3 ; i++){
            trgstart[i] = trgstart_in[i];
            trgend[i]   = trgend_in[i];

            srcstart[i] = srcstart_in[i];
            srcend[i]   = srcend_in[i];
        }
        alpha = alpha_in;
    };
#else 
/** @brief Constructor of the structure to avoid the default initialisation of the const_data_ptr and of the data_ptr */
    interp_ctx_t(const_data_ptr sdata_in, bidx_t srcstr_in, data_ptr tdata_in, bidx_t trgstr_in, bidx_t trgstart_in[3],bidx_t trgend_in[3], real_t alpha_in):sdata(sdata_in, srcstr_in), tdata(tdata_in, trgstr_in){
        for(bidx_t i = 0; i < 3 ; i++){
            trgstart[i] = trgstart_in[i];
            trgend[i]   = trgend_in[i];
        }
        alpha = alpha_in;
    };
#endif


} interp_ctx_t;

#define M_GS_MIN 0

/**
 * @brief defines a generic wavelet class
 * 
 * The target field is computed as alpha * constant field + interpolation(source field).
 * The interpolation procedure is one of the following:
 * - a copy
 * - a get/put RMA operation on an already activated window
 * - a refinement (to be provided by the child class)
 * - a coarsening (to be provided by the child class)
 * 
 */
class Wavelet {
   public:
    //................................................
    // need for empty constructor/destructor to call the virtual ones
    // explicit Wavelet(){};
    virtual ~Wavelet(){};

    //................................................
    /**
    * @name basic implemented interpolating functions
    * @{
    */
   public:
    void Copy(const level_t dlvl, const lid_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg) const;
    void Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg) const;
    void Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg, const real_t alpha, const_data_ptr data_cst) const;
    void GetRma(const level_t dlvl, const lid_t shift[3], const MemLayout*  block_src, MPI_Aint disp_src, const MemLayout*  block_trg, data_ptr data_trg, rank_t src_rank, MPI_Win win) const;
    void PutRma(const level_t dlvl, const lid_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win) const;

    real_t Criterion(/* source */ const MemLayout* const  block_src, const const_data_ptr& data,
                     /* target */ const MemLayout* const  detail_block) const;
    void   Details(/* source */ const MemLayout* const  block_src, const const_data_ptr& data,
                 /* target */ const MemLayout* const  detail_block, const data_ptr& detail, const real_t tol,
                 /* output*/ real_t*  details_max) const;
    void   SmoothOnMask(/* source */ const MemLayout* const  block_src,
                      /* target */ const MemLayout* const  block_trg, const data_ptr& data,
                      /* detail */ const MemLayout* const  detail_block, const data_ptr& detail_mask) const;
    void   OverwriteDetails(/* source */ const MemLayout* const  block_src,
                          /* target */ const MemLayout* const  block_trg, const data_ptr& data) const;
    void   StoreDetails(/* source */ const MemLayout* const  block_src, const const_data_ptr& data,
                           /* target */ const MemLayout* const  block_detail, const data_ptr& detail) const;
    /** @} */

    //................................................
   public:
    /**
    * @name Identity functions
    * @{
    */
    std::string Identity() const { return "interpolating wavelet " + std::to_string(N()) + "." + std::to_string(Nt()); }

    virtual const short_t N() const  = 0;  // return the interpolation order
    virtual const short_t Nt() const = 0;  // return the moments order
    virtual const real_t eps_const() const = 0; //return the theoretical factor in front of the epsilon

    const bool smoothed() const { return (Nt() != 0); };  // return if the wavelet details will smooth or not.

    /** @} */

    //................................................
    /**
     * @name filter management and associated lengths
     * @{
     */
   public:
    virtual const short_t len_ha() const = 0;  //!< length of filter Ha
    virtual const short_t len_ga() const = 0;  //!< length of filter Ga
    virtual const short_t len_js() const = 0;  //!< length of filter Js
    virtual const short_t len_ks() const = 0;  //!< length of filter Ks

    /**
     * @brief returns the number of gp needed for the coarsening operation = get the scaling, in front
     */
    const bidx_t nghost_front_coarsen() const { return m_max(len_ha() / 2, 0); };
    /**
     * @brief returns the number of gp needed for the coarsening operation = get the scaling, in the back
     */
    const bidx_t nghost_back_coarsen() const { return m_max(len_ha() / 2 - 1, 0); };
    /**
     * @brief returns the number of gp needed for the refinement operation
     * 
     * it overestimates the actuall needed number, but it's not the actual contrain
     */
    const bidx_t nghost_front_refine() const {
        const bidx_t n_js   = len_js() / 2;
        const bidx_t n_ks   = len_ks() / 2;
        const bidx_t n_scal = n_js - (n_js % 2);        // remove the last point if it's a detail
        const bidx_t n_det  = n_ks - (1 - (n_ks % 2));  // remove the last point if it's a detail
        // const bidx_t last_pt = m_max(n_scal, n_det - 1);
        const bidx_t last_pt = m_max(n_scal, n_det);
        // const bidx_t last_scal = last_pt - (last_pt % 2);  // remove the last point if it's a detail
        // this is how the access is actually made as we skip every detail point
        return m_max((last_pt) / 2, 0);
    };
    /**
     * @brief returns the number of gp needed for the refinement operation assuming detail = 0, in the back
     * 
     * it overestimates the actuall needed number, but it's not the actual contrain
     */
    const bidx_t nghost_back_refine() const {
        const bidx_t n_js   = len_js() / 2;
        const bidx_t n_ks   = len_ks() / 2;
        const bidx_t n_scal = n_js - (n_js % 2);        // remove the last point if it's a detail
        const bidx_t n_det  = n_ks - (1 - (n_ks % 2));  // remove the last point if it's a detail
        // const bidx_t last_pt = m_max(n_scal - 1, n_det);
        const bidx_t last_pt = m_max(n_scal, n_det);
        // const bidx_t last_scal = last_pt - (1 - (last_pt % 2));  // remove the last point if it's a detail
        // this is how the access is actually made as we skip every detail point
        return m_max((last_pt + 1) / 2, 0);
    };

    const bidx_t nghost_front_overwrite() const {
        const bidx_t n_ga   = len_ga() / 2;
        // overwrites only takes place at the first point
        const bidx_t last_pt = n_ga-1;
        return m_max(last_pt, 0);
    };
    const bidx_t nghost_back_overwrite() const {
        const bidx_t n_ga = len_ga() / 2;
        // overwrites only takes place at the first point
        const bidx_t last_pt = n_ga;
        return m_max(last_pt, 0);
    };

    /**
     * @brief number of details to compute in front of the block while doing the criterion
     * 
     * To guarantee the consitency in the ghosting and refinement in case of jump in resolution
     * we need no influence of any detail on the current block
     * 
     * This is the number of detail OUTSIDE the block that will influence my last point!
     * 
     */
    const bidx_t ndetail_citerion_extend_front() const {
        const bidx_t n_js   = len_js() / 2;
        const bidx_t n_ks   = len_ks() / 2;
        const bidx_t n_scal = m_max(n_js, 0);      //- (1 - (n_js % 2));  // remove the last point if it's a scaling
        const bidx_t n_det  = m_max(n_ks - 1, 0);  //- (n_ks % 2);        // remove the last point if it's a scaling
        // return m_sign(Nt()) * m_max(n_scal, n_det - 1);
        return m_max(n_scal, n_det);
    };

    /**
     * @brief number of details to compute at the back of the block while doing the criterion
     * 
     * To guarantee the consitency in the ghosting and refinement in case of jump in resolution
     * we need no influence of any detail on the current block
     * 
     * This is the number of detail OUTSIDE the block that will influence my last point!
     * 
     */
    const bidx_t ndetail_citerion_extend_back() const {
        const bidx_t n_js   = len_js() / 2;
        const bidx_t n_ks   = len_ks() / 2;
        const bidx_t n_scal = m_max(n_js - 1, 0);  //- (1 - (n_js % 2));  // remove the last point if it's a scaling
        const bidx_t n_det  = m_max(n_ks, 0);      //- (n_ks % 2);        // remove the last point if it's a scaling
        // return m_sign(Nt()) * m_max(n_scal-1, n_det);
        // return m_max(0, m_max(n_scal - 1, n_det));
        return m_max(n_scal, n_det);
    };

    /**
     * @brief number of details to take into account in front of the block while smoothing over a jump of resolution
     */
    const bidx_t ndetail_smooth_extend_front() const {
        // const bidx_t n_js   = len_js() / 2;
        // const bidx_t n_ks   = len_ks() / 2;
        // const bidx_t n_scal = n_js - (1 - (n_js % 2));  // remove the last point if it's a scaling
        // const bidx_t n_det  = n_ks - (n_ks % 2);        // remove the last point if it's a scaling
        // // return m_sign(Nt()) * m_max(n_scal, n_det - 1);
        // return m_max(0, m_max(n_scal, n_det - 1));
        // this is the same number as for the criterion!
        return ndetail_citerion_extend_front();
    };

    /**
     * @brief number of details to take into account at the back of the block while smoothing over a jump of resolution
     */
    const bidx_t ndetail_smooth_extend_back() const {
        // const bidx_t n_js   = len_js() / 2;
        // const bidx_t n_ks   = len_ks() / 2;
        // const bidx_t n_scal = n_js - (1 - (n_js % 2));  // remove the last point if it's a scaling
        // const bidx_t n_det  = n_ks - (n_ks % 2);        // remove the last point if it's a scaling
        // // return m_sign(Nt()) * m_max(n_scal - 1, n_det);
        // return m_max(0, m_max(n_scal - 1, n_det));
        // this is the same number as for the criterion!
        return ndetail_citerion_extend_back();
    };

    /**
     * @brief returns the number of gp needed to apply the smooth over a resolution jump
     */
    const lid_t nghost_front_criterion_smooth() const {
        // // this is the last detail I will ever need, already takes into account that we are in the front
        // const bidx_t max_detail = m_max(ndetail_citerion_extend_front(), ndetail_smooth_extend_front());
        // const bidx_t offset     = 1 * (max_detail == 0);  // we must remove 1 if max_detail ==0
        // return m_max(max_detail + (len_ga() / 2) - offset, 0);
        // get the last point, might be a detail or a scaling
        const bidx_t last_pt     = m_max(ndetail_citerion_extend_front(), ndetail_smooth_extend_front());
        const bidx_t lim_scaling = m_max(last_pt + (len_ha() / 2), last_pt - 1 + (len_ga() / 2));
        const bidx_t lim_detail  = m_max(last_pt + (len_ga() / 2), last_pt - 1 + (len_ha() / 2));
        const bool   is_detail   = (last_pt % 2) == 1;
        return (is_detail)*lim_detail + (!is_detail) * lim_scaling;
    };
    /**
     * @brief returns the number of gp needed to apply the smooth over a resolution jump
     */
    const lid_t nghost_back_criterion_smooth() const {
        const bidx_t last_pt = m_max(ndetail_citerion_extend_back(), ndetail_smooth_extend_back());
        // compute the max between computing the last point and the point just before
        const bidx_t lim_scaling = m_max(last_pt + (len_ha() / 2), last_pt - 1 + (len_ga() / 2));
        const bidx_t lim_detail  = m_max(last_pt + (len_ga() / 2), last_pt - 1 + (len_ha() / 2));
        const bool   is_detail   = (last_pt % 2) == 0;
        return (is_detail)*lim_detail + (!is_detail) * lim_scaling;
    };


    // nghosts needed for the interpolation
    const bidx_t nghost_front() const {
        const bidx_t min_wavelet = m_max(nghost_front_coarsen(),
                                         m_max(nghost_front_refine(),
                                               nghost_front_criterion_smooth()));
        return m_max(min_wavelet, M_GS_MIN);
    }
    const bidx_t nghost_back() const {
        const bidx_t min_wavelet = m_max(nghost_back_coarsen(),
                                         m_max(nghost_back_refine(),
                                               nghost_back_criterion_smooth()));
        return m_max(min_wavelet, M_GS_MIN);
    }
    
    /** @} */

    /**
     * @name Coarse sizes
     * @{
     */
    inline bidx_t CoarseNGhostFront(const bidx_t nghost_front) const {
        // we need the max between the number of scaling in front
        //      = (nghost_front / 2)
        // and the number of details that need to be computed
        //      = (nghost_front / 2) + nrefine_front()
        const bidx_t gp = ((nghost_front + 1) / 2) + nghost_front_refine();
        return gp;
    };
    inline bidx_t CoarseNGhostBack(const bidx_t nghost_back) const {
        // we need the max between the number of scaling
        //      = (interp->nghost_back+1) / 2
        const bidx_t gp_scaling = ((nghost_back + 1) / 2);
        // and the number of details
        //      = (interp->nghost_back()) / 2 +  + interp->nrefine_back()
        const bidx_t gp_detail = ((nghost_back) / 2) + nghost_back_refine();
        return m_max(gp_scaling, gp_detail);
    };
    inline bidx_t CoarseStride(const bidx_t ghost_size[2]) const {
        const bidx_t gp_front = CoarseNGhostFront(ghost_size[0]);
        const bidx_t gp_back  = CoarseNGhostBack(ghost_size[1]);
        return gp_front + M_NHALF + gp_back;
    };
    inline bidx_t CoarseSize(const bidx_t ghost_size[2]) const {
        const bidx_t stride = CoarseStride(ghost_size);
        return stride * stride * stride;
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
    *      = 3 if a is in the positive GP
    * 
    * the correct ID is returned based on the value of c
    * 
    * @param a the id in the block
    * @param interp the Wavelet used
    * @return lid_t 
    */
    // inline bidx_t CoarseFromBlock(const lid_t a) const {
    //     const bidx_t gp_front = CoarseNGhostFront();
    //     const bidx_t gp_back  = CoarseNGhostBack();
    //     const bidx_t b        = (a + M_N);
    //     const bidx_t c        = (b / M_N) + (a > M_N);
    //     const bidx_t res[4]   = {-gp_front, (a / 2) + (a % 2), M_NHALF, M_NHALF + gp_back};
    //     // return the correct choice
    //     return res[c];
    // }
    // inline void CoarseFromFine(const MemLayout* block_fine, SubBlock* block_coarse) const {
    //     bidx_t coarse_start[3], coarse_end[3];
    //     for (bidx_t id = 0; id < 3; ++id) {
    //         coarse_start[id] = CoarseFromBlock(block_fine->start(id));
    //         coarse_end[id]   = CoarseFromBlock(block_fine->end(id));
    //     }
    //     block_coarse->Reset(CoarseNGhostFront(), CoarseStride(), coarse_start, coarse_end);
    // }
    // inline void CoarseFromFine(const bidx_t id_fine[3], bidx_t* id_coarse) const {
    //     for (bidx_t id = 0; id < 3; ++id) {
    //         id_coarse[id] = CoarseFromBlock(id_fine[id]);
    //     }
    // }
    /** @} */

    //-------------------------------------------------------------------------
   protected:
    /**
     * @name Interpolation functions, to be implemented
     * @{
     */
    virtual void DoMagic_(const level_t dlvl, const bool force_copy, const lid_t shift[3], const MemLayout*  block_src, const_data_ptr data_src, const MemLayout*  block_trg, data_ptr data_trg, const real_t alpha, const_data_ptr data_cst) const;

    virtual void Coarsen_(const interp_ctx_t* const  ctx) const = 0;
    virtual void RefineZeroDetails_(const interp_ctx_t* const  ctx) const  = 0;
    virtual void OverwriteDetailsDualLifting_(const interp_ctx_t* const  ctx) const  = 0;
    // virtual void Scaling_(const interp_ctx_t* const  ctx) const                                                   = 0;
    virtual void Detail_(const interp_ctx_t* const  ctx, real_t* const  details_max) const                  = 0;
    virtual void ForwardWaveletTransform_(const interp_ctx_t* const  ctx, real_t* const  details_max) const = 0;
    virtual void Smooth_(const interp_ctx_t* const  ctx) const                                                    = 0;
    virtual void Clear_(const interp_ctx_t* const  ctx) const                                                     = 0;
    // virtual void WriteDetail_(const interp_ctx_t*  ctx) const                                     = 0;
    /** @} */

    // defined function -- might be overriden
    virtual void Copy_(const level_t dlvl, const interp_ctx_t*  ctx) const;
};

#endif  // SRC_WAVELET_HPP_
