
#ifndef SRC_INTERPOLATE_HPP_
#define SRC_INTERPOLATE_HPP_

#include <string>

#include "field.hpp"
#include "ghostblock.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"
#include "p8est.h"

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
    virtual void Copy(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg);
    virtual void Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg);
    virtual void Interpolate(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr data_src, const MemLayout* block_trg, data_ptr data_trg, const real_t alpha, data_ptr data_cst);
    virtual void GetRma(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, MPI_Aint disp_src, const MemLayout* block_trg, data_ptr data_trg, rank_t src_rank, MPI_Win win);
    virtual void PutRma(const level_t dlvl, const lid_t shift[3], const MemLayout* block_src, const data_ptr ptr_src, const MemLayout* block_trg, MPI_Aint disp_trg, rank_t trg_rank, MPI_Win win);

    virtual real_t Criterion(MemLayout* block, real_p data) const;
    void           Details(MemLayout* block, real_p data, real_t* details_max) const;

    //................................................
   public:
    virtual const sid_t N() const  = 0;
    virtual const sid_t Nt() const = 0;
    string              Identity() const { return "interpolating wavelet " + std::to_string(N()) + "." + std::to_string(Nt()); }

    /** @} */

    //................................................
    /**
     * @name filter management and ghost lengths
     * @{
     */
   public:
    virtual const sid_t   len_ha() const = 0;
    virtual const sid_t   len_ga() const = 0;
    virtual const sid_t   len_gs() const = 0;
    virtual const real_t* ha() const     = 0;  //!< ha = compute the scaling from a given level -> coarsening
    virtual const real_t* ga() const     = 0;  //!< ga = compute the details from a given level -> detail computation: dual lifting only :-)
    virtual const real_t* gs() const     = 0;  //!< gs = reconstruction from the scaling and 0 detail of the scaling coef -> refinement: should be -ga!!

    // nghosts
    const lid_t ncoarsen_front() const { return (len_ha() / 2); };
    const lid_t ncoarsen_back() const { return (len_ha() / 2 - 1); };
    const lid_t nrefine_front() const { return (len_gs() / 2 - 1); };
    const lid_t nrefine_back() const { return (len_gs() / 2); };
    const lid_t ncriterion_front() const { return (len_ga() / 2 - 1); };
    const lid_t ncriterion_back() const { return (len_ga() / 2); };

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
