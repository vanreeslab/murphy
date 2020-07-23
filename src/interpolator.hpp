
#ifndef SRC_INTERPOLATE_HPP_
#define SRC_INTERPOLATE_HPP_

#include <string>

#include "field.hpp"
#include "ghostblock.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"
#include "p8est.h"

// /**
//  * @brief define what a standard interpolation function looks like
//  */
// typedef void (*InterpFunction)(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);

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
    // sid_t* normal = nullptr;  //!< normal to the ghost zone
    /**
     * @name position pointers
     * 
     * They both refer the position (0,0,0) of the target, hence the ghostsize is assumed to be zero
     * @{
     */
    real_p sdata;  //!< refers the (0,0,0) location of the target memory, in the source memory layout
    real_p cdata;  //!< refers the (0,0,0) location of the target memory, in the constant memory layout
    real_p tdata;  //!< refers the (0,0,0) location of the target memory
    /** @} */
} interp_ctx_t;

using std::string;

/**
 * @brief defines a set of function used to interpolate
 * 
 * The memory description relies on the MemLayout object and they are working on one dimension at a time.
 * 
 */
class Interpolator {
   public:
    // need to define the destructor as virtual to be sure to pass by by the wavelet destructor
    virtual ~Interpolator(){};

    /**
    * @name accessible interpolating functions
    * @{
    */
    virtual void Copy(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);
    virtual void Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);
    virtual void Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg, const real_t alpha, real_p data_cst);
    // the RMA variant of Copy
    virtual void GetRma(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, MPI_Aint disp_src, MemLayout* block_trg, real_t* data_trg, int src_rank, MPI_Win win);
    /** @} */

    virtual real_t Criterion(MemLayout* block, real_p data) = 0;
    virtual string Identity() const                         = 0;

    /**
     * @name filter length - to be implemented
     * @{
     */
    virtual lid_t ncoarsen_front() const   = 0;
    virtual lid_t ncriterion_front() const = 0;
    virtual lid_t nrefine_front() const    = 0;
    virtual lid_t ncoarsen_back() const    = 0;
    virtual lid_t ncriterion_back() const  = 0;
    virtual lid_t nrefine_back() const     = 0;
    /** @} */

    /**
     * @name ghost length, worst case of everyfilter
     * @{
     */
    virtual lid_t nghost_front() const { return m_max(ncoarsen_front(), m_max(ncriterion_front(), nrefine_front())); }
    virtual lid_t nghost_back() const { return m_max(ncoarsen_back(), m_max(ncriterion_back(), nrefine_back())); }
    /** @} */

   protected:
    // call Copy_, Coarsen_ or Refine_
    virtual void DoMagic_(const sid_t dlvl, const bool force_copy, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg, const real_t alpha, real_p data_cst);

    // to be implemented functions
    virtual void Coarsen_(const interp_ctx_t* ctx) = 0;
    virtual void Refine_(const interp_ctx_t* ctx)  = 0;

    // defined function -- might be overriden
    virtual void Copy_(const sid_t dlvl, const interp_ctx_t* ctx);
};

#endif  // SRC_INTERPOLATE_HPP_
