
#ifndef SRC_INTERPOLATE_HPP_
#define SRC_INTERPOLATE_HPP_

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

    /**
     * @name position pointers
     * 
     * They both refer the position (0,0,0) of the target, hence the ghostsize is assumed to be zero
     * @{
     */
    real_p sdata;  //!< refers the (0,0,0) location of the target memory, in the source layout
    real_p tdata;  //!< refers the (0,0,0) location of the target memory
    /** @} */
} interp_ctx_t;


/**
 * @brief defines a set of function used to interpolate
 * 
 * The memory description relies on the MemLayout object and they are working on one dimension at a time.
 * 
 */
class Interpolator {
   public:
    virtual void Criterion(MemLayout* block, real_p data, real_t* criterion) = 0;

    virtual void Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);

   protected:
    virtual void Coarsen_(const interp_ctx_t* ctx, const lid_t dlvl) const = 0;
    virtual void Refine_(const interp_ctx_t* ctx) const                    = 0;
    virtual void Copy_(const interp_ctx_t* ctx) const                      = 0;
};

#endif  // SRC_INTERPOLATE_HPP_
