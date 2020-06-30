
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
    real_t alpha  = 0.0;      //!< the constant multiplication factor: target = alpha * constant + interpolation(source)
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
    /**
    * @brief returns a positive, single value that represents a refinement/coarsening criterion
    * 
    * @param block the memory layout to use to compute that value
    * @param data the memory adress refering to point (0,0,0)
    * @return real_t the criterion value, always >= 0
    */
    virtual real_t Criterion(MemLayout* block, real_p data) = 0;

    // different interpolating functions
    virtual void Copy(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);
    virtual void Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);
    virtual void Interpolate(const sid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg, const real_t alpha, real_p data_cst);

    virtual void DoMagic(const sid_t dlvl, const bool force_copy, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg, const real_t alpha, real_p data_cst);

    /**
    * @brief returns a string identifying the operator
    */
    virtual string Identity() const = 0;
    /**
    * @brief returns how many coarse points are needed for the computation of a fine ghost values (front and back)
    */
    virtual lid_t NGhostCoarseFront() const = 0;
    virtual lid_t NGhostCoarseBack() const  = 0;
    /**
    * @brief returns how many ghost points are needed for a block (front and back)
    */
    virtual lid_t NGhostFineFront() const = 0;
    virtual lid_t NGhostFineBack() const  = 0;

   protected:
    /**
    * @brief coarsen the information for blocks having a jump in level
    * 
    * @param ctx the interpolation context, see @ref interp_ctx_t
    */
    virtual void Coarsen_(const interp_ctx_t* ctx) = 0;
    /**
     * @brief refines by one level the information
     * 
     * @param ctx the interpolation context, see @ref interp_ctx_t
     */
    virtual void Refine_(const interp_ctx_t* ctx) = 0;
    // /**
    //  * @brief refines by one level the information, in the presence of a resolution jump
    //  *
    //  * @param ctx the interpolation context, see @ref interp_ctx_t
    //  * @param normal the normal of the ghost layer to compute
    //  */
    // virtual void RefineGhost_(const interp_ctx_t* ctx)  = 0;
    /**
     * @brief copy the information for blocks at the same level or finer
     * 
     * @param dlvl the number of level we have to coarsen
     * @param ctx the interpolation context, see @ref interp_ctx_t
     */
    virtual void Copy_(const sid_t dlvl, const interp_ctx_t* ctx) = 0;
};

#endif  // SRC_INTERPOLATE_HPP_
