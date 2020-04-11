#ifndef SRC_INTERPOLATE_HPP_
#define SRC_INTERPOLATE_HPP_

#include "field.hpp"
#include "murphy.hpp"
#include "p8est.h"
#include "memlayout.hpp"

class Interpolator {
   protected:
    lid_t srcgs_;
    lid_t trggs_;
    lid_t srcstr_;
    lid_t trgstr_;
    lid_t trgstart_[3];
    lid_t trgend_[3];
    lid_t srcstart_[3];
    lid_t srcend_[3];

   public:
    
   protected:

   /**
    * @brief interpolate the data from data_src to data_trg
    * 
    * @param dlvl the difference of level: level_src - level_trg, i.e. > 0 means coarsening, = 0 means copy and < 0 means refinement
    * @param shift the position of the trg (0,0,0) in the src framework (and resolution!)
    * @param block_src description of @ref data_src memory layout 
    * @param data_src the actual src memory
    * @param block_trg descripiton of @ref data_trg
    * @param data_trg the actual target memory
    */
    virtual void interpolate_(const lid_t dlvl, const lid_t shift[3], MemLayout* block_src, real_p data_src, MemLayout* block_trg, real_p data_trg);

    
    virtual void Coarsen_(const lid_t dlvl, const real_p sdata, real_p tdata) const = 0;
    virtual void Refine_(const real_p sdata, real_p tdata) const                    = 0;
    virtual void Copy_(const real_p sdata, real_p tdata) const                      = 0;
};

#endif  // SRC_INTERPOLATE_HPP_