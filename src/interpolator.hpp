#ifndef SRC_INTERPOLATE_HPP_
#define SRC_INTERPOLATE_HPP_

#include "field.hpp"
#include "murphy.hpp"
#include "p8est.h"
#include "subblock.hpp"

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
    virtual void operator()(SubBlock* block_src, GhostBlock* block_trg);
    /**
     * @brief reconstuct missing block_trg data from block_src data
     * 
     * @param dlvl the level difference: src_level - trg_level
     * @param shift the position of the (0,0,0) from block_trg in block_src framework (and resolution!!)
     * @param block_src 
     * @param block_trg 
     */
    virtual void operator()(const lid_t dlvl, const lid_t shift[3], SubBlock* block_src, SubBlock* block_trg);

   protected:
    virtual void Coarsen_(const lid_t dlvl, const real_p sdata, real_p tdata) const = 0;
    virtual void Refine_(const real_p sdata, real_p tdata) const                    = 0;
    virtual void Copy_(const real_p sdata, real_p tdata) const                      = 0;
};

#endif  // SRC_INTERPOLATE_HPP_