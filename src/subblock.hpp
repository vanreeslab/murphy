#ifndef SRC_SUBBLOCK_HPP_
#define SRC_SUBBLOCK_HPP_

#include "gridblock.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"

/**
 * @brief Implementation of a @ref MemLayout as any part of memory
 * 
 */
class SubBlock : public MemLayout {
   protected:
    lid_t gs_;
    lid_t stride_;
    lid_t start_[3];
    lid_t end_[3];
   public:
    SubBlock(){};
    SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]);
    /**
     * @name Memory Layout Implementation
     * 
     * @{ */
    inline lid_t gs() const override { return gs_; }
    inline lid_t stride() const override { return stride_; }
    inline lid_t start(const int id) const override { return start_[id]; }
    inline lid_t end(const int id) const override { return end_[id]; }
    /** @} */
};

#endif  // SRC_GHOSTBLOCK_HPP_