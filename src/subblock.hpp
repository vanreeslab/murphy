#ifndef SRC_SUBBLOCK_HPP_
#define SRC_SUBBLOCK_HPP_

#include "gridblock.hpp"
#include "memlayout.hpp"
#include "murphy.hpp"

/**
 * @brief Implementation of a @ref MemLayout as a sub part of a @ref GridBlock
 * 
 */
class SubBlock : public MemLayout {
   protected:
    lid_t gs_;
    lid_t stride_;
    lid_t start_[3];
    lid_t range_[3];

    // store the block I depend on
    GridBlock* origin_;

   public:
    /**
     * @name Memory Layout Implementation
     * 
     * @{ */
    inline lid_t gs() const override { return gs_; }
    inline lid_t stride() const override { return stride_; }
    inline lid_t start(const int id) const override { return start_[id]; }
    inline lid_t range(const int id) const override { return range_[id]; }
    /** @} */

    inline lid_t level() const { return origin_->level(); }
    
    GridBlock* origin() { return origin_; }
};

#endif  // SRC_GHOSTBLOCK_HPP_