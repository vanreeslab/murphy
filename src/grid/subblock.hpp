#ifndef SRC_SUBBLOCK_HPP_
#define SRC_SUBBLOCK_HPP_

#include "core/macros.hpp"
#include "core/memlayout.hpp"
#include "core/types.hpp"

/**
 * @brief Implementation of a @ref MemLayout as any part of any 3D memory block
 * 
 */
class SubBlock : public MemLayout {
   protected:
    lid_t gs_       = 0;          //!< the ghostsize
    lid_t stride_   = 0;          //!< the stride
    lid_t start_[3] = {0, 0, 0};  //!< starting index for the region of interest
    lid_t end_[3]   = {0, 0, 0};  //!< ending index for the region of interest

   public:
    explicit SubBlock(){};
    explicit SubBlock(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end);
    explicit SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]);

    virtual ~SubBlock(){};  //!< destructor needed to guarantee the call to

    void Reset(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end);
    void Reset(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]);

    /**
     * @name Memory Layout Implementation
     * 
     * @{ */
    inline bidx_t gs() const override { return gs_; }
    inline bidx_t stride() const override { return stride_; }
    inline bidx_t start(const lda_t ida) const override { return start_[ida]; }
    inline bidx_t end(const lda_t ida) const override { return end_[ida]; }
    /** @} */
};

#endif  // SRC_SUBBLOCK_HPP_
