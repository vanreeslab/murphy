#ifndef SRC_SUBBLOCK_HPP_
#define SRC_SUBBLOCK_HPP_

#include "core/macros.hpp"
#include "core/layout.hpp"
#include "core/types.hpp"

/**
 * @brief Implementation of a @ref MemLayout as any part of any 3D memory block
 * 
 */
class SubBlock : public Layout {
   protected:
    lid_t start_[3] = {0, 0, 0};  //!< starting index for the region of interest
    lid_t end_[3]   = {0, 0, 0};  //!< ending index for the region of interest

   public:
    explicit SubBlock(){};
    explicit SubBlock(const lid_t start, const lid_t end);
    explicit SubBlock(const lid_t start[3], const lid_t end[3]);

    virtual ~SubBlock(){};  //!< destructor needed to guarantee the call to

    void Reset(const lid_t start, const lid_t end);
    void Reset(const lid_t start[3], const lid_t end[3]);

    void Extend(/* param */ const real_t sign[3], const bidx_t n_front, const bidx_t n_back,
                /* output */ SubBlock* new_block);

    /**
     * @name Memory Layout Implementation
     * 
     * @{ */
    inline bidx_t start(const lda_t ida) const override { return start_[ida]; }
    inline bidx_t end(const lda_t ida) const override { return end_[ida]; }
    /** @} */
};

#endif  // SRC_SUBBLOCK_HPP_
