#ifndef SRC_SUBBLOCK_HPP_
#define SRC_SUBBLOCK_HPP_

#include "memlayout.hpp"
#include "core/macros.hpp"
#include "core/types.hpp"

/**
 * @brief Implementation of a @ref MemLayout as any part of a memory block
 * 
 */
class SubBlock : public MemLayout {
   protected:
    lid_t gs_;
    lid_t stride_;
    lid_t start_[3];
    lid_t end_[3];

   public:
    explicit SubBlock() {}
    explicit SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]);
    explicit SubBlock(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end);

    virtual ~SubBlock(){};

    void Reset(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]);
    void Reset(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end);

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

#endif  // SRC_SUBBLOCK_HPP_
