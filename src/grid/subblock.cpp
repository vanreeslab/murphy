#include "subblock.hpp"

/**
 * @brief Construct a new Sub Block with the given ghost size, stride, start and end, see @ref MemLayout
 * 
 * @param gs the ghostsize
 * @param stride the stride
 * @param start the start index in 3D
 * @param end the end index in 3D
 */
SubBlock::SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]) {
    gs_     = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end[id];
        start_[id] = start[id];
    }
}

/**
 * @brief Construct a new Sub Block with the given ghost size, stride, start and end, see @ref MemLayout
 * where the start and end is the same in each direction
 */
SubBlock::SubBlock(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end) {
    gs_     = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end;
        start_[id] = start;
    }
}

/**
 * @brief reset the SubBlock on the given values
 */
void SubBlock::Reset(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]) {
    gs_     = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end[id];
        start_[id] = start[id];
    }
}

/**
 * @brief reset the SubBlock on the given values (same start and end index in the 3 dimensions)
 */
void SubBlock::Reset(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end) {
    gs_     = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end;
        start_[id] = start;
    }
}
