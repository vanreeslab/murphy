#include "subblock.hpp"

/**
 * @brief Construct a new Sub Block with the given ghost size, stride, start and end, see @ref MemLayout
 * 
 * @param gs the ghostsize
 * @param stride the stride
 * @param start the start index in 3D
 * @param end the end index in 3D
 */
SubBlock::SubBlock(const lid_t start[3], const lid_t end[3]) {
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end[id];
        start_[id] = start[id];
    }
}

/**
 * @brief Construct a new Sub Block with the given ghost size, stride, start and end, see @ref MemLayout
 * where the start and end is the same in each direction
 */
SubBlock::SubBlock(const lid_t start, const lid_t end) {
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end;
        start_[id] = start;
    }
}

/**
 * @brief reset the SubBlock on the given values
 */
void SubBlock::Reset(const lid_t start[3], const lid_t end[3]) {
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end[id];
        start_[id] = start[id];
    }
}

/**
 * @brief reset the SubBlock on the given values (same start and end index in the 3 dimensions)
 */
void SubBlock::Reset(const lid_t start, const lid_t end) {
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end;
        start_[id] = start;
    }
}

/**
 * @brief extend the subblock and return the result in new_block
 * 
 * the extension is based on @ref sign: if <0 we extend the start index, if >0 we extend the end index
 * 
 * @param sign 
 * @param n_front number of point for the extension in front, used if sign < 0
 * @param n_back number of point for the extension at the back, used if sign > 0
 * @param new_block 
 */
void SubBlock::Extend(/* param */ const real_t sign[3], const bidx_t n_front, const bidx_t n_back,
                      /* output */ SubBlock* new_block) {
    //-------------------------------------------------------------------------
    bidx_t start[3] = {this->start(0), this->start(1), this->start(2)};
    bidx_t end[3]   = {this->end(0), this->end(1), this->end(2)};

    for (lda_t ida = 0; ida < 3; ++ida) {
        if (sign[ida] > 0.5) {
            end[ida] += n_back;  // extend outside the block
        } else if (sign[ida] < (-0.5)) {
            start[ida] -= n_front;  // extend outside the block
        }
    }
    m_log("extension from %d %d %d -> %d %d %d to %d %d %d -> %d %d %d", this->start(0), this->start(1), this->start(2), this->end(0), this->end(1), this->end(2), start[0], start[1], start[2], end[0], end[1], end[2]);
    // set the new block to the computed start/end
    new_block->Reset(start, end);
    //-------------------------------------------------------------------------
}
