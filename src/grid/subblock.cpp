#include "subblock.hpp"

/**
 * @brief Construct a new Sub Block with the given ghost size and stride, start and end are set to 0
 * 
 * @param gs the ghostsize
 * @param stride the stride
 */
SubBlock::SubBlock(const lid_t gs, const lid_t stride) {
    gs_[0]  = gs;
    gs_[1]  = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = 0;
        start_[id] = 0;
    }
}
/**
 * @brief Construct a new Sub Block with the given ghost size and stride, start and end are set to 0
 * 
 * @param gs the ghostsize: front and back
 * @param stride the stride
 */
SubBlock::SubBlock(const lid_t gs[2], const lid_t stride) {
    gs_[0]  = gs[0];
    gs_[1]  = gs[1];
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = 0;
        start_[id] = 0;
    }
}

/**
 * @brief Construct a new Sub Block with the given ghost size, stride, start and end
 * 
 * @param gs the ghost size
 * @param stride the stride
 * @param start the start index in 3D
 * @param end the end index in 3D
 */
SubBlock::SubBlock(const lid_t gs, const lid_t stride, const lid_t start[3], const lid_t end[3]) {
    gs_[0]  = gs;
    gs_[1]  = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end[id];
        start_[id] = start[id];
    }
}

/**
 * @brief Construct a new Sub Block with the given ghost size, stride, start and end
 * 
 * @param gs the ghostsizes: front and back
 * @param stride the stride
 * @param start the start index in 3D
 * @param end the end index in 3D
 */
SubBlock::SubBlock(const lid_t gs[2], const lid_t stride, const lid_t start[3], const lid_t end[3]) {
    gs_[0]  = gs[0];
    gs_[1]  = gs[1];
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
    gs_[0]  = gs;
    gs_[1]  = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end;
        start_[id] = start;
    }
}

/**
 * @brief reset the SubBlock on the given values (same start and end index in the 3 dimensions)
 */
void SubBlock::Reset(const lid_t gs, const lid_t stride, const lid_t start, const lid_t end) {
    //-------------------------------------------------------------------------
    gs_[0]  = gs;
    gs_[1]  = gs;
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        end_[id]   = end;
        start_[id] = start;
    }
    //-------------------------------------------------------------------------
}
/**
 * @brief reset the SubBlock on the given values
 */
void SubBlock::Reset(const lid_t gs[2], const lid_t stride, const lid_t start[3], const lid_t end[3]) {
    //-------------------------------------------------------------------------
    gs_[0]  = gs[0];
    gs_[1]  = gs[1];
    stride_ = stride;
    for (lda_t id = 0; id < 3; id++) {
        start_[id] = start[id];
        end_[id]   = end[id];
    }
    //-------------------------------------------------------------------------
}

// /**
//  * @brief extend the subblock and return the result in new_block
//  *
//  * the extension is based on @ref sign: if <0 we extend the start index, if >0 we extend the end index
//  *
//  * @param sign
//  * @param n_front number of point for the extension in front, used if sign < 0
//  * @param n_back number of point for the extension at the back, used if sign > 0
//  * @param new_block
//  */
// void SubBlock::Extend(/* param */ const real_t sign[3], const bidx_t n_front, const bidx_t n_back,
//                       /* output */ SubBlock* new_block) {
//     //-------------------------------------------------------------------------
//     bidx_t start[3] = {this->start(0), this->start(1), this->start(2)};
//     bidx_t end[3]   = {this->end(0), this->end(1), this->end(2)};

//     for (lda_t ida = 0; ida < 3; ++ida) {
//         if (sign[ida] > 0.5) {
//             end[ida] += n_back;  // extend outside the block
//         } else if (sign[ida] < (-0.5)) {
//             start[ida] -= n_front;  // extend outside the block
//         }
//         // else {
//         //     start[ida] -= n_front;
//         //     end[ida] += n_back;
//         // }
//     }
//     m_log("extension from %d %d %d -> %d %d %d to %d %d %d -> %d %d %d", this->start(0), this->start(1), this->start(2), this->end(0), this->end(1), this->end(2), start[0], start[1], start[2], end[0], end[1], end[2]);
//     // set the new block to the computed start/end
//     new_block->Reset(this->gs(), this->stride(), start, end);
//     //-------------------------------------------------------------------------
// }

void SubBlock::Resize(const bidx_t new_gs[2], const bidx_t new_core, SubBlock* new_block) {
    //-------------------------------------------------------------------------
    // get some info
    const bidx_t current_core = this->core();
    // compute the coarse indexes
    bidx_t new_start[3], new_end[3];
    for (bidx_t id = 0; id < 3; ++id) {
        new_start[id] = TranslateBlockLimits(start(id), current_core, new_gs[0], new_core, new_gs[1]);
        new_end[id]   = TranslateBlockLimits(end(id), current_core, new_gs[0], new_core, new_gs[1]);
    }
    // update the SubBlock
    const bidx_t new_stride = new_gs[0] + new_gs[1] + new_core;
    new_block->Reset(new_gs, new_stride, new_start, new_end);
    //-------------------------------------------------------------------------
}
