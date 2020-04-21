#include "physblock.hpp"

PhysBlock::PhysBlock(const sid_t iface, GridBlock* block) {
    //-------------------------------------------------------------------------
    // remember the block origin
    gs_     = block->gs();
    stride_ = block->stride();

    for (int id = 0; id < 3; id++) {
        start_[id] = -block->gs();
        end_[id]   = block->stride()-block->gs();
    }

    // store the face ID
    iface_ = iface;
    // overwrite in the face direction
    const sid_t dir  = iface_ / 2;
    const sid_t sign = iface_ % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus
    // see definition of Boundary for acessing
    start_[dir]      = (sign == 0) ? (-block->gs()) : 1;
    end_[dir]        = start_[dir] + block->gs();

    //-------------------------------------------------------------------------
}