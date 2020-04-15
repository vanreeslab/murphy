#include "physblock.hpp"

PhysBlock::PhysBlock(const sid_t iface, GridBlock* block) {
    //-------------------------------------------------------------------------
    // remember the block origin
    gs_     = block->gs();
    stride_ = block->stride();

    for (int id = 0; id < 3; id++) {
        normal_sign_[id] = 0;
        start_[id]       = -block->gs();
        end_[id]       = block->stride();
    }

    // overwrite in the face direction
    const sid_t dir   = iface / 2;
    const sid_t sign  = iface % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus
    end_[dir]       = block->gs();
    start_[dir]       = (sign == 0) ? (-block->gs()) : (block->stride()-block->gs());
    normal_sign_[dir] = (sign == 0) ? -1 : 1;
    //-------------------------------------------------------------------------
}