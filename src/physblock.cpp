#include "physblock.hpp"

/**
 * @brief Construct a new Phys Block on a given face id and on a given block
 * 
 * see Boundary::operator()() for the detail of the start and end values
 * 
 * @param iface the face ID (0 <= iface < 6)
 * @param block the block on which we act
 */
PhysBlock::PhysBlock(const sid_t iface, MemLayout* block) {
    //-------------------------------------------------------------------------
    // remember the block origin
    gs_     = block->gs();
    stride_ = block->stride();

    for (int id = 0; id < 3; id++) {
        start_[id] = -block->gs();
        end_[id]   = block->stride() - block->gs();
    }

    // store the face ID
    iface_ = iface;
    m_assert(0 <= iface && iface < 6, "iface must be 0<= iface < 6");
    // overwrite in the face direction
    const sid_t dir  = iface_ / 2;
    const sid_t sign = iface_ % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus
    // see definition of @ref face_start
    // if we go minus, we need to overwrite the first point as it sits on the boundary
    start_[dir] = (sign == 0) ? (-block->gs()) : 0;
    end_[dir]   = start_[dir] + block->gs() + (sign == 0) ? 1 : 0;

    //-------------------------------------------------------------------------
}