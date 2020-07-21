#include "physblock.hpp"

/**
 * @brief Construct a new Phys Block on a given face id and on a given block
 * 
 * see Boundary::operator()() for the detail of the start and end values
 * 
 * @param iface the face ID (0 <= iface < 6)
 * @param block the block on which we act
 */
PhysBlock::PhysBlock(const sid_t iface, MemLayout* block, const sid_t nghost_front[3], const sid_t nghost_back[3]) {
    //-------------------------------------------------------------------------
    // remember the block origin
    gs_     = block->gs();
    stride_ = block->stride();

    for (int id = 0; id < 3; id++) {
        start_[id] = -nghost_front[id];
        end_[id]   = block->stride() - nghost_front[id];
    }

    // store the face ID
    iface_ = iface;
    m_assert(0 <= iface && iface < 6, "iface must be 0<= iface < 6");
    // overwrite in the face direction
    const sid_t dir  = iface_ / 2;
    const sid_t sign = iface_ % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus

    // update the needed coordinates in the dir direction (we overwrite the first info if we go minus)
    start_[dir] = (sign == 0) ? start_[dir] : M_N;
    end_[dir]   = (sign == 0) ? 1 : end_[dir];
    //-------------------------------------------------------------------------
}