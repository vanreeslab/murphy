#include "physblock.hpp"

/**
 * @brief Construct a new Phys Block on a given face id and on a given block
 * 
 * see Boundary::operator()() for the detail of the start and end values
 * 
 * @param iface the face ID (0 <= iface < 6)
 * @param block the block on which we act
 */
PhysBlock::PhysBlock(const iface_t iface, const Layout* block, const lid_t nghost_front, const lid_t nghost_back) {
    m_assert(0 <= iface && iface < 6, "iface must be 0<= iface < 6");
    //-------------------------------------------------------------------------
    // remember the block origin
    // gs_     = block->gs();
    // stride_ = block->stride();

    for (int id = 0; id < 3; id++) {
        start_[id] = -nghost_front;
        end_[id]   = M_N + nghost_back;
        // m_assert((end_[id] - start_[id]) <= stride_, "the face is too big for the stride");
    }

    // store the face ID
    iface_ = iface;
    // overwrite in the face direction
    const iface_t dir  = iface_ / 2;
    const iface_t sign = iface_ % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus

    // update the needed coordinates in the dir direction (we overwrite the first info if we go minus)
    start_[dir] = (sign == 0) ? start_[dir] : (M_N);
    end_[dir]   = (sign == 0) ? 1 : end_[dir];  // 1 to ensure that we visit the 0th point
    //-------------------------------------------------------------------------
}