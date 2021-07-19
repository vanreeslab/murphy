#include "physblock.hpp"

/**
 * @brief Construct a new Phys Block on a given face id and on a given block
 * 
 * we set the scale_dir_start_/end_ all true but the start and end index are different depending on the iface
 * 
 * @param iface the face ID (0 <= iface < 6)
 * @param block the block on which we act
 */
PhysBlock::PhysBlock(const iface_t iface, const MemLayout* block) : GhostBlock(block->gs(), block->stride()),
                                                                    iface_(iface) {
    m_assert(block->gs() == M_GS, "the ghost size %d must be %d", block->gs(), M_GS);
    m_assert(block->stride() == M_STRIDE, "the ghost size %d must be %d", block->stride(), M_STRIDE);
    m_assert(0 <= iface && iface < 6, "iface must be 0<= iface < 6");
    //-------------------------------------------------------------------------
    // get the start and end index in all direction
    for (lda_t id = 0; id < 3; id++) {
        // set the start and end
        start_[id] = block->start(id);
        end_[id]   = block->end(id);

        scale_dir_start_[id] = true;
        scale_dir_end_[id]   = true;

        m_assert(start_[id] == 0, "the start %d must be %d", start_[id], 0);
        m_assert(end_[id] == M_N, "the start %d must be %d", end_[id], M_STRIDE);
    }

    // get info on dir and sign of the face
    const iface_t dir  = iface_ / 2;
    const iface_t sign = iface_ % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus

    // update the needed coordinates in the dir direction (we overwrite the first info if we go minus!)
    start_[dir] = (sign == 0) ? start_[dir] : (M_N);
    end_[dir]   = (sign == 0) ? 1 : end_[dir];  // 1 to ensure that we visit the 0th point

    // only scale in the direction of the face!
    scale_dir_start_[dir] = (sign == 0);
    scale_dir_end_[dir]   = (sign == 1);
    //-------------------------------------------------------------------------
}