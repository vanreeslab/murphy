#include "physblock.hpp"

/**
 * @brief Construct a new Phys Block on a given face id and on a given block
 * 
 * we set the is_ghost_dir_start_/end_ all true but the start and end index are different depending on the iface
 * 
 * @param ghost_len the pointer to the ghost_len information
 * @param iface the face ID (0 <= iface < 6)
 * @param block the block on which we act
 */
PhysBlock::PhysBlock(const bidx_t (*const ghost_len)[2],
                     const iface_t iface, const MemSpan* block) : GhostBlock(ghost_len),
                                                                  iface_(iface) {
    m_assert(0 <= iface && iface < 6, "iface must be 0<= iface < 6");
    //-------------------------------------------------------------------------
    // get the start and end index in all direction
#pragma unroll 3
    for (lda_t id = 0; id < 3; ++id) {
        // set the start and end
        start[id] = block->start[id];
        end[id]   = block->end[id];

        is_ghost_dir_start_[id] = true;
        is_ghost_dir_end_[id]   = true;

        m_assert(start[id] == 0, "the start %d must be %d", start[id], 0);
        m_assert(end[id] == M_N, "the start %d must be %d", end[id], M_N);
    }

    // get info on dir and sign of the face
    const iface_t dir  = iface_ / 2;
    const iface_t sign = iface_ % 2;  // sign = 1, -> we go plus, sign = 0 -> we go minus

    // update the needed coordinates in the dir direction (we overwrite the first info if we go minus!)
    start[dir] = (sign == 0) ? start[dir] : (M_N);
    end[dir]   = (sign == 0) ? 1 : end[dir];  // 1 to ensure that we visit the 0th point

    // only scale in the direction of the face!
    is_ghost_dir_start_[dir] = (sign == 0);
    is_ghost_dir_end_[dir]   = (sign == 1);
    //-------------------------------------------------------------------------
}