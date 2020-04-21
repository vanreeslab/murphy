#include "ghostblock.hpp"
#include "p8est_bits.h"

// GhostBlock::GhostBlock(GridBlock* me, const qdrt_t* ngh, const real_t ngh_tree_offset[3], real_p data) {
//     m_assert(m_isaligned(data), "memory is not aligned, which is not good");
//     //-------------------------------------------------------------------------
//     data_src_ = data;
//     // get the position of the current quad:
//     real_t u_size     = (1.0 / P8EST_ROOT_LEN);
//     real_t pos_ngh[3] = {ngh_tree_offset[0] + ngh->x * u_size,
//                          ngh_tree_offset[1] + ngh->y * u_size,
//                          ngh_tree_offset[2] + ngh->z * u_size};

//     m_verb("is my quad a neighbor?? %d ",p8est_quadrant_is_extended(ngh));
//     // setup the indexesd
//     GhostBlock_(me, ngh, pos_ngh);
//     //-------------------------------------------------------------------------
// }

// /**
//  * @brief Initialize the GhostBlock if ngh is a locally own quadrant
//  * 
//  * @param me the current block
//  * @param ngh a locally own quadrant
//  */
// GhostBlock::GhostBlock(GridBlock* me, const qdrt_t* ngh) {
//     //-------------------------------------------------------------------------
//     block_src_ = reinterpret_cast<GridBlock*>(ngh->p.user_data);

//     m_verb("is my quad a neighbor?? %d ",p8est_quadrant_is_extended(ngh));
//     // setup the indexes
//     GhostBlock_(me, ngh, block_src_->xyz());
//     //-------------------------------------------------------------------------
// }

GhostBlock::GhostBlock(GridBlock* me, const sid_t ngh_level, const real_t ngh_pos[3]) {
    //-------------------------------------------------------------------------
    // store my origin
    // origin_ = me;
    // get the ghost size and the
    gs_     = me->gs();
    stride_ = me->stride();
    // set the level gap > 0 if the neighbor if finer
    dlvl_ = ngh_level - me->level();

    // compute the lenght of a quadrant at the neighbor's level
    real_t len_ngh = m_quad_len(ngh_level);

    for (int id = 0; id < 3; id++) {
        // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
        real_t shift_pos = me->xyz(id) - ngh_pos[id];
        shift_[id]       = (lid_t)(shift_pos / len_ngh * M_N);
        // m_assert(floor(shift_pos / len_ngh) == (shift_pos / len_ngh), "incomplete division, this is baaad: %f / %f",shift_pos , len_ngh);
        // the start = the position of my neighbor in my frame, bounded to 0
        lid_t start_idx = (lid_t)((ngh_pos[id] - me->xyz(id)) / me->hgrid(id));
        start_[id]      = m_max(start_idx, -me->gs());
        // the end = min of how many my nieghbor can give to me and how many I can receive
        lid_t end_idx = (lid_t)((ngh_pos[id] + len_ngh - me->xyz(id)) / me->hgrid(id));
        end_[id]      = m_min(end_idx, me->end(id) + me->gs());
        // m_assert(((real_t)end_[id]) == end, "the end has to be an integer");
    }

    m_verb("shift = %d %d %d",shift_[0],shift_[1],shift_[2]);

    if (dlvl_ < 0) {
        m_verb("dbg my pos = %f %f %f", me->xyz(0), me->xyz(1), me->xyz(2));
        m_verb("dbg ngh pos = %f %f %f", ngh_pos[0], ngh_pos[1], ngh_pos[2]);
    }

    //-------------------------------------------------------------------------
}