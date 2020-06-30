#include "ghostblock.hpp"
#include "p8est_bits.h"

/**
 * @brief Construct a new Ghost Block. Computes the area that is at the intersection between a neighboring block and me.
 * 
 * @warning this computation assumes that the nieghbor is a full block of size M_N. In practise, this assumption should always be valid
 * 
 * @param me the current block
 * @param ngh_level the neighbor level
 * @param ngh_pos the origin position of my neighbor
 */
GhostBlock::GhostBlock(GridBlock* me, const sid_t ngh_level, const real_t ngh_pos[3], const bool ngh_to_me) {
    //-------------------------------------------------------------------------
    // get the ghost size and the
    gs_     = me->gs();
    stride_ = me->stride();

    // compute the lenght of a quadrant at the neighbor's level
    real_t len_ngh = m_quad_len(ngh_level);
    real_t len_me  = m_quad_len(me->level());

    if (ngh_to_me) {
        // set the level gap > 0 if the neighbor if finer
        dlvl_ = ngh_level - me->level();
        for (int id = 0; id < 3; id++) {
            // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
            real_t shift_pos = me->xyz(id) - ngh_pos[id];
            shift_[id]       = (lid_t)(shift_pos / len_ngh * M_N);
            // the start = the position of my neighbor in my frame, bounded to 0
            lid_t start_idx = (lid_t)((ngh_pos[id] - me->xyz(id)) / me->hgrid(id));
            start_[id]      = m_max(start_idx, -me->gs());
            // the end = min of how many my nieghbor can give to me and how many I can receive
            lid_t end_idx = (lid_t)((ngh_pos[id] + len_ngh - me->xyz(id)) / me->hgrid(id));
            end_[id]      = m_min(end_idx, me->end(id) + me->gs());
        }
    } else {
        // set the level gap > 0 if the I am finer
        dlvl_ = me->level() - ngh_level;
        for (int id = 0; id < 3; id++) {
            // the shift = (the neighbor position - my position) expressed in my number of point
            real_t shift_pos = ngh_pos[id] - me->xyz(id);
            shift_[id]       = (lid_t)(shift_pos / me->hgrid(id));
            // the start = the position of me in my neighbor's frame, bounded to 0
            lid_t start_idx = (lid_t)((me->xyz(id) - ngh_pos[id]) / len_me * M_N);
            start_[id]      = m_max(start_idx, -me->gs());
            // the end = min of how many my neighbor can give to me and how many I can receive
            lid_t end_idx = (lid_t)((me->xyz(id) + len_ngh - ngh_pos[id]) / len_me * M_N);
            end_[id]      = m_min(end_idx, me->end(id) + me->gs());
        }
    }
    //-------------------------------------------------------------------------
}

// GhostBlock::GhostBlock(GridBlock* me, GhostBlock* gb, const sid_t ngh_level, const real_t ngh_pos[3]) {
//     //-------------------------------------------------------------------------
//     // get the ghost size and the
//     gs_     = me->gs();
//     stride_ = me->stride();
//     // set the level gap > 0 if the neighbor if finer
//     dlvl_ = ngh_level - me->level();

//     // compute the lenght of a quadrant at the neighbor's level
//     real_t len_ngh = m_quad_len(ngh_level);

//     for (int id = 0; id < 3; id++) {
//         // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
//         real_t shift_pos = me->xyz(id) - ngh_pos[id];
//         shift_[id]       = (lid_t)(shift_pos / len_ngh * M_N);
//         // m_assert(floor(shift_pos / len_ngh) == (shift_pos / len_ngh), "incomplete division, this is baaad: %f / %f",shift_pos , len_ngh);
//         // the start = the position of my neighbor in my frame, bounded to 0
//         lid_t start_idx = (lid_t)((ngh_pos[id] - me->xyz(id)) / me->hgrid(id));
//         start_[id]      = m_max(start_idx, -me->gs());
//         // the end = min of how many my nieghbor can give to me and how many I can receive
//         lid_t end_idx = (lid_t)((ngh_pos[id] + len_ngh - me->xyz(id)) / me->hgrid(id));
//         end_[id]      = m_min(end_idx, me->end(id) + me->gs());
//         // m_assert(((real_t)end_[id]) == end, "the end has to be an integer");
//     }
//     //-------------------------------------------------------------------------
// }

// GhostBlock* GhostBlock::CreateOpposite(GridBlock* me,const real_t ngh_pos[3]){
//     //-------------------------------------------------------------------------
//     // create an empty version of a new GhostBlock
//     GhostBlock* gb = new GhostBlock();

//     // set the gs to be the same as mine
//     gs
//     //-------------------------------------------------------------------------
// }