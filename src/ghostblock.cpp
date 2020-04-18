#include "ghostblock.hpp"

/**
 * @brief Initialize the GhostBlock if ngh is a ghost quadrant
 * 
 * @param me the current block
 * @param ngh a ghost quadrant
 * @param ngh_tree_offset the position of the left corner of the tree it belongs
 * @param data the associated memory location
 */
GhostBlock::GhostBlock(GridBlock* me, const qdrt_t* ngh, const real_t ngh_tree_offset[3], real_p data) {
    m_assert(m_isaligned(data), "memory is not aligned, which is not good");
    //-------------------------------------------------------------------------
    data_src_ = data;
    // get the position of the current quad:
    real_t u_size     = (1.0 / P8EST_ROOT_LEN);
    real_t pos_ngh[3] = {ngh_tree_offset[0] + ngh->x * u_size,
                         ngh_tree_offset[1] + ngh->y * u_size,
                         ngh_tree_offset[2] + ngh->z * u_size};
    // setup the indexesd
    GhostBlock_(me, ngh, ngh_tree_offset);
    //-------------------------------------------------------------------------
}

/**
 * @brief Initialize the GhostBlock if ngh is a locally own quadrant
 * 
 * @param me the current block
 * @param ngh a locally own quadrant
 */
GhostBlock::GhostBlock(GridBlock* me, const qdrt_t* ngh) {
    //-------------------------------------------------------------------------
    block_src_ = reinterpret_cast<GridBlock*>(ngh->p.user_data);
    // setup the indexes
    GhostBlock_(me, ngh, block_src_->xyz());
    //-------------------------------------------------------------------------
}

void GhostBlock::GhostBlock_(GridBlock* me, const qdrt_t* ngh, const real_t pos_ngh[3]) {
    //-------------------------------------------------------------------------
    // store my origin
    // origin_ = me;
    // get the ghost size and the
    gs_     = me->gs();
    stride_ = me->stride();
    // set the level gap > 0 if the neighbor if finer
    dlvl_ = ngh->level - me->level();

    // compute the lenght of a quadrant at the neighbor's level
    real_t len_ngh = m_quad_len(ngh->level);

    for (int id = 0; id < 3; id++) {
        // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
        real_t shift_pos = me->xyz(id) - pos_ngh[id];
        shift_[id]       = (lid_t)(shift_pos / len_ngh * M_N);
        // m_assert(floor(shift_pos / len_ngh) == (shift_pos / len_ngh), "incomplete division, this is baaad: %f / %f",shift_pos , len_ngh);
        // the start = the position of my neighbor in my frame, bounded to 0
        lid_t start_idx = (lid_t)((pos_ngh[id] - me->xyz(id)) / me->hgrid(id));
        start_[id]      = m_max(start_idx, -me->gs());
        // the end = min of how many my nieghbor can give to me and how many I can receive
        lid_t end_idx = (lid_t)((pos_ngh[id] + len_ngh - me->xyz(id)) / me->hgrid(id));
        end_[id]      = m_min(end_idx, me->end(id) + me->gs());
        // m_assert(((real_t)end_[id]) == end, "the end has to be an integer");
    }


    if(dlvl_ < 0){
        m_verb("dbg my pos = %f %f %f",me->xyz(0),me->xyz(1),me->xyz(2));
        m_verb("dbg ngh pos = %f %f %f",pos_ngh[0],pos_ngh[1],pos_ngh[2]);

    }

    //-------------------------------------------------------------------------
}