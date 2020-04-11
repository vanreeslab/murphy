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
    m_begin;
    m_assert(m_isaligned(data), "memory is not aligned, which is not good");
    //-------------------------------------------------------------------------
    origin_   = me;
    data_src_ = data;
    // get the position of the current quad:
    real_t u_size     = (1.0 / P8EST_ROOT_LEN);
    real_t pos_ngh[3] = {ngh_tree_offset[0] + ngh->x * u_size,
                         ngh_tree_offset[1] + ngh->y * u_size,
                         ngh_tree_offset[2] + ngh->z * u_size};
    // setup the indexes
    GhostBlock_(me, ngh, ngh_tree_offset);
    //-------------------------------------------------------------------------
    m_end;
}

/**
 * @brief Initialize the GhostBlock if ngh is a locally own quadrant
 * 
 * @param me the current block
 * @param ngh a locally own quadrant
 */
GhostBlock::GhostBlock(GridBlock* me, const qdrt_t* ngh) {
    m_begin;
    //-------------------------------------------------------------------------
    origin_    = me;
    block_src_ = reinterpret_cast<GridBlock*>(ngh->p.user_data);
    // setup the indexes
    GhostBlock_(me, ngh, block_src_->xyz());
    //-------------------------------------------------------------------------
    m_end;
}

void GhostBlock::GhostBlock_(GridBlock* me, const qdrt_t* ngh, const real_t pos_ngh[3]) {
    m_begin;
    //-------------------------------------------------------------------------
    // get the ghost size and the
    gs_     = M_GS;
    stride_ = M_STRIDE;
    // set the level gap > 0 if the neighbor if finer
    dlvl_ = ngh->level - me->level();

    // compute the lenght of a quadrant at the neighbor's level
    real_t len_ngh = m_quad_len(ngh->level );    

    for (int id = 0; id < 3; id++) {
        // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
        shift_[id] = ((me->xyz(id) - pos_ngh[id]) / len_ngh) * M_N;
        // the start = the position of my neighbor in my frame, bounded to 0
        start_[id] = m_max(pos_ngh[id] - me->xyz(id), 0.0) / me->hgrid(id);
        // the range = how many my nieghbor can give to me - how many I can receive
        range_[0] = m_min(len_ngh / me->hgrid(id), M_STRIDE - start_[0]);
    }
    //-------------------------------------------------------------------------
    m_end;
}