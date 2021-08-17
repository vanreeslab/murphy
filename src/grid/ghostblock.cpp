#include "ghostblock.hpp"
// /**
//  * @brief Translate the ghost limits from one block to another
//  * 
//  * if the current index is in a ghost region (the front one or the back one), returns the limit of the new ghost region
//  * if the current index is in the block, return the index scaled to the new length of the block
//  * 
//  * @warning if the new index is NOT an integer, we ceil it! (sometimes we use an index of 1 because we want to visit the index 0, we should preserve that)
//  * 
//  * @param c_id current index in the current reference
//  * @param c_core the number of points inside the current block (typically M_N)
//  * @param n_front the number of ghost points in front of the block in the new block
//  * @param n_core the number of points inside the new block
//  * @param n_back the number of ghost points at the back of the block in the new block
//  * @return bidx_t the new ghost index
//  */
// static bidx_t TranslateBlockLimits(const bidx_t c_id, const bidx_t c_core,
//                                    const bidx_t n_front, const bidx_t n_core, const bidx_t n_back) {
//     // if ratio > 1, we scale up, if ratio < 1, we scale down
//     const real_t ratio = static_cast<real_t>(c_core) / static_cast<real_t>(n_core);
//     // get where we are (in front, in the block or at the back)
//     const bidx_t b = (c_id + c_core);
//     const bidx_t c = (b / c_core) + static_cast<bidx_t>(c_id > c_core);
//     // compute every possibility
//     const bidx_t res[4] = {-n_front, static_cast<bidx_t>(ceil(c_id / ratio)), n_core, n_core + n_back};
//     // check that we didn't screw up the indexes
//     // m_assert(!(c == 1 && !m_fequal(static_cast<real_t>(c_id) / ratio, res[1])), "we cannot translate the id %d from a core of %d to a core of %d: %f == %d", c_id, c_core, n_core, static_cast<real_t>(c_id) / ratio, res[1]);
//     // return the correct choice
//     return res[c];
// }

//==============================================================================

/**
 * @brief returns the MemSpan attached to my ghost region
 * 
 * The ghost region is changing its size and is driven by scale_dir_start_ given at construction
 * 
 * @return MemSpan the CURRENT region covered by the ghosts
 */
MemSpan GhostBlock::GetSpan() const {
    // create an empty Span
    MemSpan me;
    // populate is given the current ghost size
#pragma unroll
    for (lda_t ida = 0; ida < 3; ++ida) {
        me.start[ida] = start[ida] - static_cast<bidx_t>(scale_dir_start_[ida]) * (*ghost_len_)[0];
        me.end[ida]   = end[ida] + static_cast<bidx_t>(scale_dir_end_[ida]) * (*ghost_len_)[1];
    }
    return me;
}

void GhostBlock::GetCoarseSpan(const MemLayout* layout, const MemLayout* layout_coarse, MemSpan* span_coarse) const {
    // get the current span (warning cannot be the start or end arrays)
    MemSpan me = GetSpan();

    // we assume that the new layout contains the new desired ghost sizes!
#pragma unroll 3
    for (lda_t ida = 0; ida < 3; ++ida) {
        span_coarse->start[ida] = layout->TranslateLimit(layout_coarse, me.start[ida]);
        span_coarse->end[ida]   = layout->TranslateLimit(layout_coarse, me.end[ida]);
    }

    // // get the coarse ghost sizes
    // const bidx_t coarse_gs[2] = {interp->CoarseNGhostFront((*ghost_len_)[0]),
    //                              interp->CoarseNGhostBack((*ghost_len_)[1])};
    // this->Resize(coarse_gs, M_NHALF, block_coarse);
}

void GhostBlock::GetCoarseLength(const MemLayout* layout, const MemLayout* layout_coarse, const bidx_t fine_len[3], bidx_t coarse_len[3]) const {
    // we assume that the new layout contains the new desired ghost sizes!
#pragma unroll 3
    for (lda_t ida = 0; ida < 3; ++ida) {
        coarse_len[ida] = layout->TranslateLimit(layout_coarse, fine_len[ida]);
    }

    // // get the coarse ghost sizes
    // const bidx_t coarse_gs[2] = {interp->CoarseNGhostFront((*ghost_len_)[0]),
    //                              interp->CoarseNGhostBack((*ghost_len_)[1])};
    // this->Resize(coarse_gs, M_NHALF, block_coarse);
}