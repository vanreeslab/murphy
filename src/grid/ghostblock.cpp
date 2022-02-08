#include "ghostblock.hpp"
/**
 * @brief returns the MemSpan attached to my ghost region
 * 
 * The ghost region is changing its size and is driven by is_ghost_dir_*_ given at construction
 * 
 * @return MemSpan the CURRENT region covered by the ghosts
 */
MemSpan GhostBlock::GetSpan() const {
    // create an empty Span
    MemSpan me;
    // populate is given the current ghost size
#pragma unroll
    for (lda_t ida = 0; ida < 3; ++ida) {
        me.start[ida] = start[ida] - static_cast<bidx_t>(is_ghost_dir_start_[ida]) * (*ghost_len_)[0];
        me.end[ida]   = end[ida] + static_cast<bidx_t>(is_ghost_dir_end_[ida]) * (*ghost_len_)[1];
    }
    return me;
}

/**
 * @brief returns a MemSpan corresponding to the coarse representation of my MemSpan
 * 
 */
void GhostBlock::GetCoarseSpan(const MemLayout* layout, const MemLayout* layout_coarse, MemSpan* span_coarse) const {
    // get the current span (warning cannot be the start or end arrays)
    MemSpan me = GetSpan();

    // we assume that the new layout contains the new desired ghost sizes!
#pragma unroll 3
    for (lda_t ida = 0; ida < 3; ++ida) {
        span_coarse->start[ida] = layout->TranslateLimit(layout_coarse, me.start[ida]);
        span_coarse->end[ida]   = layout->TranslateLimit(layout_coarse, me.end[ida]);
    }
}

/**
 * @brief Compute the coarse indexes based on fine ones 
 * 
 * @param layout 
 * @param layout_coarse 
 * @param fine_len 
 * @param coarse_len 
 */
void GhostBlock::GetCoarseLength(const MemLayout* layout, const MemLayout* layout_coarse, const bidx_t fine_len[3], bidx_t coarse_len[3]) const {
    // we assume that the new layout contains the new desired ghost sizes!
#pragma unroll 3
    for (lda_t ida = 0; ida < 3; ++ida) {
        coarse_len[ida] = layout->TranslateLimit(layout_coarse, fine_len[ida]);
    }
}