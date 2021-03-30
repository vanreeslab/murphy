#ifndef SRC_GHOST_BLOCK_HPP_
#define SRC_GHOST_BLOCK_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/boundary.hpp"
#include "grid/subblock.hpp"

/**
 * @brief GhostBlock: a @ref SubBlock that will be used to compute ghost points.
 * 
 * It stores the relationship between two blocks, identified as source or target.
 * Although it depends on the context, the source is usually my neighbor and the target is usually me (it can be the opposite, see reversed children!)
 * 
 * @warning Because every block owns the first and the last point, there is twice the same information (might be on different levels)
 * the way to handle that depends on the difference of levels, but the rule is simple: we always trust the finer information (cfr wavelet needs)
 * 
 */
template <typename T>
class GhostBlock : public SubBlock {
   protected:
    level_t dlvl_;      //!< delta level = source level - target level (see @ref Wavelet)
    bidx_t  shift_[3];  //!< target point (0,0,0) in the framework of the source  (see @ref Wavelet::interpolate_())
    T       data_src_;  //!< provide a mean to access the source data (in practice a GridBlock* or a displacement MPI_Aint)
    rank_t  rank_rma_;  //!< the source rank of the info

   public:
    /**
    * @brief Construct a new Ghost Block object
    * 
    * We compute the intersection of two blocks and how the source one can give information to the coarse one.
    * 
    * In the case of adjacent fine source blocks, they might give redundant information to the same coarse neighbor.
    * This is ok as the sync between the fine source blocks will happen before that.
    * 
    * 
    * @param src_info_lvl the level at which the information is living, might not be the same as the level at which the information is provided
    * @param src_lvl the level of the source block
    * @param src_pos the position of the source block
    * @param src_hgrid the grid spacing of the source block
    * @param src_len the length (in the physical domain) of the source block
    * @param restrict_start true if the source must be shorter than usual, i.e. don't take the first point
    * @param trg_lvl the level of the target block
    * @param trg_pos the position of the target block
    * @param trg_hgrid the grid spacing of the target block
    * @param src_overlap_trg the overlapping of the source over the target
    * @param trg_min the minimun position of the target block
    * @param trg_max the max position of the target
    * @param trg_gs the ghost size of the target block
    * @param trg_stride the stide of the target block
    * @param rank_rma the rank needed for RMA calls
    */
    GhostBlock(/* source info */ const level_t src_lvl, const real_t src_pos[3], const real_t src_hgrid[3], const real_t src_len[3],
               /* extend info */ const bool restrict_start[3], const bool extend_end[3],
               /* target info */ const level_t trg_lvl, const real_t trg_pos[3], const real_t trg_hgrid[3],
               /* target block */ const bidx_t trg_min[3], const bidx_t trg_max[3], const bidx_t trg_gs, const bidx_t trg_stride, const rank_t rank_rma) {
        //-------------------------------------------------------------------------
        // store the source rank, the stride and the ghost size
        rank_rma_ = rank_rma;
        gs_       = trg_gs;
        stride_   = trg_stride;

        // set the level gap > 0 if the source if finer
        dlvl_ = src_lvl - trg_lvl;
        for (lda_t id = 0; id < 3; id++) {
            // the shift = (trg position - source position) expressed in the number of point of the source
            real_t shift_pos = trg_pos[id] - src_pos[id];
            shift_[id]       = (lid_t)(shift_pos / src_hgrid[id]);

            // the source start and end must be changed to match the last/first grid point policy
            // start: the first point on the source is included so we need to remove 1h (the target h!!!) if the source level < trg level
            const real_t src_start = src_pos[id] + trg_hgrid[id] * restrict_start[id];
            // end: the last point is never included, so we need to add it if the source level > target level
            const real_t src_end = src_pos[id] + src_len[id] + trg_hgrid[id] * extend_end[id];

            // // if the target is lover or the same level (dlvl > 0), nothing changes.
            // // if the target is higher (dlvl < 0), the target thrust its points more than the source
            // //      - is_src_right: the source start must be shifted by trg_hgrid (to avoid overwriting the last point)
            // //      - is_trg_right: the source end doesn't change (the first point will not be overwritten)
            // // if the source is not shrinked, the start - end will behave as usual and the last target point
            // const bool   trg_shrink   = (dlvl_ < 0);
            // const bool   is_src_right = src_pos[id] > trg_pos[id];
            // // const bool   is_trg_right = trg_pos[id] > src_pos[id];
            // const real_t sstart       = src_pos[id] + trg_hgrid[id] * (!trg_shrink && is_src_right);
            // const real_t send         = src_pos[id] + src_len[id] + trg_hgrid[id] * (trg_shrink && is_trg_right);

            // the start = the position of the souce in the trg frame, bounded to the ghost point
            const lid_t start_idx = (lid_t)((src_start - trg_pos[id]) / trg_hgrid[id]);
            start_[id]            = m_max(start_idx, trg_min[id]);

            // the end = min of how many the src can give to the trg and how many the trg can receive
            const lid_t end_idx = (lid_t)((src_end - trg_pos[id]) / trg_hgrid[id]);
            end_[id]            = m_min(end_idx, trg_max[id]);
        }
        m_assert(start_[0] <= end_[0], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[0], end_[0], shift_[0], src_pos[0], trg_pos[0]);
        m_assert(start_[1] <= end_[1], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[1], end_[1], shift_[1], src_pos[1], trg_pos[1]);
        m_assert(start_[2] <= end_[2], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[2], end_[2], shift_[2], src_pos[2], trg_pos[2]);

        m_log("ghost created: from [%d %d %d] to [%d %d %d]", start_[0], start_[1], start_[2], end_[0], end_[1], end_[2]);
        //-------------------------------------------------------------------------
    }

    // return arguments
    sid_t        dlvl() const { return dlvl_; }
    lid_t        shift(const int id) const { return shift_[id]; }
    const lid_t* shift() const { return shift_; }
    rank_t       rank() const { return rank_rma_; }

    T data_src() {
        return data_src_;
    }
    T* data_src_ptr() {
        return &data_src_;
    }
    void data_src(T data) {
        data_src_ = data;
    }
};

#endif  // SRC_GHOST_BLOCK_HPP_
