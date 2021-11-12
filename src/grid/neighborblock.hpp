#ifndef SRC_NEIGHBORBLOCK_HPP_
#define SRC_NEIGHBORBLOCK_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/ghostblock.hpp"

constexpr void face_sign(const bool mask_bool, const iface_t iface, iface_t* face_dir, real_t sign[3]) {
    //-------------------------------------------------------------------------
    const real_t  mask = static_cast<real_t>(mask_bool);
    const iface_t dir  = (iface / 2) * mask;
    // store the dir
    (*face_dir) += mask * dir;
    // update the sign
    sign[dir] += mask * (((iface % 2) == 1) ? 1.0 : -1.0);
    sign[(dir + 1) % 3] += 0.0;
    sign[(dir + 2) % 3] += 0.0;
    //-------------------------------------------------------------------------
}

constexpr void edge_sign(const bool mask_bool, const iface_t iedge, iface_t* edge_dir, real_t sign[3]) {
    /*
    the plane convention for the sign variable convention for the sign
    2 +--------------+ 3
      |              |
      |              |
      |dir2          |
      |              |
    0 +--------------+ 1
        dir1
    */
    //-------------------------------------------------------------------------
    const real_t  mask = static_cast<real_t>(mask_bool);
    const iface_t dir  = (iedge / 4) * mask;                  // this is the direction of the edge
    const iface_t dir1 = static_cast<iface_t>(dir == 0);      // dir1 in the plane: dir1 = x if dir = y or z, or y if dir = x
    const iface_t dir2 = 2 - static_cast<iface_t>(dir == 2);  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
    // store the dir
    (*edge_dir) += mask * dir;
    // update the sign
    sign[dir] += 0.0;
    sign[dir1] += mask * (((iedge % 4) % 2) == 1 ? +1.0 : -1.0);
    sign[dir2] += mask * (((iedge % 4) / 2) == 1 ? +1.0 : -1.0);
    //-------------------------------------------------------------------------
}

constexpr void corner_sign(const bool mask_bool, const iface_t icorner, real_t sign[3]) {
    //-------------------------------------------------------------------------
    const real_t mask = static_cast<real_t>(mask_bool);
    sign[0] += mask * ((icorner % 2) == 1 ? +1.0 : -1.0);
    sign[1] += mask * (((icorner % 4) / 2) == 1 ? +1.0 : -1.0);
    sign[2] += mask * ((icorner / 4) == 1 ? +1.0 : -1.0);
    //-------------------------------------------------------------------------
}

/**
 * @brief Given a face, edge or a corner returns the outgoing normal (sign)
 * 
 * @param ibidule for a face (`ibidule<6`), an edge (`6<= ibidule < 18`) or a corner (`18<= ibidule < 26`). can also be -1, we return 0.0 then
 * @param sign the sign of the outgoing normal
 */
static void GhostGetSign(const iface_t ibidule, real_t sign[3]) {
    //-------------------------------------------------------------------------
    iface_t dir = 0;
    for (iface_t i = 0; i < 3; ++i) {
        sign[i] = 0.0;
    }
    face_sign((0 <= ibidule) && (ibidule < 6), ibidule, &dir, sign);
    edge_sign((6 <= ibidule) && (ibidule < 18), ibidule - 6, &dir, sign);
    corner_sign((18 <= ibidule) && (ibidule < 26), ibidule - 18, sign);

    m_assert(sign[0] == 0.0 || sign[0] == 1.0 || sign[0] == -1.0, "wrong sign value: %e", sign[0]);
    m_assert(sign[1] == 0.0 || sign[1] == 1.0 || sign[1] == -1.0, "wrong sign value: %e", sign[1]);
    m_assert(sign[2] == 0.0 || sign[2] == 1.0 || sign[2] == -1.0, "wrong sign value: %e", sign[2]);
    //-------------------------------------------------------------------------
};

/**
 * @brief NeighborBlock: a @ref SubBlock that will be used to compute ghost points.
 * 
 * It stores the relationship between two blocks, identified as source or target.
 * Although it depends on the context, the source is usually my neighbor and the target is usually me (it can be the opposite, see reversed children!)
 * 
 * @warning Because every block owns the first and the last point, there is twice the same information (might be on different levels)
 * the way to handle that depends on the difference of levels, but the rule is simple: we always trust the finer information (cfr wavelet needs)
 * 
 */
template <typename T>
class NeighborBlock : public GhostBlock {
   private:
    const iface_t ibidule_;  //!< the id of the face, edge or corner that is responsible for this ghostblock

    // other block
    short_t dlvl_     = 0;          //!< delta level = source level - target level (see @ref Wavelet)
    bidx_t  shift_[3] = {0, 0, 0};  //!< target point (0,0,0) in the framework of the source  (see @ref Wavelet::interpolate_())
    T       data_src_;              //!< provide a mean to access the source data (in practice a GridBlock* or a displacement MPI_Aint)

    const iblock_t cum_block_id_;  //!< cummulative id of the corresponding block or id to use to retrieve information on another rank
    const rank_t   rank_rma_;      //!< the source rank of the info

   public:
    /**
     * @brief Construct a new Neighbor Block object
     * 
     * @param gs the ghost size that corresponds to the block (this is NOT the size of the ghost layer!!)
     * @param stride the stride that corresponds to the block
     * @param ghost_len the pointer to the ghost_len information
     * @param ibidule the id of the ghost
     * @param block_id the cummulative id of the underlying neighboring block (used to access the data)
     * @param rank the rank where to find the underlying neighboring block (-1 if local)
     */
    NeighborBlock(const bidx_t (*const ghost_len)[2],
                  const iface_t ibidule, const iblock_t block_id, const rank_t rank = -1) : GhostBlock(ghost_len),
                                                                                            ibidule_(ibidule),
                                                                                            cum_block_id_(block_id),
                                                                                            rank_rma_(rank) {
        //---------------------------------------------------------------------
        //---------------------------------------------------------------------
    }

    // return arguments
    M_INLINE short_t       dlvl() const { return dlvl_; }
    M_INLINE bidx_t        shift(const int id) const { return shift_[id]; }
    M_INLINE const bidx_t* shift() const { return shift_; }
    M_INLINE rank_t        rank() const { return rank_rma_; }
    M_INLINE iblock_t      cum_block_id() const { return cum_block_id_; };
    M_INLINE iface_t       ibidule() const { return ibidule_; };

    M_INLINE T data_src() const noexcept {
        return data_src_;
    }
    M_INLINE T* data_src_ptr() noexcept {
        return &data_src_;
    }
    M_INLINE void data_src(T data) {
        data_src_ = data;
    }
    M_INLINE MemSpan SourceSpan() const noexcept {
        return data_src_->BlockSpan();
    }
    M_INLINE MemLayout SourceLayout() const noexcept {
        return data_src_->BlockLayout();
    }

    /**
    * @brief Computes the ghost region as the intersection between two blocks
    * 
    * We compute the intersection of two blocks and how the source one can give information to the coarse one.
    * 
    * In the case of adjacent fine source blocks, they might give redundant information to the same coarse neighbor.
    * This is ok as the sync between the fine source blocks will happen before that.
    * 
    * 
    * @param src_lvl the level of the source block
    * @param src_pos the position of the source block
    * @param src_hgrid the grid spacing of the source block
    * @param src_len the length (in the physical domain) of the source block
    * @param trg_lvl the level of the target block
    * @param trg_pos the position of the target block
    * @param trg_hgrid the grid spacing of the target block
    * @param trg_min the minimun position of the target block
    * @param trg_max the max position of the target
    */
    void Intersect(/* source info */ const level_t src_lvl, const real_t src_pos[3], const real_t src_hgrid[3], const real_t src_len[3],
                   /* target info */ const level_t trg_lvl, const real_t trg_pos[3], const real_t trg_hgrid[3],
                   /* target index */ const bidx_t trg_ghost_len[2], const bidx_t trg_core) {  //}, const bidx_t trg_gs, const bidx_t trg_stride) {  //, const rank_t rank_rma) {
        //---------------------------------------------------------------------
        const bidx_t trg_min = -trg_ghost_len[0];
        const bidx_t trg_max = trg_core + trg_ghost_len[1];
        // set the level gap > 0 if the source if finer
        dlvl_ = src_lvl - trg_lvl;
        for (lda_t id = 0; id < 3; id++) {
            // the shift = (trg position - source position) expressed in the number of point of the source
            real_t shift_pos = trg_pos[id] - src_pos[id];
            shift_[id]       = static_cast<bidx_t>(shift_pos / src_hgrid[id]);

            // the source start and end must be changed to match the last/first grid point policy
            // start: the first point on the source is included so we need to remove 1h (the target h!!!) if the source level < trg level
            const real_t src_start = src_pos[id];
            const real_t src_end   = src_pos[id] + src_len[id];

            // the start = the position of the souce in the trg frame, bounded to the ghost point
            const bidx_t start_idx = m_max(static_cast<bidx_t>((src_start - trg_pos[id]) / trg_hgrid[id]), trg_min);
            // the end = min of how many the src can give to the trg and how many the trg can receive
            const bidx_t end_idx = m_min(static_cast<bidx_t>((src_end - trg_pos[id]) / trg_hgrid[id]), trg_max);

            // check that if the start is negative, it's at least the whole ghosting region
            m_assert(!(start_idx < 0 && start_idx != trg_min), "if the start index = %d is negative, if must cover the whole ghost region", start_idx);
            m_assert(!(end_idx > trg_core && end_idx != trg_max), "if the start index = %d is negative, if must cover the whole ghost region", end_idx);
            // if one of the indexes is outside of the block, we reset it to the closest boundary and register it for ghost scaling

            // determine if we need to change the start or the end index in the specified direction 
            is_ghost_dir_start_[id] = (start_idx == trg_min);
            is_ghost_dir_end_[id]   = (end_idx == trg_max);

            // if yes, register the closest block boundary instead
            start[id] = is_ghost_dir_start_[id] ? 0 : start_idx;
            end[id]   = is_ghost_dir_end_[id] ? trg_core : end_idx;
        }
        m_assert(start[0] <= end[0], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start[0], end[0], shift_[0], src_pos[0], trg_pos[0]);
        m_assert(start[1] <= end[1], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start[1], end[1], shift_[1], src_pos[1], trg_pos[1]);
        m_assert(start[2] <= end[2], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start[2], end[2], shift_[2], src_pos[2], trg_pos[2]);
        //---------------------------------------------------------------------
    }
};

// specify the functions when we have a MPI_Aint
template <>
M_INLINE MemSpan NeighborBlock<MPI_Aint>::SourceSpan() const noexcept {
    return MemSpan(0, M_N);
}
// specify the function for the BlockSpan function
template <>
M_INLINE MemLayout NeighborBlock<MPI_Aint>::SourceLayout() const noexcept {
    return MemLayout(M_LAYOUT_BLOCK, M_GS, M_N);
}

#endif  // SRC_GHOST_BLOCK_HPP_
