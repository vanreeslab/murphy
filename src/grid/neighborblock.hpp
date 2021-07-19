#ifndef SRC_NEIGHBORBLOCK_HPP_
#define SRC_NEIGHBORBLOCK_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "ghostblock.hpp"

constexpr void face_sign(const bool mask, const iface_t iface, iface_t* face_dir, real_t sign[3]) {
    //-------------------------------------------------------------------------
    const iface_t dir = iface / 2;
    // store the dir
    (*face_dir) += static_cast<real_t>(mask) * dir;
    // update the sign
    sign[dir] += static_cast<real_t>(mask) * (((iface % 2) == 1) ? 1.0 : -1.0);
    sign[(dir + 1) % 3] += 0.0;
    sign[(dir + 2) % 3] += 0.0;
    //-------------------------------------------------------------------------
}
constexpr void edge_sign(const bool mask, const iface_t iedge, iface_t* edge_dir, real_t sign[3]) {
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
    iface_t dir = iedge / 4;  // this is the direction of the edge
    // iface_t dir1 = (dir == 0) ? 1 : 0;  // dir1 in the plane: dir1 = x if dir = y or z, or y if dir = x
    // iface_t dir2 = (dir == 2) ? 1 : 2;  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
    iface_t dir1 = static_cast<iface_t>(dir == 0);      // dir1 in the plane: dir1 = x if dir = y or z, or y if dir = x
    iface_t dir2 = 2 - static_cast<iface_t>(dir == 2);  // dir2 in the plane: dir2 = y if dir=z, = z if dir=x or dir = y
    // store the dir
    (*edge_dir) += static_cast<real_t>(mask) * dir;
    // update the sign
    sign[dir] += 0.0;
    sign[dir1] += static_cast<real_t>(mask) * (((iedge % 4) % 2) == 1 ? +1.0 : -1.0);
    sign[dir2] += static_cast<real_t>(mask) * (((iedge % 4) / 2) == 1 ? +1.0 : -1.0);
    //-------------------------------------------------------------------------
}
constexpr void corner_sign(const bool mask, const iface_t icorner, real_t sign[3]) {
    //-------------------------------------------------------------------------
    sign[0] += static_cast<real_t>(mask) * ((icorner % 2) == 1 ? +1.0 : -1.0);
    sign[1] += static_cast<real_t>(mask) * (((icorner % 4) / 2) == 1 ? +1.0 : -1.0);
    sign[2] += static_cast<real_t>(mask) * ((icorner / 4) == 1 ? +1.0 : -1.0);
    //-------------------------------------------------------------------------
}

/**
 * @brief Given a face, edge or a corner returns the outgoing normal (sign)
 * 
 * @param ibidule for a face (`ibidule<6`), an edge (`6<= ibidule < 18`) or a corner (`18<= ibidule < 26`). can also be -1, we return 0.0 then
 * @param sign the sign of the outgoing normal
 */
static void GhostGetSign(const iface_t ibidule, real_t sign[3]) {
    // m_assert(0 <= ibidule && ibidule < 26, "ibidule must be 0 <= %d < 26", ibidule);
    //-------------------------------------------------------------------------
    // check depending on the plane, the edge of the corner
    // if (ibidule < 6) {
    //     iface_t dir;
    //     face_sign(ibidule, &dir, sign);
    //     m_assert(fabs(sign[0]) + fabs(sign[1]) + fabs(sign[2]) == 1, "we cannot have more than 1 nonzero sign: %f %f %f (bidule=%d)", sign[0], sign[1], sign[2], ibidule);
    // } else if (ibidule < 18) {
    //     iface_t dir;
    //     edge_sign(ibidule - 6, &dir, sign);
    //     m_assert(fabs(sign[0]) + fabs(sign[1]) + fabs(sign[2]) == 2, "we cannot have more than 1 nonzero sign: %f %f %f (bidule=%d, dir = %d)", sign[0], sign[1], sign[2], ibidule, dir);
    // } else {
    //     corner_sign(ibidule - 18, sign);
    //     m_assert(fabs(sign[0]) + fabs(sign[1]) + fabs(sign[2]) == 3, "we cannot have more than 1 nonzero sign: %f %f %f (bidule=%d)", sign[0], sign[1], sign[2], ibidule);
    // }
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
    // indexing
    // bool scale_dir_start_[3];  //!< indicate if we have to scale in a given direction at the start index (faces scale in 1 direction, edges in 2 and )
    // bool scale_dir_end_[3];    //!< indicate if we have to scale in a given direction at the end index (faces scale in 1 direction, edges in 2 and )

    const iface_t ibidule_;  //!< the id of the face, edge or corner that is responsible for this ghostblock

    // other block
    short_t dlvl_     = 0;          //!< delta level = source level - target level (see @ref Wavelet)
    bidx_t  shift_[3] = {0, 0, 0};  //!< target point (0,0,0) in the framework of the source  (see @ref Wavelet::interpolate_())
    T       data_src_;              //!< provide a mean to access the source data (in practice a GridBlock* or a displacement MPI_Aint)

    const iblock_t cum_block_id_;  //!< cummulative id of the corresponding block or id to use to retrieve information on another rank
    const rank_t   rank_rma_;      //!< the source rank of the info

   public:
    /**
     * @brief Construct a new Ghost Block object
     * 
     * @param gs the ghost size that corresponds to the block (this is NOT the size of the ghost layer!!)
     * @param stride the stride that corresponds to the block
     * @param ibidule the id of the ghost
     * @param block_id the cummulative id of the underlying neighboring block (used to access the data)
     * @param rank the rank where to find the underlying neighboring block (-1 if local)
     */
    NeighborBlock(const bidx_t gs, const bidx_t stride, const iface_t ibidule, const iblock_t block_id, const rank_t rank = -1) : GhostBlock(gs, stride),
                                                                                                                                  ibidule_(ibidule),
                                                                                                                                  cum_block_id_(block_id),
                                                                                                                                  rank_rma_(rank) {
        //---------------------------------------------------------------------
        m_assert(this->core() == M_N, "the core = %d should be %d", this->core(), M_N);
        m_assert(this->gs(0) == M_GS && this->gs(1) == M_GS, "the gs = %d %d should be %d", this->gs(0), this->gs(1), M_GS);
        // given the ibidule, determine if I have to scale and in which direction
        // real_t sign[3];
        // GhostGetSign(ibidule, sign);
        // // I scale at the start of the block if the sign is negative
        // for (lda_t ida = 0; ida < 3; ++ida) {
        //     scale_dir_start_[ida] = sign[ida] < (-0.5);  // scale the start if the sign is negative
        //     scale_dir_end_[ida]   = 0.5 < sign[ida];     // scale the end if the sign is positive
        // }
        //---------------------------------------------------------------------
    }

    // return arguments
    [[nodiscard]] short_t      dlvl() const { return dlvl_; }
    [[nodiscard]] bidx_t       shift(const int id) const { return shift_[id]; }
    [[nodiscard]] const lid_t* shift() const { return shift_; }
    [[nodiscard]] rank_t       rank() const { return rank_rma_; }
    [[nodiscard]] iblock_t     cum_block_id() const { return cum_block_id_; };
    [[nodiscard]] iface_t      ibidule() const { return ibidule_; };

    [[nodiscard]] T data_src() {
        return data_src_;
    }
    [[nodiscard]] T* data_src_ptr() {
        return &data_src_;
    }
    void data_src(T data) {
        data_src_ = data;
    }

    // MemLayout return functions
    // [[nodiscard]] bidx_t gs() const override {
    //     return gs_;
    // };
    // [[nodiscard]] bidx_t stride() const override {
    //     return stride_;
    // };
    // [[nodiscard]] bidx_t start(const lda_t ida) const override final {
    //     return start_[ida] - scale_dir_start_[ida] * ghost_len_[0];
    // };
    // [[nodiscard]] bidx_t end(const lda_t ida) const override final {
    //     return start_[ida] + scale_dir_end_[ida] * ghost_len_[1];
    // };
    // [[nodiscard]] const bidx_t* ghost_len() const { return ghost_len_; };

    // /**
    //  * @brief transforms indexes to a coarse layout
    //  *
    //  * @param interp
    //  * @param id_fine
    //  * @param id_coarse
    //  */
    // void ToCoarse(const Wavelet* interp, const bidx_t id_fine[3], bidx_t id_coarse[3]) {
    //     // get the coarse ghost sizes
    //     const bidx_t n_front = interp->CoarseNGhostFront(ghost_len_[0]);
    //     const bidx_t n_back  = interp->CoarseNGhostBack(ghost_len_[1]);

    //     // compute the coarse indexes
    //     for (bidx_t id = 0; id < 3; ++id) {
    //         id_coarse[id] = TranslateBlockLimits(id_fine[id], M_N, n_front, M_NHALF, n_back);
    //     }
    // }

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

            // m_assert(!(scale_dir_start_[id] && (end_idx - start_idx) != M_GS), "if we have to scale in this direction, the ghost must be MAXIMUM. here: %d - %d vs %d", end_idx, start_idx, M_GS);
            // m_assert(!(scale_dir_end_[id] && (end_idx - start_idx) != M_GS), "if we have to scale in this direction, the ghost must be MAXIMUM. here: %d - %d vs %d", end_idx, start_idx, M_GS);

            // if we need to scale, we save the other index (if in the start, save the end and vice-versa)
            // start_[id] = scale_dir_start_[id] ? end_idx : start_idx;
            // end_[id]   = scale_dir_end_[id] ? start_idx : end_idx;

            // check that if the start is negative, it's at least the whole ghosting region
            m_assert(!(start_idx < 0 && start_idx != trg_min), "if the start index = %d is negative, if must cover the whole ghost region", start_idx);
            m_assert(!(end_idx > trg_core && end_idx != trg_max), "if the start index = %d is negative, if must cover the whole ghost region", end_idx);
            // if one of the indexes is outside of the block, we reset it to the closest boundary and register it for ghost scaling

            // determine if we need to scale the direction in the start or the end
            scale_dir_start_[id] = (start_idx == trg_min);
            scale_dir_end_[id]   = (end_idx == trg_max);

            // if yes, register the closest block boundary instead
            start_[id] = scale_dir_start_[id] ? 0 : start_idx;
            end_[id]   = scale_dir_end_[id] ? trg_core : end_idx;
        }
        m_assert(start_[0] <= end_[0], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[0], end_[0], shift_[0], src_pos[0], trg_pos[0]);
        m_assert(start_[1] <= end_[1], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[1], end_[1], shift_[1], src_pos[1], trg_pos[1]);
        m_assert(start_[2] <= end_[2], "the starting index must be < the ending index, here: %d <= %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[2], end_[2], shift_[2], src_pos[2], trg_pos[2]);
        //---------------------------------------------------------------------
    }
};

#endif  // SRC_GHOST_BLOCK_HPP_
