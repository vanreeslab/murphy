#ifndef SRC_GHOST_BLOCK_HPP_
#define SRC_GHOST_BLOCK_HPP_

#include "boundary.hpp"
#include "defs.hpp"
#include "subblock.hpp"

/**
 * @brief GhostBlock: a @ref SubBlock that will be used to compute ghost points.
 * 
 * 
 * It adds to the @ref SubBloc: the delta level, the shift between the blocks and memory pointer.
 * A GhostBlock store either a local @ref GridBlock pointer or the buffer memory it has to use
 * 
 */
template <typename T>
class GhostBlock : public SubBlock {
   protected:
    level_t dlvl_;      //!< delta level = neighbor level - my level (see @ref Wavelet)
    lid_t   shift_[3];  //!<  my (0,0,0) in the framework of my neighbour  (see @ref Wavelet::interpolate_())
    rank_t  rank_rma_;  //!< the source rank of the info
    T       data_src_;  //!< provide a mean to access the source data (in practice a GridBlock* or a displacement MPI_Aint)

   public:
   /**
    * @brief Construct a new Ghost Block object
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
    * @param trg_gs the ghost size of the target block
    * @param trg_stride the stide of the target block
    * @param rank_rma the rank needed for RMA calls
    */
    GhostBlock(const level_t src_lvl, const real_t src_pos[3], const real_t src_hgrid[3], const real_t src_len[3],
               const level_t trg_lvl, const real_t trg_pos[3], const real_t trg_hgrid[3],
               const lid_t trg_min[3], const lid_t trg_max[3], const lid_t trg_gs, const lid_t trg_stride, const rank_t rank_rma) {
        //-------------------------------------------------------------------------
        // store the source rank, the stride and the ghost size
        rank_rma_ = rank_rma;
        gs_       = trg_gs;
        stride_   = trg_stride;

        // set the level gap > 0 if the neighbor if finer
        dlvl_ = src_lvl - trg_lvl;
        for (int id = 0; id < 3; id++) {
            // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
            real_t shift_pos = trg_pos[id] - src_pos[id];
            shift_[id]       = (lid_t)(shift_pos / src_hgrid[id]);
            // the start = the position of the souce in the trg frame, bounded to the ghost point
            lid_t start_idx = (lid_t)((src_pos[id] - trg_pos[id]) / trg_hgrid[id]);
            start_[id]      = m_max(start_idx, trg_min[id]);
            // the end = min of how many the src can give to the trg and how many the trg can receive
            lid_t end_idx = (lid_t)((src_pos[id] + src_len[id] - trg_pos[id]) / trg_hgrid[id]);
            end_[id]      = m_min(end_idx, trg_max[id]);
        }
        m_assert(start_[0] < end_[0], "the starting index must be < the ending index, here: %d < %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[0], end_[0], shift_[0], src_pos[0], trg_pos[0]);
        m_assert(start_[1] < end_[1], "the starting index must be < the ending index, here: %d < %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[1], end_[1], shift_[1], src_pos[1], trg_pos[1]);
        m_assert(start_[2] < end_[2], "the starting index must be < the ending index, here: %d < %d, the shift = %d: the src_pos = %f, trg_pos = %f", start_[2], end_[2], shift_[2], src_pos[2], trg_pos[2]);
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
