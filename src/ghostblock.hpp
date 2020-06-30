#ifndef SRC_GHOST_BLOCK_HPP_
#define SRC_GHOST_BLOCK_HPP_

#include "boundary.hpp"
#include "murphy.hpp"
#include "subblock.hpp"

/**
 * @brief GhostBlock: a @ref SubBlock that will be used to compute ghost points.
 * 
 * 
 * It adds to the @ref SubBloc: the delta level, the shift between the blocks and memory pointer.
 * A GhostBlock store either a local @ref GridBlock pointer or the buffer memory it has to use
 * 
 */
class GhostBlock : public SubBlock {
   protected:
    sid_t dlvl_;      //!< delta level = neighbor level - my level (see @ref Interpolator)
    lid_t shift_[3];  //!<  my (0,0,0) in the framework of my neighbour  (see @ref Interpolator::interpolate_())

    union {
        real_p     data_src_;   //!<  the ghost source comes from another rank
        GridBlock* block_src_;  //!< the ghost source comes from the same rank, so we store the @ref GridBlock
    };

   public:
    GhostBlock(GridBlock* me, const sid_t ngh_level, const real_t ngh_pos[3], const sid_t nghost[3], const bool ngh_to_me);

    // return arguments
    sid_t        dlvl() const { return dlvl_; }
    lid_t        shift(const int id) const { return shift_[id]; }
    const lid_t* shift() const { return shift_; }
    real_p       data_src() { return data_src_; }
    GridBlock*   block_src() { return block_src_; }

    void data_src(real_p data) { data_src_ = data; }
    void block_src(GridBlock* block) { block_src_ = block; }
};

#endif  // SRC_GHOST_BLOCK_HPP_
