#ifndef SRC_GHOST_BLOCK_HPP_
#define SRC_GHOST_BLOCK_HPP_

#include "boundary.hpp"
#include "murphy.hpp"
#include "subblock.hpp"

/**
 * @brief GhostBlock: a @ref SubBlock that will be used to compute ghost points
 * 
 * A GhostBlock store either ano
 * 
 */
class GhostBlock : public SubBlock {
   protected:
    sid_t dlvl_;      //!< delta level = neighbor level - my level (see @ref Interpolator)
    lid_t shift_[3];  //!<  my (0,0,0) in the framework of my neighbour  (see @ref Interpolator::interpolate_())

    // GridBlock* origin_;

    bool isghost_ = false;
    union {
        real_p     data_src_;   //!<  the ghost source comes from another rank
        GridBlock* block_src_;  //!< the ghost source comes from the same rank, so we store the @ref GridBlock
    };

   public:
    sid_t        dlvl() const { return dlvl_; }
    lid_t        shift(const int id) const { return shift_[id]; }
    const lid_t* shift() const { return shift_; }
    // GridBlock*   origin() { return origin_; }

    real_p     data_src() { return data_src_; }
    GridBlock* block_src() { return block_src_; }

    GhostBlock(GridBlock* me, const qdrt_t* ngh);
    GhostBlock(GridBlock* me, const qdrt_t* ngh, const real_t ngh_tree_offset[3], real_p data);

   protected:
    void GhostBlock_(GridBlock* me, const qdrt_t* ngh, const real_t ngh_pos[3]);
};

#endif  // SRC_GHOST_BLOCK_HPP_
