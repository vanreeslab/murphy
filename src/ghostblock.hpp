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
template <typename T>
class GhostBlock : public SubBlock {
   protected:
    sid_t dlvl_;      //!< delta level = neighbor level - my level (see @ref Interpolator)
    lid_t shift_[3];  //!<  my (0,0,0) in the framework of my neighbour  (see @ref Interpolator::interpolate_())

    rank_t rank_src_;  //!< the source rank of the info
    T      data_src_;  //!< provide a mean to access the source data (in practice a GridBlock* or a displacement MPI_Aint)

   public:
    /**
    * @brief Construct a new Ghost Block< T>:: Ghost Block object
    * 
    * @warning this computation assumes that the nieghbor is a full block of size M_N. In practise, this assumption should always be valid
    * 
    * @tparam T 
    * @param me 
    * @param ngh_level 
    * @param ngh_pos 
    * @param nghost_front 
    * @param nghost_back 
    * @param ngh_to_me 
    * @param rank_src 
    */
    GhostBlock(GridBlock* me, const sid_t ngh_level, const real_t ngh_pos[3], const sid_t nghost_front[3], const sid_t nghost_back[3], const bool ngh_to_me, const rank_t rank_src) {
        //-------------------------------------------------------------------------
        // store the source rank
        rank_src_ = rank_src;

        // get the ghost size and the
        gs_     = me->gs();
        stride_ = me->stride();

        // compute the lenght of a quadrant at the neighbor's level
        real_t len_ngh = m_quad_len(ngh_level);
        real_t len_me  = m_quad_len(me->level());

        if (ngh_to_me) {
            // set the level gap > 0 if the neighbor if finer
            dlvl_ = ngh_level - me->level();
            for (int id = 0; id < 3; id++) {
                // the shift = (my position - the neighbor position) expressed in the number of point in my neighbor
                real_t shift_pos = me->xyz(id) - ngh_pos[id];
                shift_[id]       = (lid_t)(shift_pos / len_ngh * M_N);
                // the start = the position of my neighbor in my frame, bounded to 0
                lid_t start_idx = (lid_t)((ngh_pos[id] - me->xyz(id)) / me->hgrid(id));
                start_[id]      = m_max(start_idx, -nghost_front[id]);
                // the end = min of how many my nieghbor can give to me and how many I can receive
                lid_t end_idx = (lid_t)((ngh_pos[id] + len_ngh - me->xyz(id)) / me->hgrid(id));
                end_[id]      = m_min(end_idx, me->end(id) + nghost_back[id]);
            }
        } else {
            // set the level gap > 0 if the I am finer
            dlvl_ = me->level() - ngh_level;
            for (int id = 0; id < 3; id++) {
                // the shift = (the neighbor position - my position) expressed in my number of point
                real_t shift_pos = ngh_pos[id] - me->xyz(id);
                shift_[id]       = (lid_t)(shift_pos / me->hgrid(id));
                // the start = the position of me in my neighbor's frame, bounded to 0
                lid_t start_idx = (lid_t)((me->xyz(id) - ngh_pos[id]) / len_me * M_N);
                start_[id]      = m_max(start_idx, -nghost_front[id]);
                // the end = min of how many my neighbor can give to me and how many I can receive
                lid_t end_idx = (lid_t)((me->xyz(id) + len_ngh - ngh_pos[id]) / len_me * M_N);
                end_[id]      = m_min(end_idx, me->end(id) + nghost_back[id]);
            }
        }
        //-------------------------------------------------------------------------
    }

    GhostBlock(GridBlock* me, const sid_t ngh_level, const real_t ngh_pos[3], const sid_t nghost_front[3], const sid_t nghost_back[3], const bool ngh_to_me) : GhostBlock(me, ngh_level, ngh_pos, nghost_front, nghost_back, ngh_to_me, -1) {}

    // return arguments
    sid_t        dlvl() const { return dlvl_; }
    lid_t        shift(const int id) const { return shift_[id]; }
    const lid_t* shift() const { return shift_; }
    rank_t       rank() const { return rank_src_; }

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
