#ifndef SRC_GHOSTBLOCK_HPP_
#define SRC_GHOSTBLOCK_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "core/memspan.hpp"
#include "wavelet/wavelet.hpp"

/**
 * @brief defines a MemSpan which happens to be used for ghosting purposes
 * 
 * -> heritage is protected to prevent anybody from calling GhostBlock::start
 */
class GhostBlock : protected MemSpan {
   protected:
    // indexing stuffs
    // Keep track of the ghost length the Ghostblock has to provide and in which direction. 
    // Examples: If my Ghostblock is a face, the is_ghost_dir_* indicates in which direction the indices have to be changed : [0,0,-2] -> [24, 24, 0] 
    bool is_ghost_dir_start_[3] = {false};   //!< indicate if we have to change the start index in a given direction (faces will change 1 direction, edges 2 and corner 3)
    bool is_ghost_dir_end_[3]   = {false};   //!< indicate if we have to change the end index in a given direction (faces will change 1 direction, edges 2 and corner 3)
    const bidx_t (*const ghost_len_)[2];  //!< ghost length information (shared by all the ghost, comes from the GridBlock

   public:
    GhostBlock() = delete;
    /**
     * @brief Construct a new Ghost Block object with the desired ghost_len
     */
    GhostBlock(const bidx_t (*const ghost_len)[2]) : ghost_len_{ghost_len} {};

    MemSpan GetSpan() const;

    void GetCoarseSpan(const MemLayout* layout, const MemLayout* layout_coarse, MemSpan* span_coarse) const;
    void GetCoarseLength(const MemLayout* layout, const MemLayout* layout_coarse, const bidx_t fine_len[3], bidx_t coarse_len[3]) const;
};

#endif  //SRC_GHOSTBLOCK_HPP_