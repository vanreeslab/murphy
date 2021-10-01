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
    bool scale_dir_start_[3] = {false};   //!< indicate if we have to scale in a given direction at the start index (faces scale in 1 direction, edges in 2 and )
    bool scale_dir_end_[3]   = {false};   //!< indicate if we have to scale in a given direction at the end index (faces scale in 1 direction, edges in 2 and )
    const bidx_t (*const ghost_len_)[2];  //!< ghost length informatin (shared by all the ghost, comes from the GridBlock

   public:
    GhostBlock() = delete;
    GhostBlock(const bidx_t (*const ghost_len)[2]) : ghost_len_{ghost_len} {};

    MemSpan GetSpan() const;

    void GetCoarseSpan(const MemLayout* layout, const MemLayout* layout_coarse, MemSpan* span_coarse) const;
    void GetCoarseLength(const MemLayout* layout, const MemLayout* layout_coarse, const bidx_t fine_len[3], bidx_t coarse_len[3]) const;

    // [[nodiscard]] bidx_t start(const lda_t ida) const{
    //     return start[ida] - static_cast<bidx_t>(scale_dir_start_[ida]) * (*ghost_len_)[0];
    // };
    // [[nodiscard]] bidx_t end(const lda_t ida) const{
    //     return end[ida] + static_cast<bidx_t>(scale_dir_end_[ida]) * (*ghost_len_)[1];
    // };

    // /**
    //  * @brief transforms indexes to a coarse layout
    //  *
    //  * @param interp
    //  * @param id_fine
    //  * @param id_coarse
    //  */
    // void ToCoarse(const Wavelet* interp, const bidx_t id_fine[3], bidx_t id_coarse[3]) {
    //     // get the coarse ghost sizes
    //     const bidx_t n_front = interp->CoarseNGhostFront((*ghost_len_)[0]);
    //     const bidx_t n_back  = interp->CoarseNGhostBack((*ghost_len_)[1]);

    //     // compute the coarse indexes
    //     for (bidx_t id = 0; id < 3; ++id) {
    //         id_coarse[id] = TranslateBlockLimits(id_fine[id], M_N, n_front, M_NHALF, n_back);
    //     }
    // }
};

#endif  //SRC_GHOSTBLOCK_HPP_