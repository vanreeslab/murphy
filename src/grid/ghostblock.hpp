#ifndef SRC_GHOSTBLOCK_HPP_
#define SRC_GHOSTBLOCK_HPP_

#include "core/macros.hpp"
#include "core/types.hpp"
#include "grid/subblock.hpp"
#include "wavelet/wavelet.hpp"

class GhostBlock{
   protected:
    // indexing stuffs
    bool scale_dir_start_[3] = {false};   //!< indicate if we have to scale in a given direction at the start index (faces scale in 1 direction, edges in 2 and )
    bool scale_dir_end_[3]   = {false};   //!< indicate if we have to scale in a given direction at the end index (faces scale in 1 direction, edges in 2 and )
    const bidx_t (*const ghost_len_)[2];  //!< ghost length informatin (shared by all the ghost, comes from the GridBlock

   public:
    GhostBlock() = delete;
    GhostBlock(const bidx_t gs, const bidx_t stride, const bidx_t (*const ghost_len)[2]) : SubBlock(gs, stride),
                                                                                           ghost_len_{ghost_len} {};

    [[nodiscard]] bidx_t start(const lda_t ida) const{
        return start_[ida] - static_cast<bidx_t>(scale_dir_start_[ida]) * (*ghost_len_)[0];
    };
    [[nodiscard]] bidx_t end(const lda_t ida) const{
        return end_[ida] + static_cast<bidx_t>(scale_dir_end_[ida]) * (*ghost_len_)[1];
    };

    /**
     * @brief set block_coarse as a coarser layout the current block
     * 
     * @param interp 
     * @param block_coarse 
     */
    void ToCoarse(const Wavelet* interp, SubBlock* block_coarse) {
        // get the coarse ghost sizes
        const bidx_t coarse_gs[2] = {interp->CoarseNGhostFront((*ghost_len_)[0]),
                                     interp->CoarseNGhostBack((*ghost_len_)[1])};
        this->Resize(coarse_gs, M_NHALF, block_coarse);
    }

    /**
     * @brief transforms indexes to a coarse layout
     * 
     * @param interp 
     * @param id_fine 
     * @param id_coarse 
     */
    void ToCoarse(const Wavelet* interp, const bidx_t id_fine[3], bidx_t id_coarse[3]) {
        // get the coarse ghost sizes
        const bidx_t n_front = interp->CoarseNGhostFront((*ghost_len_)[0]);
        const bidx_t n_back  = interp->CoarseNGhostBack((*ghost_len_)[1]);

        // compute the coarse indexes
        for (bidx_t id = 0; id < 3; ++id) {
            id_coarse[id] = TranslateBlockLimits(id_fine[id], M_N, n_front, M_NHALF, n_back);
        }
    }
};

#endif  //SRC_GHOSTBLOCK_HPP_