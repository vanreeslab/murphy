#ifndef SRC_CORE_MEMSPAN_HPP_
#define SRC_CORE_MEMSPAN_HPP_

#include "core/types.hpp"

struct MemLayout {
    size_t stride[2];  //!< the strides of the data: fastest rotating stride in [0], user stride in [1]
    size_t gs[2];      //!< the front ghost sizes: fastest rotating ghost size in [0], user size in [1]
    size_t n_elem;     //!< the total number of elements in the memory layout

    explicit MemLayout() = delete;
    explicit MemLayout(const lda_t n_dim, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept;
};

struct MemSpan {
    bidx_t start[3] = {0, 0, 0};  //!< starting index for the region of interest
    bidx_t end[3]   = {0, 0, 0};  //!< ending index for the region of interest

    explicit MemSpan() = delete;
    explicit MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept;
    explicit MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept;
};

#endif