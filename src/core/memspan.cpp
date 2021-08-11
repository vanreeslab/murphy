#include "core/memspan.hpp"

//==============================================================================
MemLayout::MemLayout(const lda_t dim, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept {
    //--------------------------------------------------------------------------
    // we want to ensure that the offset is aligned in memory.
    // the offset is given by
    //      gs_front[0] + stride[0] *(....)
    // as the raw pointer is aligned, we must have gs_front and stride aligned!

    // 1. get the ghost sizes, only the front one matters here!
    gs[0] = MemPadSize(n_gs_front, sizeof(real_t));
    gs[1] = n_gs_front;

    // 2. get the strides
    const bidx_t gs_back    = (n_gs_back > -1) ? n_gs_back : n_gs_front;
    const bidx_t stride_usr = n_gs_front + n_block + gs_back;  // we use here the new gs_back!
    stride[0]               = MemPadSize(stride_usr, sizeof(real_t));
    stride[1]               = stride_usr;

    // compute the number of elements as the maximum offset that will be used + 1
    // the offset is given by
    //      (gs[0] + id) + stride[0] * ((gs[1] + id) + stride[1] * (gs[1] + id))
    // where id = (n_block + gs_back)-1
    n_elem = 0;
    for (lda_t ida = 0; ida < (dim - 1); ++ida) {
        n_elem = (stride[1]-1) + stride[1] * n_elem ;
    }
    // the last dimension
    n_elem = (gs[0]+n_block+gs_back) + stride[0] * n_elem;
    //--------------------------------------------------------------------------
}

//==============================================================================
MemSpan::MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept {
    //--------------------------------------------------------------------------
#pragma unroll
    for (lda_t ida = 0; ida < 3; ++ida) {
        start[ida] = in_start;
        end[ida]   = in_end;
    }
    //--------------------------------------------------------------------------
}

MemSpan::MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept {
    //--------------------------------------------------------------------------
#pragma unroll
    for (lda_t ida = 0; ida < 3; ++ida) {
        start[ida] = in_start[ida];
        end[ida]   = in_end[ida];
    }
    //--------------------------------------------------------------------------
}