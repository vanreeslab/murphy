#include "core/memspan.hpp"
#include "macros.hpp"

static bidx_t MemPadSize(const bidx_t size, const bidx_t size_type) {
    m_assert((M_ALIGNMENT % size_type) == 0, "the alignement must be a mutiple of %d", size_type);
    //--------------------------------------------------------------------------
    const bidx_t chunk = M_ALIGNMENT / size_type;  //alignement in terms of T
    const bidx_t mod   = (size % chunk);
    // pad and return
    return size - mod + (mod > 0) * chunk;
    //--------------------------------------------------------------------------
};

MemLayout::MemLayout(const m_layout_t dim, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back = -1) noexcept {
    //--------------------------------------------------------------------------
    // we want to ensure that the offset is aligned in memory.
    // the offset is given by
    //      gs_front[0] + stride[0] *(....)
    // as the raw pointer is aligned, we must have gs_front and stride aligned!

    // 1. get the ghost sizes and the shift to apply
    gs    = n_gs_front;
    shift = MemPadSize(gs, sizeof(real_t)) - gs;
    m_assert(MemPadSize(gs, sizeof(real_t)) >= gs, "the padded size %ld must be >= the actual size %ld", MemPadSize(gs, sizeof(real_t)), gs);

    // 2. get the strides
    const bidx_t gs_back    = (n_gs_back > -1) ? n_gs_back : n_gs_front;
    const bidx_t stride_usr = n_gs_front + n_block + gs_back;  // we use here the new gs_back!
    stride[0]               = MemPadSize(stride_usr, sizeof(real_t));
    stride[1]               = stride_usr;

    // compute the number of elements as the product of the respective strides + the shift
    n_elem = 1;
    for (lda_t ida = 0; ida < (dim - 1); ++ida) {
        n_elem = stride[1] * n_elem;
    }
    // fastest dimension
    n_elem = shift + stride[0] * n_elem;
    //--------------------------------------------------------------------------
}

//==============================================================================
MemSpan::MemSpan(const bidx_t in_start, const bidx_t in_end) noexcept {
    //--------------------------------------------------------------------------
#pragma unroll 3
    for (lda_t ida = 0; ida < 3; ++ida) {
        start[ida] = in_start;
        end[ida]   = in_end;
    }
    //--------------------------------------------------------------------------
}

MemSpan::MemSpan(const bidx_t in_start[3], const bidx_t in_end[3]) noexcept {
    //--------------------------------------------------------------------------
#pragma unroll 3
    for (lda_t ida = 0; ida < 3; ++ida) {
        start[ida] = in_start[ida];
        end[ida]   = in_end[ida];
    }
    //--------------------------------------------------------------------------
}