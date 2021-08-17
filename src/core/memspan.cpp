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

MemLayout::MemLayout(const m_layout_t dim, const bidx_t n_gs_front, const bidx_t n_block, const bidx_t n_gs_back) noexcept {
    //--------------------------------------------------------------------------
    // we want to ensure that the offset is aligned in memory.
    // the offset is given by
    //      gs_front[0] + stride[0] *(....)
    // as the raw pointer is aligned, we must have gs_front and stride aligned!

    // 1. get the ghost sizes and the shift to apply
    gs    = n_gs_front;
    shift = MemPadSize(gs, sizeof(real_t)) - gs;
    m_assert(MemPadSize(gs, sizeof(real_t)) >= gs, "the padded size %d must be >= the actual size %d", MemPadSize(gs, sizeof(real_t)), gs);

    // 2. get the strides
    const bidx_t gs_back    = (n_gs_back > -1) ? n_gs_back : n_gs_front;
    const bidx_t stride_usr = n_gs_front + n_block + gs_back;  // we use here the new gs_back!
    stride[0]               = MemPadSize(stride_usr, sizeof(real_t));
    stride[1]               = stride_usr;
    block                   = n_block;

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

//==============================================================================
/**
 * @brief translate the ghost limits from the current layout to another one
 * 
 * if the current index is in a ghost region (the front one or the back one), returns the limit of the new ghost region (known from the new layout)
 * if the current index is in the block, return the index scaled to the new length of the block
 * 
 * @warning if the new index is NOT an integer, we ceil it! (sometimes we use an index of 1 because we want to visit the index 0, we should preserve that)
 * 
 * @param new_layout 
 * @param id index in the current reference 
 * @return bidx_t 
 */
bidx_t MemLayout::TranslateLimit(const MemLayout* new_layout, const bidx_t c_id) const noexcept {
    //--------------------------------------------------------------------------
    const bidx_t n_front = new_layout->gs;
    const bidx_t n_core  = new_layout->block;
    const bidx_t n_back  = new_layout->stride[1] - (new_layout->block + new_layout->gs);
    const bidx_t c_core  = block;

    // let's go for some dark magic
    //..........................................................................
    // if ratio > 1, we scale up, if ratio < 1, we scale down
    const real_t ratio = static_cast<real_t>(c_core) / static_cast<real_t>(n_core);
    // get where we are (in front, in the block or at the back)
    const bidx_t b = (c_id + c_core);
    const bidx_t c = (b / c_core) + static_cast<bidx_t>(c_id > c_core);
    // compute every possibility
    const bidx_t res[4] = {-n_front, static_cast<bidx_t>(ceil(c_id / ratio)), n_core, n_core + n_back};
    // return the correct choice
    return res[c];

    //--------------------------------------------------------------------------
}