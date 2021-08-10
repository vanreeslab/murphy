#include "core/data.hpp"

//==============================================================================
static size_t MemPadSize(const size_t size, const size_t size_type) {
    m_assert((M_ALIGNMENT % size_type) == 0, "the alignement must be a mutiple of %d", size_type);
    //--------------------------------------------------------------------------
    const size_t chunk = M_ALIGNMENT / size_type;  //alignement in terms of T
    const size_t mod   = (size % chunk);
    // pad and return
    return size - mod + (mod > 0) * chunk;
    //--------------------------------------------------------------------------
};

//==============================================================================
MemPtr::~MemPtr() {
#ifndef NDEBUG
    m_assert(!is_allocated, "you must free the memory first!");
#endif
};

void MemPtr::Allocate(const size_t element_size) {
    ptr = reinterpret_cast<real_t*>(aligned_alloc(M_ALIGNMENT, element_size * sizeof(real_t)));
#ifndef NDEBUG
    size         = element_size;
    is_allocated = true;
#endif
};

void MemPtr::Free() {
    free(ptr);
#ifndef NDEBUG
    m_assert(is_allocated, "the memory must have been allocated before being free'ed");
    is_allocated = false;
#endif
};

//==============================================================================
MemLayout::MemLayout(const bidx_t in_n_gs_front, const bidx_t in_n_block, const bidx_t in_n_gs_back = -1) noexcept {
    //--------------------------------------------------------------------------
    n_g_frt = in_n_gs_front;
    n_block = in_n_block;
    n_g_bck = (in_n_gs_back > -1) ? in_n_gs_back : in_n_gs_front;
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

//==============================================================================
MemData::MemData(const MemLayout& layout, const MemPtr& ptr) noexcept {
    //--------------------------------------------------------------------------
    // we want to ensure that the offset is aligned in memory.
    // the offset is given by
    //      gs_front[0] + stride[0] *(....)
    // as the raw pointer is aligned, we must have gs_front and stride aligned!

    // 1. get the ghost sizes, only the front one matters here!
    const bidx_t gs[2] = {MemPadSize(layout.n_g_frt, sizeof(real_t)),
                          layout.n_g_frt};

    // 2. get the strides
    const bidx_t stride = layout.n_g_frt + layout.n_block + layout.n_g_bck;
    stride_[0]          = MemPadSize(stride, sizeof(real_t));
    stride_[1]          = stride;

    // get the shifted pointer
    const size_t offset = gs[0] + stride_[0] * (gs[1] + stride_[1] * gs[1]);
    data_               = ptr.ptr + offset;
    m_assert(m_isaligned(data_), "the value of data must be aligned!");
    //--------------------------------------------------------------------------
}