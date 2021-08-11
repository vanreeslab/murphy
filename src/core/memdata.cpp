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
