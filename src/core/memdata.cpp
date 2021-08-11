#include "core/memdata.hpp"

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
