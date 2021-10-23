#include "core/memdata.hpp"

#include <cstring>

MemPtr::~MemPtr() {
#ifndef NDEBUG
    m_assert(!is_allocated, "you must free the memory first!");
#endif
};

void MemPtr::Allocate(const size_t element_size) {
    //--------------------------------------------------------------------------
    size = element_size;
    ptr  = reinterpret_cast<real_t*>(m_calloc(element_size * sizeof(real_t)));
#ifndef NDEBUG
    is_allocated = true;
#endif
    //--------------------------------------------------------------------------
};

void MemPtr::Free() {
    //--------------------------------------------------------------------------
#ifndef NDEBUG
    m_assert(is_allocated, "the memory must have been allocated before being free'ed");
    is_allocated = false;
#endif
    free(ptr);
    size = 0;
    //--------------------------------------------------------------------------
};

void MemPtr::MemSetZero() {
    //--------------------------------------------------------------------------
    std::memset(ptr, 0, size * sizeof(real_t));
    //--------------------------------------------------------------------------
}