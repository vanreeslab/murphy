#ifndef SRC_CORE_MEMDATA_HPP_
#define SRC_CORE_MEMDATA_HPP_

#include <functional>

#include "core/macros.hpp"
#include "core/memspan.hpp"
#include "core/types.hpp"

/**
 * @brief defines a raw pointer that can be allocated and free'ed
 * 
 */
struct MemPtr {
#ifndef NDEBUG
    size_t size         = 0;
    bool   is_allocated = false;
#endif
    real_t* ptr = nullptr;

    explicit MemPtr() noexcept {};
    ~MemPtr();

    void Allocate(const size_t element_size);
    void Free();
};

class MemData {
   private:
    real_t* const data_;
    size_t const  stride_[2];

   public:
    explicit MemData(const MemLayout& layout, const MemPtr& ptr) noexcept
        : stride_{layout.stride[0], layout.stride[1]},
          data_(ptr.ptr +
                (layout.gs[0] + layout.stride[0] *
                                    (layout.gs[1] + layout.stride[1] * layout.gs[1]))) {
        m_assert(m_isaligned(data_), "the value of data must be aligned!");
    };

    // get now the fancy read, write, etc
    __restrict real_t* RWrite(const lda_t ida = 0) const {
        return data_;
    }
    __restrict const real_t* RRead(const lda_t ida = 0) const {
        return data_;
    }

    real_t* Write(const lda_t ida = 0) const {
        return data_;
    }

    const real_t* Read(const lda_t ida = 0) const {
        return data_;
    }

    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Accessor(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0) const noexcept {
        return i0 + stride_[0] * (i1 + stride_[1] * (i2 + stride_[1] * ida));
    }
};

#endif