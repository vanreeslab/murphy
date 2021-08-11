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
    const size_t  stride_[2];

    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Boffset_(const lda_t ida) const noexcept { return stride_[0] * stride_[1] * stride_[1] * ida; };

   public:
    explicit MemData() = delete;
    explicit MemData(const MemPtr& ptr, const MemLayout& layout) noexcept
        : stride_{layout.stride[0], layout.stride[1]},
          data_(ptr.ptr + layout.shift +
                (layout.gs + layout.stride[0] * (layout.gs + layout.stride[1] * layout.gs))) {
        m_assert(m_isaligned(data_), "the value of data must be aligned!");
        m_assert(ptr.size >= layout.n_elem, "the size of the pointer = %ld must be >= the layout nelem = %ld", ptr.size, layout.n_elem);
    };

    __restrict real_t* RWrite(const lda_t ida = 0) const noexcept {
        return data_ + Boffset_(ida);
    }
    __restrict const real_t* RRead(const lda_t ida = 0) const noexcept {
        return data_ + Boffset_(ida);
    }
    real_t* Write(const lda_t ida = 0) const noexcept {
        return data_ + Boffset_(ida);
    }
    const real_t* Read(const lda_t ida = 0) const noexcept {
        return data_ + Boffset_(ida);
    }

    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Accessor(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        return i0 + stride_[0] * (i1 + stride_[1] * i2);
    }
};

class ConstMemData {
   private:
    const real_t* const data_;
    const size_t        stride_[2];
    
    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Boffset_(const lda_t ida) const noexcept { return stride_[0] * stride_[1] * stride_[1] * ida; };

   public:
    explicit ConstMemData() = delete;
    explicit ConstMemData(const MemPtr& ptr, const MemLayout& layout) noexcept
        : stride_{layout.stride[0], layout.stride[1]},
          data_(ptr.ptr + layout.shift +
                (layout.gs + layout.stride[0] * (layout.gs + layout.stride[1] * layout.gs))) {
        m_assert(m_isaligned(data_), "the value of data must be aligned!");
        m_assert(ptr.size >= layout.n_elem, "the size of the pointer = %ld must be >= the layout nelem = %ld", ptr.size, layout.n_elem);
    };

    __restrict const real_t* RRead(const lda_t ida = 0) const noexcept {
        return data_ + Boffset_(ida);
    }
    const real_t* Read(const lda_t ida = 0) const noexcept {
        return data_ + Boffset_(ida);
    }

    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Accessor(const bidx_t i0, const bidx_t i1, const bidx_t i2) const noexcept {
        return i0 + stride_[0] * (i1 + stride_[1] * i2);
    }
};

#endif