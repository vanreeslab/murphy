#ifndef SRC_CORE_DATA_HPP
#define SRC_CORE_DATA_HPP

#include "core/types.hpp"
#include "core/macros.hpp"

#include <functional>

template <typename T>
class m_Ptr {
    T* ptr_ = std::nullptr;
#ifndef NDEBUG
    size_t size_         = 0;
    bool   is_allocated_ = false;
#endif

    explicit m_Ptr() noexcept {};
    ~m_Ptr() noexcept {
#ifndef NDEBUG
        m_assert(!is_allocated_, "you must free the memory first!");
#endif
    }

    [[nodiscard]] T* ptr() noexcept const {
        return ptr_;
    }
#ifndef NDEBUG
    [[nodiscard]] size_t size() noexcept const {
        return size_;
    }
#endif

    void Allocate(const size_t size) {
        ptr_ = reinterpret_cast<T*>(m_calloc(size * sizeof(T)));
#ifndef NDEBUG
        size_         = size;
        is_allocated_ = true;
#endif
    }

    void Free() {
    }
};

struct m_Layout {
    // we could store the number of dimensions, but I would keep the gs size to 3
    // this allow us to not loop over the dimensions to compute the offset (see here-under)
    size_t stride = 0;
    size_t gs[3]  = {0, 0, 0};
};

// defines a function that takes 3 arguments and returns an offset
using accessor_t = std::function<bidx_t(const bidx_t, const bidx_t, const bidx_t, const bidx_t)>;

template <typename T>
class m_Data {
   private:
    T*     data_;
    size_t stride_mem_;
    size_t stride_data_;

   public:
    explicit m_Data(const m_Layout& layout, const m_Ptr& ptr) noexcept {
        // some definitions
        m_assert((M_ALIGNMENT % sizeof(T)) == 0, "the alignement must be a mutiple of %d", sizeof(T));
        const size_t chunk = M_ALIGNMENT / sizeof(T);  //alignement in terms of T
        const size_t mod   = (layout.stride % chunk);

        // get the strides
        stride_data_ = layou.stride;
        stride_mem_  = layout.stride - mod + (mod > 0) * chunk;

        // @PB if we want to check that the size of the pointer is wide enough then we need to know how many dimensions
        // is there in the layout. Why not, it can be nice for debugging purposes

        // get the sifted pointer
        const size_t offset = layout.gs[0] + stride_mem_ * (layout.gs[1] + stride_data_ * layout.gs[2]);
        m_assert((offset % chunk) == 0, "the offset must be a modulo of the chunk!");
        data_ = layout.ptr + offset;
    }

    // @PB the user should now be able to do:
    // accessor_t idx = data.Accessor()
    // ...
    // data[id(i0,i1,i2)] = 4.59;
    [[nodiscard]] __attribute__((always_inline)) inline bidx_t Accessor(const bidx_t i0, const bidx_t i1, const bidx_t i2, const bidx_t ida = 0) const noexcept {
        return i0 + stride_mem_ * (i1 + stride_data_ * (i2 + stride_data_ * ida));
    }

    // get now the fancy read, write, etc
    __restrict T* RWrite(const lda_t ida = 0) const {
        return data_;
    }
    __restrict T* RRead(const lda_t ida = 0) const {
        return data_;
    }

    // etc...
};

#endif